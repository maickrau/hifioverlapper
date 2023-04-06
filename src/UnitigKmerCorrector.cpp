#include <queue>
#include <cassert>
#include <mutex>
#include "UnitigKmerCorrector.h"

UnitigKmerCorrector::UnitigKmerCorrector(size_t k) :
	kmerSize(k),
	unitigs(k)
{
}

void UnitigKmerCorrector::build(const ReadpartIterator& partIterator)
{
	std::mutex mutex;
	std::cerr << "collect k-mers" << std::endl;
	partIterator.iterateHashes([this, &mutex](const ReadInfo& read, const SequenceCharType& seq, const SequenceLengthType& poses, const std::string& rawSeq, const std::vector<size_t>& positions, const std::vector<HashType>& hashes)
	{
		std::lock_guard lock { mutex };
		std::pair<size_t, bool> last { std::numeric_limits<size_t>::max(), true };
		for (size_t i = 0; i < hashes.size(); i++)
		{
			std::pair<size_t, bool> thisKmer = unitigs.getNode(hashes[i]);
			if (last.first != std::numeric_limits<size_t>::max())
			{
				unitigs.addEdge(last, thisKmer, kmerSize - (positions[i] - positions[i-1]));
			}
			last = thisKmer;
		}
	});
	std::cerr << unitigs.numHashes() << " hashes" << std::endl;
	std::cerr << "build unitigs" << std::endl;
	unitigs.buildUnitigGraph();
	std::cerr << unitigs.numUnitigs() << " unitigs" << std::endl;
	std::cerr << "build unitig sequences and collect read paths" << std::endl;
	partIterator.iterateHashes([this, &mutex](const ReadInfo& read, const SequenceCharType& seq, const SequenceLengthType& poses, const std::string& rawSeq, const std::vector<size_t>& positions, const std::vector<HashType>& hashes)
	{
		std::lock_guard lock { mutex };
		for (size_t i = 0; i < hashes.size(); i++)
		{
			std::pair<size_t, bool> thisKmer = unitigs.getNodeOrNull(hashes[i]);
			unitigs.addKmerSequence(thisKmer, rawSeq, positions[i], positions[i]+kmerSize);
		}
		reads.emplace_back();
		reads.back().name = read.readName.first;
		std::tie(reads.back().leftClip, reads.back().rightClip, reads.back().unitigPath) = unitigs.getPath(hashes);
		if (reads.back().unitigPath.size() == 0)
		{
			reads.back().leftHanger = rawSeq;
		}
		else
		{
			reads.back().leftHanger = rawSeq.substr(0, positions[0]);
			reads.back().rightHanger = rawSeq.substr(positions.back() + kmerSize);
		}
	});
	std::cerr << reads.size() << " reads" << std::endl;
	std::cerr << "finalize" << std::endl;
	unitigs.finalizeSequences();
	std::cerr << unitigs.totalBps() << " bp" << std::endl;
}

class DistantKmerComparer
{
public:
	bool operator()(const std::pair<size_t, std::pair<size_t, bool>> left, const std::pair<size_t, std::pair<size_t, bool>> right) const
	{
		return left.first > right.first;
	}
};

std::vector<std::pair<size_t, bool>> UnitigKmerCorrector::getUniqueReplacementPath(const std::pair<size_t, bool> start, const std::pair<size_t, bool> end, const std::vector<bool>& allowedNodes, const VectorWithDirection<std::vector<std::pair<size_t, bool>>>& allowedEdges, const std::vector<size_t>& localToUnitig, size_t maxLength) const
{
	if (end == start) return std::vector<std::pair<size_t, bool>>{};
	assert(end.first < allowedEdges.size());
	assert(start.first < allowedEdges.size());
	phmap::flat_hash_set<std::pair<size_t, bool>> bwVisited;
	std::priority_queue<std::pair<size_t, std::pair<size_t, bool>>, std::vector<std::pair<size_t, std::pair<size_t, bool>>>, DistantKmerComparer> checkStack;
	checkStack.emplace(0, reverse(end));
	while (checkStack.size() > 0)
	{
		auto topandlen = checkStack.top();
		checkStack.pop();
		size_t distance = topandlen.first;
		std::pair<size_t, bool> top = topandlen.second;
		if (distance > maxLength) continue;
		if (bwVisited.count(top) == 1) continue;
		bwVisited.insert(top);
		if (bwVisited.size() > 1000) return std::vector<std::pair<size_t, bool>>{};
		for (auto edge : allowedEdges[top])
		{
			if (edge == start) return std::vector<std::pair<size_t, bool>>{};
			if (!allowedNodes[edge.first]) continue;
			std::pair<size_t, bool> fromglobal { localToUnitig[top.first], top.second };
			std::pair<size_t, bool> toglobal { localToUnitig[edge.first], edge.second };
			size_t extra = unitigs.unitigMinusEdgeLength(fromglobal, toglobal);
			checkStack.emplace(distance + extra, edge);
		}
	}
	std::vector<std::pair<size_t, bool>> result;
	result.push_back(start);
	phmap::flat_hash_set<std::pair<size_t, bool>> visited;
	while (result.back() != end)
	{
		std::pair<size_t, bool> uniqueOption { std::numeric_limits<size_t>::max(), false };
		if (allowedEdges[result.back()].size() == 0) return std::vector<std::pair<size_t, bool>>{};
		for (auto edge : allowedEdges[result.back()])
		{
			if (bwVisited.count(reverse(edge)) == 0) continue;
			if (uniqueOption.first == std::numeric_limits<size_t>::max())
			{
				uniqueOption = edge;
			}
			else
			{
				return std::vector<std::pair<size_t, bool>>{};
			}
		}
		if (uniqueOption.first == std::numeric_limits<size_t>::max()) return std::vector<std::pair<size_t, bool>>{};
		assert(visited.count(uniqueOption) == 0);
		result.push_back(uniqueOption);
		visited.insert(uniqueOption);
	}
	assert(result.size() > 1);
	assert(result[0] == start);
	assert(result.back() == end);
	return result;
}

std::pair<std::string, bool> UnitigKmerCorrector::getCorrectedSequence(size_t readIndex, const std::vector<size_t>& context, size_t minAmbiguousCoverage, size_t minSafeCoverage) const
{
	if (reads[readIndex].unitigPath.size() == 1) return std::make_pair(getRawSequence(readIndex), false);
	phmap::flat_hash_map<size_t, size_t> unitigToLocal;
	std::vector<size_t> localToUnitig;
	std::vector<size_t> localCoverage;
	phmap::flat_hash_map<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>, size_t> localEdgeCoverage;
	std::vector<bool> hasFwCoverage;
	std::vector<bool> hasBwCoverage;
	for (size_t read : context)
	{
		std::pair<size_t, bool> lastLocal { std::numeric_limits<size_t>::max(), false };
		for (auto node : reads[read].unitigPath)
		{
			if (unitigToLocal.count(node.first) == 0)
			{
				unitigToLocal[node.first] = localCoverage.size();
				localToUnitig.push_back(node.first);
				localCoverage.push_back(0);
				hasFwCoverage.push_back(false);
				hasBwCoverage.push_back(false);
			}
			std::pair<size_t, bool> thislocal { unitigToLocal.at(node.first), node.second };
			localCoverage[thislocal.first] += 1;
			if (thislocal.second)
			{
				hasFwCoverage[thislocal.first] = true;
			}
			else
			{
				hasBwCoverage[thislocal.first] = true;
			}
			if (lastLocal.first != std::numeric_limits<size_t>::max()) localEdgeCoverage[canon(lastLocal, thislocal)] += 1;
			lastLocal = thislocal;
		}
	}
	std::vector<bool> safeNode;
	std::vector<bool> ambiguousNode;
	VectorWithDirection<std::vector<std::pair<size_t, bool>>> safeEdges;
	VectorWithDirection<std::vector<std::pair<size_t, bool>>> ambiguousEdges;
	safeNode.resize(localCoverage.size(), false);
	ambiguousNode.resize(localCoverage.size(), false);
	safeEdges.resize(localCoverage.size());
	ambiguousEdges.resize(localCoverage.size());
	for (size_t i = 0; i < localCoverage.size(); i++)
	{
		if (!hasFwCoverage[i]) continue;
		if (!hasBwCoverage[i]) continue;
		if (localCoverage[i] >= minSafeCoverage) safeNode[i] = true;
		if (localCoverage[i] >= minAmbiguousCoverage) ambiguousNode[i] = true;
	}
	for (auto pair : localEdgeCoverage)
	{
		if (!ambiguousNode[pair.first.first.first]) continue;
		if (!ambiguousNode[pair.first.second.first]) continue;
		if (pair.second >= minSafeCoverage)
		{
			safeEdges[pair.first.first].emplace_back(pair.first.second);
			safeEdges[reverse(pair.first.second)].emplace_back(reverse(pair.first.first));
		}
		if (pair.second >= minAmbiguousCoverage)
		{
			ambiguousEdges[pair.first.first].emplace_back(pair.first.second);
			ambiguousEdges[reverse(pair.first.second)].emplace_back(reverse(pair.first.first));
		}
	}
	std::vector<std::pair<size_t, bool>> correctedLocalPath;
	std::pair<size_t, bool> lastLocal { std::numeric_limits<size_t>::max(), false };
	size_t lastMatch = std::numeric_limits<size_t>::max();
	for (size_t i = 0; i < reads[readIndex].unitigPath.size(); i++)
	{
		std::pair<size_t, bool> thisLocal { unitigToLocal.at(reads[readIndex].unitigPath[i].first), reads[readIndex].unitigPath[i].second };
		if (!safeNode[thisLocal.first]) continue;
		if (lastMatch == std::numeric_limits<size_t>::max())
		{
			correctedLocalPath.insert(correctedLocalPath.end(), reads[readIndex].unitigPath.begin(), reads[readIndex].unitigPath.begin()+i+1);
			lastMatch = i;
			lastLocal = thisLocal;
			continue;
		}
		size_t maxLength = 0;
		for (size_t j = lastMatch+1; j <= i; j++)
		{
			maxLength += unitigs.unitigMinusEdgeLength(reads[readIndex].unitigPath[j-1], reads[readIndex].unitigPath[j]);
		}
		std::vector<std::pair<size_t, bool>> uniqueReplacement = getUniqueReplacementPath(lastLocal, thisLocal, ambiguousNode, ambiguousEdges, localToUnitig, maxLength+500);
		if (uniqueReplacement.size() == 0)
		{
			correctedLocalPath.insert(correctedLocalPath.end(), reads[readIndex].unitigPath.begin()+lastMatch+1, reads[readIndex].unitigPath.begin()+i+1);
			lastMatch = i;
			lastLocal = thisLocal;
			continue;
		}
		bool allSafe = true;
		for (size_t j = 0; j < uniqueReplacement.size(); j++)
		{
			if (!safeNode[uniqueReplacement[j].first]) allSafe = false;
		}
		if (allSafe)
		{
			for (size_t j = 1; j < uniqueReplacement.size(); j++)
			{
				bool found = false;
				for (auto edge : safeEdges[uniqueReplacement[j-1]])
				{
					if (edge == uniqueReplacement[j])
					{
						found = true;
						break;
					}
				}
				if (!found)
				{
					allSafe = false;
					break;
				}
			}
		}
		if (!allSafe)
		{
			correctedLocalPath.insert(correctedLocalPath.end(), reads[readIndex].unitigPath.begin()+lastMatch+1, reads[readIndex].unitigPath.begin()+i+1);
			lastMatch = i;
			lastLocal = thisLocal;
			continue;
		}
		for (size_t j = 0; j < uniqueReplacement.size(); j++)
		{
			uniqueReplacement[j].first = localToUnitig[uniqueReplacement[j].first];
		}
		assert(uniqueReplacement[0] == reads[readIndex].unitigPath[lastMatch]);
		assert(uniqueReplacement.back() == reads[readIndex].unitigPath[i]);
		correctedLocalPath.insert(correctedLocalPath.end(), uniqueReplacement.begin()+1, uniqueReplacement.end());
		lastMatch = i;
		lastLocal = thisLocal;
	}
	correctedLocalPath.insert(correctedLocalPath.end(), reads[readIndex].unitigPath.begin()+lastMatch+1, reads[readIndex].unitigPath.end());
	bool changed = false;
	if (correctedLocalPath.size() != reads[readIndex].unitigPath.size())
	{
		changed = true;
	}
	else
	{
		for (size_t i = 0; i < correctedLocalPath.size(); i++)
		{
			if (correctedLocalPath[i] != reads[readIndex].unitigPath[i]) changed = true;
		}
	}
	assert(correctedLocalPath.size() >= 2);
	assert(correctedLocalPath[0] == reads[readIndex].unitigPath[0]);
	assert(correctedLocalPath.back() == reads[readIndex].unitigPath.back());
	std::string result = reads[readIndex].leftHanger + unitigs.getSequence(correctedLocalPath, reads[readIndex].leftClip, reads[readIndex].rightClip) + reads[readIndex].rightHanger;
	return std::make_pair(result, changed);
}

std::string UnitigKmerCorrector::getRawSequence(size_t index) const
{
	std::string result;
	result = reads[index].leftHanger + unitigs.getSequence(reads[index].unitigPath, reads[index].leftClip, reads[index].rightClip) + reads[index].rightHanger;
	return result;
}

size_t UnitigKmerCorrector::numReads() const
{
	return reads.size();
}

const std::string& UnitigKmerCorrector::getName(size_t index) const
{
	return reads[index].name;
}
