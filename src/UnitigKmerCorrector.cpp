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

std::vector<std::pair<size_t, bool>> UnitigKmerCorrector::getUniqueReplacementPath(const std::pair<size_t, bool> start, const std::pair<size_t, bool> end, const phmap::flat_hash_set<size_t>& allowedNodes, const phmap::flat_hash_map<std::pair<size_t, bool>, phmap::flat_hash_set<std::pair<size_t, bool>>>& allowedEdges) const
{
	if (end == start) return std::vector<std::pair<size_t, bool>>{};
	phmap::flat_hash_set<std::pair<size_t, bool>> bwVisited;
	std::priority_queue<std::pair<size_t, std::pair<size_t, bool>>, std::vector<std::pair<size_t, std::pair<size_t, bool>>>, DistantKmerComparer> checkStack;
	checkStack.emplace(0, reverse(end));
	while (checkStack.size() > 0)
	{
		auto topandlen = checkStack.top();
		checkStack.pop();
		size_t distance = topandlen.first;
		std::pair<size_t, bool> top = topandlen.second;
		// if (distance > maxLength) continue;
		if (bwVisited.count(top) == 1) continue;
		bwVisited.insert(top);
		if (bwVisited.size() > 1000) return std::vector<std::pair<size_t, bool>>{};
		if (allowedEdges.count(top) == 0) return std::vector<std::pair<size_t, bool>>{};
		for (auto edge : allowedEdges.at(top))
		{
			if (edge == start) return std::vector<std::pair<size_t, bool>>{};
			if (allowedNodes.count(edge.first) == 0) continue;
			// size_t extra = kmerSize - reads.getOverlap(top, edge);
			size_t extra = 1;
			checkStack.emplace(distance + extra, edge);
		}
	}
	std::vector<std::pair<size_t, bool>> result;
	result.push_back(start);
	phmap::flat_hash_set<std::pair<size_t, bool>> visited;
	while (result.back() != end)
	{
		std::pair<size_t, bool> uniqueOption { std::numeric_limits<size_t>::max(), false };
		if (allowedEdges.count(result.back()) == 0) return std::vector<std::pair<size_t, bool>>{};
		for (auto edge : allowedEdges.at(result.back()))
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
	phmap::flat_hash_map<size_t, size_t> unitigCoverage;
	phmap::flat_hash_map<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>, size_t> edgeCoverage;
	phmap::flat_hash_set<size_t> hasFwCoverage;
	phmap::flat_hash_set<size_t> hasBwCoverage;
	for (size_t read : context)
	{
		for (size_t i = 0; i < reads[read].unitigPath.size(); i++)
		{
			unitigCoverage[reads[read].unitigPath[i].first] += 1;
			if (reads[read].unitigPath[i].second)
			{
				hasFwCoverage.insert(reads[read].unitigPath[i].first);
			}
			else
			{
				hasBwCoverage.insert(reads[read].unitigPath[i].first);
			}
		}
		for (size_t i = 1; i < reads[read].unitigPath.size(); i++)
		{
			edgeCoverage[canon(reads[read].unitigPath[i-1], reads[read].unitigPath[i])] += 1;
		}
	}
	phmap::flat_hash_set<size_t> safeNode;
	phmap::flat_hash_set<size_t> ambiguousNode;
	phmap::flat_hash_map<std::pair<size_t, bool>, phmap::flat_hash_set<std::pair<size_t, bool>>> safeEdges;
	phmap::flat_hash_map<std::pair<size_t, bool>, phmap::flat_hash_set<std::pair<size_t, bool>>> ambiguousEdges;
	for (auto pair : unitigCoverage)
	{
		if (hasFwCoverage.count(pair.first) == 0) continue;
		if (hasBwCoverage.count(pair.first) == 0) continue;
		if (pair.second >= minSafeCoverage) safeNode.insert(pair.first);
		if (pair.second >= minAmbiguousCoverage) ambiguousNode.insert(pair.first);
	}
	for (auto pair : edgeCoverage)
	{
		if (hasFwCoverage.count(pair.first.first.first) == 0) continue;
		if (hasBwCoverage.count(pair.first.first.first) == 0) continue;
		if (hasFwCoverage.count(pair.first.second.first) == 0) continue;
		if (hasBwCoverage.count(pair.first.second.first) == 0) continue;
		if (pair.second >= minSafeCoverage)
		{
			safeEdges[pair.first.first].emplace(pair.first.second);
			safeEdges[reverse(pair.first.second)].emplace(reverse(pair.first.first));
		}
		if (pair.second >= minAmbiguousCoverage)
		{
			ambiguousEdges[pair.first.first].emplace(pair.first.second);
			ambiguousEdges[reverse(pair.first.second)].emplace(reverse(pair.first.first));
		}
	}
	std::vector<std::pair<size_t, bool>> correctedLocalPath;
	size_t lastMatch = std::numeric_limits<size_t>::max();
	for (size_t i = 0; i < reads[readIndex].unitigPath.size(); i++)
	{
		if (safeNode.count(reads[readIndex].unitigPath[i].first) == 0) continue;
		if (lastMatch == std::numeric_limits<size_t>::max())
		{
			correctedLocalPath.insert(correctedLocalPath.end(), reads[readIndex].unitigPath.begin(), reads[readIndex].unitigPath.begin()+i+1);
			lastMatch = i;
			continue;
		}
		std::vector<std::pair<size_t, bool>> uniqueReplacement = getUniqueReplacementPath(reads[readIndex].unitigPath[lastMatch], reads[readIndex].unitigPath[i], ambiguousNode, ambiguousEdges);
		if (uniqueReplacement.size() == 0)
		{
			correctedLocalPath.insert(correctedLocalPath.end(), reads[readIndex].unitigPath.begin()+lastMatch+1, reads[readIndex].unitigPath.begin()+i+1);
			lastMatch = i;
			continue;
		}
		bool allSafe = true;
		for (size_t i = 0; i < uniqueReplacement.size(); i++)
		{
			if (safeNode.count(uniqueReplacement[i].first) == 0) allSafe = false;
		}
		for (size_t i = 1; i < uniqueReplacement.size(); i++)
		{
			if (safeEdges[uniqueReplacement[i-1]].count(uniqueReplacement[i]) == 0) allSafe = false;
		}
		if (!allSafe)
		{
			correctedLocalPath.insert(correctedLocalPath.end(), reads[readIndex].unitigPath.begin()+lastMatch+1, reads[readIndex].unitigPath.begin()+i+1);
			lastMatch = i;
			continue;
		}
		assert(uniqueReplacement[0] == reads[readIndex].unitigPath[lastMatch]);
		assert(uniqueReplacement.back() == reads[readIndex].unitigPath[i]);
		correctedLocalPath.insert(correctedLocalPath.end(), uniqueReplacement.begin()+1, uniqueReplacement.end());
		lastMatch = i;
	}
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
	correctedLocalPath.insert(correctedLocalPath.end(), reads[readIndex].unitigPath.begin()+lastMatch+1, reads[readIndex].unitigPath.end());
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
