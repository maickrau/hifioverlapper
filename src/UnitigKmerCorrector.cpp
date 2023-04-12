#include <queue>
#include <cassert>
#include <mutex>
#include "UnitigKmerCorrector.h"

UnitigKmerCorrector::UnitigKmerCorrector(size_t k) :
	kmerSize(k),
	unitigs(k)
{
}

std::pair<size_t, bool> findSimpleBubbleEnd(const VectorWithDirection<std::vector<std::pair<size_t, bool>>>& edges, std::pair<size_t, bool> startNode)
{
	if (edges[startNode].size() < 2) return std::make_pair(std::numeric_limits<size_t>::max(), true);
	std::pair<size_t, bool> bubbleEnd { std::numeric_limits<size_t>::max(), true };
	phmap::flat_hash_set<size_t> usedAlleles;
	for (auto edge : edges[startNode])
	{
		std::pair<size_t, bool> pos = edge;
		while (edges[pos].size() == 1 && edges[reverse(pos)].size() == 1)
		{
			if (usedAlleles.count(pos.first) == 1) return std::make_pair(std::numeric_limits<size_t>::max(), true);
			usedAlleles.insert(pos.first);
			pos = edges[pos][0];
		}
		if (bubbleEnd.first == std::numeric_limits<size_t>::max())
		{
			bubbleEnd = pos;
		}
		else if (pos != bubbleEnd) return std::make_pair(std::numeric_limits<size_t>::max(), true);
	}
	if (edges[reverse(bubbleEnd)].size() != edges[startNode].size()) return std::make_pair(std::numeric_limits<size_t>::max(), true);
	if (bubbleEnd.first == startNode.first) return std::make_pair(std::numeric_limits<size_t>::max(), true);
	return bubbleEnd;
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

std::vector<std::pair<size_t, bool>> UnitigKmerCorrector::getUniqueReplacementPath(const std::pair<size_t, bool> start, const std::pair<size_t, bool> end, const std::vector<bool>& allowedNodes, const VectorWithDirection<std::vector<std::pair<size_t, bool>>>& allowedEdges, const std::vector<size_t>& localToGlobal, size_t maxLength) const
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
			std::pair<size_t, bool> fromglobal { localToGlobal[top.first], top.second };
			std::pair<size_t, bool> toglobal { localToGlobal[edge.first], edge.second };
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

UnitigKmerCorrector::LocalGraph UnitigKmerCorrector::getLocalGraph(const std::vector<size_t>& context, size_t minAmbiguousCoverage, size_t minSafeCoverage) const
{
	LocalGraph result;
	std::vector<size_t> localCoverage;
	phmap::flat_hash_map<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>, size_t> localEdgeCoverage;
	std::vector<bool> hasFwCoverage;
	std::vector<bool> hasBwCoverage;
	for (size_t read : context)
	{
		std::pair<size_t, bool> lastLocal { std::numeric_limits<size_t>::max(), false };
		for (auto node : reads[read].unitigPath)
		{
			if (result.globalToLocal.count(node.first) == 0)
			{
				result.globalToLocal[node.first] = localCoverage.size();
				result.localToGlobal.push_back(node.first);
				localCoverage.push_back(0);
				hasFwCoverage.push_back(false);
				hasBwCoverage.push_back(false);
			}
			std::pair<size_t, bool> thislocal { result.globalToLocal.at(node.first), node.second };
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
	result.safeNode.resize(localCoverage.size(), false);
	result.ambiguousNode.resize(localCoverage.size(), false);
	result.safeEdges.resize(localCoverage.size());
	result.ambiguousEdges.resize(localCoverage.size());
	for (size_t i = 0; i < localCoverage.size(); i++)
	{
		if (!hasFwCoverage[i]) continue;
		if (!hasBwCoverage[i]) continue;
		if (localCoverage[i] >= minSafeCoverage) result.safeNode[i] = true;
		if (localCoverage[i] >= minAmbiguousCoverage) result.ambiguousNode[i] = true;
	}
	for (auto pair : localEdgeCoverage)
	{
		if (!result.ambiguousNode[pair.first.first.first]) continue;
		if (!result.ambiguousNode[pair.first.second.first]) continue;
		if (pair.second >= minSafeCoverage)
		{
			result.safeEdges[pair.first.first].emplace_back(pair.first.second);
			result.safeEdges[reverse(pair.first.second)].emplace_back(reverse(pair.first.first));
		}
		if (pair.second >= minAmbiguousCoverage)
		{
			result.ambiguousEdges[pair.first.first].emplace_back(pair.first.second);
			result.ambiguousEdges[reverse(pair.first.second)].emplace_back(reverse(pair.first.first));
		}
	}
	return result;
}

std::pair<std::string, bool> UnitigKmerCorrector::getCorrectedSequence(size_t readIndex, const std::vector<size_t>& context, size_t minAmbiguousCoverage, size_t minSafeCoverage) const
{
	if (reads[readIndex].unitigPath.size() == 1) return std::make_pair(getRawSequence(readIndex), false);
	LocalGraph graph = getLocalGraph(context, minAmbiguousCoverage, minSafeCoverage);
	std::vector<std::pair<size_t, bool>> correctedLocalPath;
	std::pair<size_t, bool> lastLocal { std::numeric_limits<size_t>::max(), false };
	size_t lastMatch = std::numeric_limits<size_t>::max();
	for (size_t i = 0; i < reads[readIndex].unitigPath.size(); i++)
	{
		std::pair<size_t, bool> thisLocal { graph.globalToLocal.at(reads[readIndex].unitigPath[i].first), reads[readIndex].unitigPath[i].second };
		if (!graph.safeNode[thisLocal.first]) continue;
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
		std::vector<std::pair<size_t, bool>> uniqueReplacement = getUniqueReplacementPath(lastLocal, thisLocal, graph.ambiguousNode, graph.ambiguousEdges, graph.localToGlobal, maxLength+500);
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
			if (!graph.safeNode[uniqueReplacement[j].first]) allSafe = false;
		}
		if (allSafe)
		{
			for (size_t j = 1; j < uniqueReplacement.size(); j++)
			{
				bool found = false;
				for (auto edge : graph.safeEdges[uniqueReplacement[j-1]])
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
			uniqueReplacement[j].first = graph.localToGlobal[uniqueReplacement[j].first];
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

void UnitigKmerCorrector::assignReadsToAlleles(const std::vector<size_t>& context, const std::vector<size_t>& localToGlobal, std::vector<std::vector<size_t>>& result, std::vector<std::pair<size_t, bool>> alleles) const
{
	assert(alleles.size() >= 2);
	result.resize(alleles.size()+1);
	phmap::flat_hash_map<size_t, size_t> nodeToAllele;
	for (size_t i = 0; i < alleles.size(); i++)
	{
		assert(nodeToAllele.count(localToGlobal[alleles[i].first]) == 0);
		nodeToAllele[localToGlobal[alleles[i].first]] = i+1;
	}
	for (size_t read : context)
	{
		size_t readAlleleCount = 0;
		for (std::pair<size_t, bool> node : reads[read].unitigPath)
		{
			if (nodeToAllele.count(node.first) == 0) continue;
			result[nodeToAllele.at(node.first)].push_back(read);
			readAlleleCount += 1;
			if (readAlleleCount >= 2) break;
		}
		if (readAlleleCount >= 2)
		{
			for (std::pair<size_t, bool> node : reads[read].unitigPath)
			{
				if (nodeToAllele.count(node.first) == 0) continue;
				assert(result[nodeToAllele.at(node.first)].size() >= 1);
				assert(result[nodeToAllele.at(node.first)].back() == read);
				result[nodeToAllele.at(node.first)].pop_back();
				readAlleleCount -= 1;
				if (readAlleleCount == 0) break;
			}
			result[0].push_back(read);
		}
	}
	for (size_t i = 0; i < result.size(); i++)
	{
		std::sort(result[i].begin(), result[i].end());
	}
}

size_t countMatches(const std::vector<size_t>& sortedLeft, const std::vector<size_t>& sortedRight)
{
	size_t result = 0;
	size_t leftPos = 0;
	size_t rightPos = 0;
	while (leftPos < sortedLeft.size() && rightPos < sortedRight.size())
	{
		if (sortedLeft[leftPos] == sortedRight[rightPos])
		{
			result += 1;
			leftPos += 1;
			rightPos += 1;
			continue;
		}
		if (sortedLeft[leftPos] < sortedRight[rightPos])
		{
			leftPos += 1;
			continue;
		}
		if (sortedLeft[leftPos] > sortedRight[rightPos])
		{
			rightPos += 1;
			continue;
		}
	}
	return result;
}

bool checkPhasingValidity(const std::vector<std::vector<size_t>>& left, const std::vector<std::vector<size_t>>& right, size_t minAmbiguousCoverage, size_t minSafeCoverage)
{
	if (left.size() != right.size()) return false; // todo: implement variable number alleles. algorithm exists but need to write
	if (left[0].size() + right[0].size() >= minAmbiguousCoverage) return false;
	std::vector<size_t> uniqueRightMatch;
	uniqueRightMatch.resize(right.size(), 0);
	size_t nonMatchReads = 0;
	for (size_t i = 1; i < left.size(); i++)
	{
		size_t uniqueMatch = 0;
		for (size_t j = 1; j < right.size(); j++)
		{
			size_t matches = countMatches(left[i], right[j]);
			if (matches >= minSafeCoverage)
			{
				if (uniqueMatch != 0) return false;
				if (uniqueRightMatch[j] != 0) return false;
				uniqueMatch = j;
				uniqueRightMatch[j] = i;
			}
			else
			{
				nonMatchReads += matches;
			}
		}
	}
	if (nonMatchReads + left[0].size() + right[0].size() >= minAmbiguousCoverage) return false;
	return true;
}

void UnitigKmerCorrector::forbidOtherHaplotypes(phmap::flat_hash_set<size_t>& forbiddenReads, size_t readIndex, const std::vector<std::vector<size_t>>& leftAlleles, const std::vector<std::vector<size_t>>& rightAlleles) const
{
	//todo: implement variable number alleles. algorithm exists but need to write
	size_t leftIndex = 0;
	size_t rightIndex = 0;
	for (size_t i = 1; i < leftAlleles.size(); i++)
	{
		for (auto read : leftAlleles[i])
		{
			if (read == readIndex)
			{
				leftIndex = i;
				break;
			}
		}
		if (leftIndex != 0) break;
	}
	for (size_t i = 1; i < rightAlleles.size(); i++)
	{
		for (auto read : rightAlleles[i])
		{
			if (read == readIndex)
			{
				rightIndex = i;
				break;
			}
		}
		if (rightIndex != 0) break;
	}
	if (leftIndex == 0) return;
	if (rightIndex == 0) return;
	for (size_t i = 0; i < leftAlleles.size(); i++)
	{
		if (i == leftIndex) continue;
		for (auto read : leftAlleles[i])
		{
			assert(read != readIndex);
			forbiddenReads.insert(read);
		}
	}
	for (size_t i = 0; i < rightAlleles.size(); i++)
	{
		if (i == rightIndex) continue;
		for (auto read : rightAlleles[i])
		{
			assert(read != readIndex);
			forbiddenReads.insert(read);
		}
	}
}

std::vector<size_t> UnitigKmerCorrector::filterDifferentHaplotypesOut(size_t readIndex, const std::vector<size_t>& context, size_t minAmbiguousCoverage, size_t minSafeCoverage) const
{
	LocalGraph graph = getLocalGraph(context, minAmbiguousCoverage, minSafeCoverage);
	std::vector<std::vector<std::vector<size_t>>> bubbleAlleles;
	std::vector<uint8_t> coverageInKeyRead;
	coverageInKeyRead.resize(graph.size(), 0);
	for (auto node : reads[readIndex].unitigPath)
	{
		size_t localIndex = graph.globalToLocal.at(node.first);
		if (coverageInKeyRead[localIndex] < 2) coverageInKeyRead[localIndex] += 1;
	}
	for (size_t i = 0; i < graph.size(); i++)
	{
		if (coverageInKeyRead[i] != 1) continue;
		auto foundBubbleEnd = findSimpleBubbleEnd(graph.ambiguousEdges, std::make_pair(i, true));
		if (foundBubbleEnd.first != std::numeric_limits<size_t>::max() && foundBubbleEnd.first < i && coverageInKeyRead[foundBubbleEnd.first] == 1)
		{
			bubbleAlleles.emplace_back();
			assignReadsToAlleles(context, graph.localToGlobal, bubbleAlleles.back(), graph.ambiguousEdges[std::make_pair(i, true)]);
		}
		foundBubbleEnd = findSimpleBubbleEnd(graph.ambiguousEdges, std::make_pair(i, false));
		if (foundBubbleEnd.first != std::numeric_limits<size_t>::max() && foundBubbleEnd.first < i && coverageInKeyRead[foundBubbleEnd.first] == 1)
		{
			bubbleAlleles.emplace_back();
			assignReadsToAlleles(context, graph.localToGlobal, bubbleAlleles.back(), graph.ambiguousEdges[std::make_pair(i, false)]);
		}
	}
	phmap::flat_hash_set<size_t> forbiddenReads;
	for (size_t i = 0; i < bubbleAlleles.size(); i++)
	{
		if (bubbleAlleles[i][0].size() >= minAmbiguousCoverage) continue;
		for (size_t j = 0; j < i; j++)
		{
			if (bubbleAlleles[j][0].size() >= minAmbiguousCoverage) continue;
			bool bubbleValidlyPhased = checkPhasingValidity(bubbleAlleles[i], bubbleAlleles[j], minAmbiguousCoverage, minSafeCoverage);
			if (!bubbleValidlyPhased) continue;
			forbidOtherHaplotypes(forbiddenReads, readIndex, bubbleAlleles[i], bubbleAlleles[j]);
		}
	}
	assert(forbiddenReads.count(readIndex) == 0);
	std::vector<size_t> fixedContext;
	for (size_t read : context)
	{
		if (forbiddenReads.count(read) == 1) continue;
		fixedContext.push_back(read);
	}
	if (fixedContext.size() < context.size())
	{
		return filterDifferentHaplotypesOut(readIndex, fixedContext, minAmbiguousCoverage, minSafeCoverage);
	}
	return fixedContext;
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

size_t UnitigKmerCorrector::LocalGraph::size() const
{
	return safeNode.size();
}