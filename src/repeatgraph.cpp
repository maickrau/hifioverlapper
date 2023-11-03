#include <mutex>
#include <iostream>
#include <vector>
#include <string>
#include "CommonUtils.h"
#include "ReadStorage.h"
#include "MatchIndex.h"

const uint64_t firstBitUint64_t = 1ull << 63ull;
const uint64_t maskUint64_t = firstBitUint64_t-1;

size_t popcount(uint64_t x);

class MatchGroup
{
public:
	class Match
	{
	public:
		uint32_t leftStart;
		uint32_t rightStart;
		uint32_t length;
	};
	size_t leftRead;
	size_t rightRead;
	bool rightFw;
	size_t leftStart;
	size_t leftEnd;
	size_t rightStart;
	size_t rightEnd;
	std::vector<Match> matches;
};

void writeGraph(std::string outputFileName, const std::vector<size_t>& nodeCoverage, const std::vector<size_t>& nodeLength, const phmap::flat_hash_map<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>, size_t>& edgeCoverage, const size_t minCoverage, const size_t k)
{
	std::ofstream graph { outputFileName };
	size_t countNodes = 0;
	size_t countEdges = 0;
	for (size_t i = 0; i < nodeCoverage.size(); i++)
	{
		if (nodeCoverage[i] < minCoverage) continue;
		graph << "S\tnode_" << i << "\t*\tLN:i:" << (nodeLength[i]+k-1) << "\tll:f:" << nodeCoverage[i] << "\tFC:i:" << nodeCoverage[i] * (nodeLength[i]+k-1) << std::endl;
		countNodes += 1;
	}
	for (auto pair : edgeCoverage)
	{
		if (pair.second < minCoverage) continue;
		graph << "L\tnode_" << pair.first.first.first << "\t" << (pair.first.first.second ? "+" : "-") << "\tnode_" << pair.first.second.first << "\t" << (pair.first.second.second ? "+" : "-") << "\t" << (k-1) << "M\tec:i:" << pair.second << std::endl;
		countEdges += 1;
	}
	std::cerr << countNodes << " graph nodes" << std::endl;
	std::cerr << countEdges << " graph edges" << std::endl;
}

uint64_t find(std::vector<uint64_t>& result, const size_t pos)
{
	assert((pos & firstBitUint64_t) == 0);
	if (result[pos] == (pos & maskUint64_t))
	{
		assert(result[pos] == pos);
		return pos;
	}
	uint64_t foundPos = result[pos];
	while (result[foundPos & maskUint64_t] != (foundPos & maskUint64_t))
	{
		uint64_t nextFoundPos = result[foundPos & maskUint64_t];
		if (foundPos & firstBitUint64_t) nextFoundPos ^= firstBitUint64_t;
		foundPos = nextFoundPos;
	}
	assert(result[foundPos & maskUint64_t] == (foundPos & maskUint64_t));
	uint64_t finalPos = foundPos;
	foundPos = result[pos];
	result[pos] = finalPos;
	while (result[foundPos & maskUint64_t] != (foundPos & maskUint64_t))
	{
		uint64_t nextFoundPos = result[foundPos & maskUint64_t];
		if (foundPos & firstBitUint64_t) finalPos ^= firstBitUint64_t;
		result[foundPos & maskUint64_t] = finalPos;
		foundPos = nextFoundPos;
	}
	return result[pos];
}

void mergeSegments(std::vector<uint64_t>& result, const size_t leftPos, const size_t rightPos, const bool fw)
{
	auto leftParent = find(result, leftPos);
	auto rightParent = find(result, rightPos);
	assert((result[leftParent & maskUint64_t] & maskUint64_t) == (leftParent & maskUint64_t));
	assert((result[rightParent & maskUint64_t] & maskUint64_t) == (rightParent & maskUint64_t));
	if ((leftParent & maskUint64_t) == (rightParent & maskUint64_t))
	{
		assert(((leftParent & firstBitUint64_t) == 0) ^ ((rightParent & firstBitUint64_t) == 0) ^ fw);
		return;
	}
	result[rightParent & maskUint64_t] = (leftParent & maskUint64_t) + ((((leftParent & firstBitUint64_t) == 0) ^ ((rightParent & firstBitUint64_t) == 0) ^ fw) ? 0 : firstBitUint64_t);
}

void mergeSegments(const std::vector<size_t>& readLengths, const std::vector<size_t>& countBeforeRead, std::vector<uint64_t>& result, const std::vector<RankBitvector>& breakpoints, const size_t leftRead, const bool leftFw, const size_t leftStart, const size_t leftEnd, const size_t rightRead, const bool rightFw, const size_t rightStart, const size_t rightEnd)
{
	assert(leftFw);
	assert(breakpoints[leftRead].get(leftStart));
	assert(breakpoints[leftRead].get(leftEnd));
	size_t leftFirst = breakpoints[leftRead].getRank(leftStart);
	size_t leftLast = breakpoints[leftRead].getRank(leftEnd);
	assert(leftLast > leftFirst);
	size_t rightFirst;
	size_t rightLast;
	if (rightFw)
	{
		assert(breakpoints[rightRead].get(rightStart));
		assert(breakpoints[rightRead].get(rightEnd));
		rightFirst = breakpoints[rightRead].getRank(rightStart);
		rightLast = breakpoints[rightRead].getRank(rightEnd);
	}
	else
	{
		assert(breakpoints[rightRead].get(readLengths[rightRead]-rightStart));
		assert(breakpoints[rightRead].get(readLengths[rightRead]-rightEnd));
		rightLast = breakpoints[rightRead].getRank(readLengths[rightRead]-rightStart);
		rightFirst = breakpoints[rightRead].getRank(readLengths[rightRead]-rightEnd);
	}
	assert(rightLast > rightFirst);
	assert(rightLast-rightFirst == leftLast-leftFirst);
	for (size_t i = 0; i < leftLast-leftFirst; i++)
	{
		if (rightFw)
		{
			mergeSegments(result, countBeforeRead[leftRead] + leftFirst + i, countBeforeRead[rightRead] + rightFirst+i, true);
		}
		else
		{
			mergeSegments(result, countBeforeRead[leftRead] + leftFirst + i, countBeforeRead[rightRead] + rightLast-i-1, false);
		}
	}
}

std::vector<uint64_t> mergeSegments(const std::vector<size_t>& readLengths, const std::vector<MatchGroup>& matches, const std::vector<RankBitvector>& breakpoints, const size_t countBreakpoints)
{
	assert(breakpoints.size() == readLengths.size());
	std::vector<uint64_t> result;
	assert(countBreakpoints >= 2*readLengths.size());
	result.resize(countBreakpoints - readLengths.size());
	std::vector<size_t> countBeforeRead;
	countBeforeRead.resize(readLengths.size(), 0);
	size_t countBeforeNow = 0;
	for (size_t i = 0; i < breakpoints.size(); i++)
	{
		countBeforeRead[i] = countBeforeNow;
		size_t index = 0;
		for (size_t j = 1; j < breakpoints[i].size(); j++)
		{
			if (!breakpoints[i].get(j)) continue;
			result[countBeforeNow + index] = countBeforeNow + index;
			index += 1;
		}
		countBeforeNow += index;
	}
	assert(countBeforeNow == result.size());
	for (size_t groupi = 0; groupi < matches.size(); groupi++)
	{
		for (size_t posi = 0; posi < matches[groupi].matches.size(); posi++)
		{
			mergeSegments(readLengths, countBeforeRead,result, breakpoints, matches[groupi].leftRead, true, matches[groupi].matches[posi].leftStart, matches[groupi].matches[posi].leftStart+matches[groupi].matches[posi].length, matches[groupi].rightRead, matches[groupi].rightFw, matches[groupi].matches[posi].rightStart, matches[groupi].matches[posi].rightStart+matches[groupi].matches[posi].length);
		}
	}
	for (size_t i = 0; i < result.size(); i++)
	{
		find(result, i);
	}
	return result;
}

phmap::flat_hash_map<uint64_t, size_t> getSegmentToNode(const std::vector<uint64_t>& segments)
{
	size_t nextIndex = 0;
	phmap::flat_hash_map<uint64_t, size_t> result;
	for (size_t i = 0; i < segments.size(); i++)
	{
		assert((segments[i] & maskUint64_t) < segments.size());
		if ((segments[i] & maskUint64_t) != i) continue;
		assert(segments[i] == i);
		assert(result.count(i) == 0);
		result[i] = nextIndex;
		nextIndex += 1;
	}
	return result;
}

phmap::flat_hash_map<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>, size_t> getEdgeCoverages(const std::vector<size_t>& readLengths, const phmap::flat_hash_map<uint64_t, size_t>& segmentToNode, const std::vector<uint64_t>& segments, const std::vector<RankBitvector>& breakpoints)
{
	phmap::flat_hash_map<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>, size_t> edgeCoverage;
	size_t i = 0;
	for (size_t readi = 0; readi < breakpoints.size(); readi++)
	{
		i += 1;
		for (size_t readpos = 1; readpos+1 < breakpoints[readi].size(); readpos++)
		{
			if (!breakpoints[readi].get(readpos)) continue;
			std::pair<size_t, bool> edgeFrom;
			edgeFrom.first = segmentToNode.at(segments[i-1] & maskUint64_t);
			edgeFrom.second = (segments[i-1] & firstBitUint64_t) == 0;
			std::pair<size_t, bool> edgeTo;
			edgeTo.first = segmentToNode.at(segments[i] & maskUint64_t);
			edgeTo.second = (segments[i] & firstBitUint64_t) == 0;
			auto key = canon(edgeFrom, edgeTo);
			edgeCoverage[key] += 1;
			i += 1;
		}
	}
	assert(i == segments.size());
	return edgeCoverage;
}

std::pair<bool, bool> extendBreakpointsFwFw(const std::vector<size_t>& readLengths, std::vector<RankBitvector>& breakpoints, size_t leftRead, size_t leftStart, size_t leftEnd, size_t rightRead, size_t rightStart, size_t rightEnd)
{
	assert(leftRead != rightRead);
	assert(leftEnd - leftStart == rightEnd - rightStart);
	std::vector<uint64_t>& leftBits = breakpoints[leftRead].getBits();
	std::vector<uint64_t>& rightBits = breakpoints[rightRead].getBits();
	size_t pos = 0;
	bool addedAnyLeft = false;
	bool addedAnyRight = false;
	while (pos <= leftEnd - leftStart)
	{
		size_t leftIndex = (leftStart + pos) / 64;
		size_t leftOffset = (leftStart + pos) % 64;
		size_t rightIndex = (rightStart + pos) / 64;
		size_t rightOffset = (rightStart + pos) % 64;
		size_t availableLeft = 64 - leftOffset;
		size_t availableRight = 64 - rightOffset;
		size_t take = std::min(availableLeft, availableRight);
		take = std::min(take, (leftEnd-leftStart) - pos + 1);
		assert(take <= 64);
		assert(take >= 1);
		uint64_t leftGot = leftBits[leftIndex] >> leftOffset;
		uint64_t rightGot = rightBits[rightIndex] >> rightOffset;
		uint64_t mask = 0xFFFFFFFFFFFFFFFFull;
		if (take < 64) mask = (1ull << (take)) - 1;
		assert(popcount(mask) == take);
		assert((mask & 1) == (1));
		leftGot &= mask;
		rightGot &= mask;
		if ((leftGot ^ rightGot) == 0)
		{
			pos += take;
			continue;
		}
		if ((leftGot & rightGot) != rightGot)
		{
			addedAnyLeft = true;
		}
		if ((leftGot & rightGot) != leftGot)
		{
			addedAnyRight = true;
		}
		uint64_t addThese = leftGot | rightGot;
		leftBits[leftIndex] |= (addThese << leftOffset);
		rightBits[rightIndex] |= (addThese << rightOffset);
		pos += take;
	}
	assert(pos == leftEnd - leftStart + 1);
	return std::make_pair(addedAnyLeft, addedAnyRight);
}

std::pair<bool, bool> extendBreakpointsFwBw(const std::vector<size_t>& readLengths, std::vector<RankBitvector>& breakpoints, size_t leftRead, size_t leftStart, size_t leftEnd, size_t rightRead, size_t rightStart, size_t rightEnd)
{
	assert(leftRead != rightRead);
	assert(leftEnd - leftStart == rightEnd - rightStart);
	bool addedAnyLeft = false;
	bool addedAnyRight = false;
	size_t leftIndex = leftStart;
	size_t rightIndex = readLengths[rightRead] - rightStart;
	for (size_t i = 0; i <= leftEnd-leftStart; i++)
	{
		bool leftbit = breakpoints[leftRead].get(leftIndex);
		bool rightbit = breakpoints[rightRead].get(rightIndex);
		if (leftbit && !rightbit)
		{
			breakpoints[rightRead].set(rightIndex, true);
			addedAnyRight = true;
		}
		else if (!leftbit && rightbit)
		{
			breakpoints[leftRead].set(leftIndex, true);
			addedAnyLeft = true;
		}
		leftIndex += 1;
		rightIndex -= 1;
	}
	return std::make_pair(addedAnyLeft, addedAnyRight);
}

std::pair<size_t, size_t> getMatchSpan(const std::vector<size_t>& readLengths, const MatchGroup& group, const MatchGroup::Match& pos, size_t readI)
{
	if (group.leftRead == readI)
	{
		assert(group.rightRead != readI);
		return std::make_pair(pos.leftStart, pos.leftStart+pos.length);
	}
	else
	{
		assert(group.rightRead == readI);
		if (group.rightFw)
		{
			return std::make_pair(pos.rightStart, pos.rightStart+pos.length);
		}
		else
		{
			return std::make_pair(readLengths[readI] - (pos.rightStart+pos.length), readLengths[readI] - pos.rightStart);
		}
	}
}

std::pair<size_t, size_t> getLinearTo2dIndices(const size_t linearKey, const RankBitvector& kmermatchToGroup, const std::vector<size_t>& numKmerMatchesBeforeGroupStart)
{
	size_t groupi = kmermatchToGroup.getRank(linearKey);
	if (kmermatchToGroup.get(linearKey)) groupi += 1;
	assert(groupi >= 1);
	groupi -= 1;
	assert(linearKey >= numKmerMatchesBeforeGroupStart[groupi]);
	assert(groupi == numKmerMatchesBeforeGroupStart.size()-1 || numKmerMatchesBeforeGroupStart[groupi+1] > linearKey);
	size_t posi = linearKey - numKmerMatchesBeforeGroupStart[groupi];
	return std::make_pair(groupi, posi);
}

void buildIntervalTree(const std::vector<size_t>& readLengths, const RankBitvector& kmermatchToGroup, const std::vector<size_t>& numKmerMatchesBeforeGroupStart, const std::vector<MatchGroup>& matches, std::vector<size_t>& relevantMatches, std::vector<uint32_t>& relevantMatchTreeMaxRight, std::vector<uint32_t>& relevantMatchTreeMinLeft, const size_t readi)
{
	std::vector<size_t> sortedMatches = relevantMatches;
	std::sort(sortedMatches.begin(), sortedMatches.end(), [&matches, &readLengths, readi, &kmermatchToGroup, &numKmerMatchesBeforeGroupStart](size_t left, size_t right)
	{
		size_t leftGroup, leftPos;
		size_t rightGroup, rightPos;
		std::tie(leftGroup, leftPos) = getLinearTo2dIndices(left & maskUint64_t, kmermatchToGroup, numKmerMatchesBeforeGroupStart);
		std::tie(rightGroup, rightPos) = getLinearTo2dIndices(right & maskUint64_t, kmermatchToGroup, numKmerMatchesBeforeGroupStart);
		return getMatchSpan(readLengths, matches[leftGroup], matches[leftGroup].matches[leftPos], readi).first < getMatchSpan(readLengths, matches[rightGroup], matches[rightGroup].matches[rightPos], readi).first;
	});
	size_t test = sortedMatches.size();
	{
		std::vector<size_t> tmp;
		std::swap(relevantMatches, tmp);
	}
	relevantMatches.resize(sortedMatches.size(), std::numeric_limits<size_t>::max());
	size_t treeMaxDepth = 0;
	while (test != 0)
	{
		treeMaxDepth += 1;
		test >>= 1;
	}
	size_t treeIndex = (1 << (treeMaxDepth-1))-1; // bottom-left-most index
	// in-order traversal
	for (size_t i = 0; i < sortedMatches.size(); i++)
	{
		assert(relevantMatches[treeIndex] == std::numeric_limits<size_t>::max());
		relevantMatches[treeIndex] = sortedMatches[i];
		if (i == sortedMatches.size()-1) break;
		// has right-child
		if (treeIndex*2+2 < sortedMatches.size())
		{
			treeIndex = treeIndex*2+2; // go one right-down
			while (treeIndex*2+1 < sortedMatches.size()) // go as much left-down as possible
			{
				treeIndex = treeIndex*2+1;
			}
		}
		else
		{
			// go as much up-left as possible
			while (treeIndex % 2 == 0)
			{
				treeIndex = (treeIndex - 1) / 2;
			}
			assert(treeIndex % 2 == 1);
			treeIndex = (treeIndex - 1) / 2; // one up-right
		}
	}
	relevantMatchTreeMaxRight.resize(sortedMatches.size(), 0);
	relevantMatchTreeMinLeft.resize(sortedMatches.size(), 0);
	for (size_t i = sortedMatches.size()-1; i < sortedMatches.size(); i--)
	{
		size_t thisgroup, thispos;
		std::tie(thisgroup, thispos) = getLinearTo2dIndices(relevantMatches[i] & maskUint64_t, kmermatchToGroup, numKmerMatchesBeforeGroupStart);
		relevantMatchTreeMaxRight[i] = getMatchSpan(readLengths, matches[thisgroup], matches[thisgroup].matches[thispos], readi).second;
		relevantMatchTreeMinLeft[i] = getMatchSpan(readLengths, matches[thisgroup], matches[thisgroup].matches[thispos], readi).first;
		if (i * 2 + 1 < sortedMatches.size())
		{
			size_t childgroup, childpos;
			std::tie(childgroup, childpos) = getLinearTo2dIndices(relevantMatches[i*2+1] & maskUint64_t, kmermatchToGroup, numKmerMatchesBeforeGroupStart);
			relevantMatchTreeMaxRight[i] = std::max(relevantMatchTreeMaxRight[i], relevantMatchTreeMaxRight[i*2+1]);
			relevantMatchTreeMinLeft[i] = std::min(relevantMatchTreeMinLeft[i], relevantMatchTreeMinLeft[i*2+1]);
			assert(getMatchSpan(readLengths, matches[childgroup], matches[childgroup].matches[childpos], readi).first <= getMatchSpan(readLengths, matches[childgroup], matches[childgroup].matches[childpos], readi).first);
		}
		if (i * 2 + 2 < sortedMatches.size())
		{
			size_t childgroup, childpos;
			std::tie(childgroup, childpos) = getLinearTo2dIndices(relevantMatches[i*2+2] & maskUint64_t, kmermatchToGroup, numKmerMatchesBeforeGroupStart);
			relevantMatchTreeMaxRight[i] = std::max(relevantMatchTreeMaxRight[i], relevantMatchTreeMaxRight[i*2+2]);
			relevantMatchTreeMinLeft[i] = std::min(relevantMatchTreeMinLeft[i], relevantMatchTreeMinLeft[i*2+2]);
			assert(getMatchSpan(readLengths, matches[childgroup], matches[childgroup].matches[childpos], readi).first <= getMatchSpan(readLengths, matches[childgroup], matches[childgroup].matches[childpos], readi).first);
		}
	}
}

template <typename F>
void iterateIntervalTreeAlnMatches(const std::vector<size_t>& readLengths, const std::vector<size_t>& relevantMatches, const std::vector<MatchGroup>& matches, const RankBitvector& kmermatchToGroup, const std::vector<size_t>& numKmerMatchesBeforeGroupStart, const std::vector<uint32_t>& relevantMatchTreeMaxRight, const std::vector<uint32_t>& relevantMatchTreeMinLeft, const size_t readi, const size_t intervalStart, const size_t intervalEnd, F callback)
{
	if (relevantMatchTreeMaxRight[0] < intervalStart) return;
	if (relevantMatchTreeMinLeft[0] > intervalEnd) return;
	// weird stack structure so no function call recursion, weird pos / stack insert to minimize stack pushes
	std::vector<size_t> checkStack;
	size_t pos = 0;
	while (true)
	{
		size_t aln = relevantMatches[pos];
		size_t rawAln = aln & maskUint64_t;
		bool isRight = (aln & firstBitUint64_t) == firstBitUint64_t;
		size_t groupi, posi;
		std::tie(groupi, posi) = getLinearTo2dIndices(rawAln, kmermatchToGroup, numKmerMatchesBeforeGroupStart);
		std::pair<size_t, size_t> interval;
		if (isRight)
		{
			interval = getMatchSpan(readLengths, matches[groupi], matches[groupi].matches[posi], matches[groupi].rightRead);
		}
		else
		{
			interval = getMatchSpan(readLengths, matches[groupi], matches[groupi].matches[posi], matches[groupi].leftRead);
		}
		if (interval.first <= intervalEnd && interval.second >= intervalStart) callback(aln);
		if (relevantMatchTreeMinLeft[pos] > intervalEnd)
		{
			if (checkStack.size() == 0) return;
			pos = checkStack.back();
			checkStack.pop_back();
			continue;
		}
		if (relevantMatchTreeMaxRight[pos] < intervalStart)
		{
			if (checkStack.size() == 0) return;
			pos = checkStack.back();
			checkStack.pop_back();
			continue;
		}
		bool descendRight = true;
		bool descendLeft = true;
		if (pos * 2 + 1 >= relevantMatches.size()) descendLeft = false;
		if (pos * 2 + 2 >= relevantMatches.size()) descendRight = false;
		if (descendLeft && relevantMatchTreeMaxRight[pos*2+1] < intervalStart) descendLeft = false;
		if (descendLeft && relevantMatchTreeMinLeft[pos*2+1] > intervalEnd) descendLeft = false;
		if (descendRight && relevantMatchTreeMaxRight[pos*2+2] < intervalStart) descendRight = false;
		if (descendRight && relevantMatchTreeMinLeft[pos*2+2] > intervalEnd) descendRight = false;
		if (descendRight && descendLeft)
		{
			checkStack.push_back(pos*2+2);
			pos = pos*2+1;
		}
		else if (descendRight && !descendLeft)
		{
			pos = pos*2+2;
		}
		else if (descendLeft && !descendRight)
		{
			pos = pos*2+1;
		}
		else
		{
			assert(!descendLeft);
			assert(!descendRight);
			if (checkStack.size() == 0) return;
			pos = checkStack.back();
			checkStack.pop_back();
		}
	}
}

std::vector<RankBitvector> extendBreakpoints(const std::vector<size_t>& readLengths, const std::vector<MatchGroup>& matches)
{
	std::vector<RankBitvector> breakpoints;
	breakpoints.resize(readLengths.size());
	for (size_t i = 0; i < readLengths.size(); i++)
	{
		breakpoints[i].resize(readLengths[i]+1);
		breakpoints[i].set(0, true);
		breakpoints[i].set(readLengths[i], true);
	}
	std::vector<std::vector<uint64_t>> relevantMatches;
	relevantMatches.resize(readLengths.size());
	size_t firstFwBwMatch = std::numeric_limits<size_t>::max();
	RankBitvector kmermatchToGroup;
	std::vector<size_t> numKmerMatchesBeforeGroupStart;
	numKmerMatchesBeforeGroupStart.resize(matches.size(), 0);
	size_t totalMatchCount = 0;
	for (size_t groupi = 0; groupi < matches.size(); groupi++)
	{
		totalMatchCount += matches[groupi].matches.size();
	}
	kmermatchToGroup.resize(totalMatchCount);
	assert(matches.size() <= (size_t)std::numeric_limits<uint32_t>::max()/2);
	size_t kmerNumber = 0;
	for (size_t groupi = 0; groupi < matches.size(); groupi++)
	{
		numKmerMatchesBeforeGroupStart[groupi] = kmerNumber;
		kmermatchToGroup.set(kmerNumber, true);
		if (matches[groupi].rightFw)
		{
			assert(firstFwBwMatch == std::numeric_limits<size_t>::max());
		}
		else
		{
			if (firstFwBwMatch == std::numeric_limits<size_t>::max())
			{
				firstFwBwMatch = kmerNumber;
			}
			assert(kmerNumber >= firstFwBwMatch);
		}
		for (size_t posi = 0; posi < matches[groupi].matches.size(); posi++)
		{
			assert(matches[groupi].matches.size() <= (size_t)std::numeric_limits<uint32_t>::max());
			size_t key = kmerNumber + posi;
			relevantMatches[matches[groupi].leftRead].emplace_back(key);
			relevantMatches[matches[groupi].rightRead].emplace_back(key + firstBitUint64_t);
			breakpoints[matches[groupi].leftRead].set(matches[groupi].matches[posi].leftStart, true);
			breakpoints[matches[groupi].leftRead].set(matches[groupi].matches[posi].leftStart + matches[groupi].matches[posi].length, true);
			if (matches[groupi].rightFw)
			{
				breakpoints[matches[groupi].rightRead].set(matches[groupi].matches[posi].rightStart, true);
				breakpoints[matches[groupi].rightRead].set(matches[groupi].matches[posi].rightStart + matches[groupi].matches[posi].length, true);
			}
			else
			{
				breakpoints[matches[groupi].rightRead].set(readLengths[matches[groupi].rightRead] - matches[groupi].matches[posi].rightStart, true);
				breakpoints[matches[groupi].rightRead].set(readLengths[matches[groupi].rightRead] - (matches[groupi].matches[posi].rightStart + matches[groupi].matches[posi].length), true);
			}
		}
		kmerNumber += matches[groupi].matches.size();
	}
	assert(kmerNumber == totalMatchCount);
	kmermatchToGroup.buildRanks();
	std::vector<std::vector<uint32_t>> relevantMatchTreeMaxRight;
	std::vector<std::vector<uint32_t>> relevantMatchTreeMinLeft;
	relevantMatchTreeMaxRight.resize(relevantMatches.size());
	relevantMatchTreeMinLeft.resize(relevantMatches.size());
	for (size_t i = 0; i < relevantMatches.size(); i++)
	{
		buildIntervalTree(readLengths, kmermatchToGroup, numKmerMatchesBeforeGroupStart, matches, relevantMatches[i], relevantMatchTreeMaxRight[i], relevantMatchTreeMinLeft[i], i);
	}
	std::vector<bool> inQueue;
	inQueue.resize(totalMatchCount, true);
	std::vector<size_t> fwfwCheckQueue;
	std::vector<size_t> fwbwCheckQueue;
	size_t nextUnqueued = 0;
	while (fwfwCheckQueue.size() >= 1 || fwbwCheckQueue.size() >= 1 || nextUnqueued < totalMatchCount)
	{
		size_t matchkey;
		size_t groupi, posi;
		std::pair<bool, bool> addedAny;
		if (fwfwCheckQueue.size() >= 1)
		{
			matchkey = fwfwCheckQueue.back();
			fwfwCheckQueue.pop_back();
		}
		else if (nextUnqueued < firstFwBwMatch && nextUnqueued < totalMatchCount)
		{
			matchkey = nextUnqueued;
			nextUnqueued += 1;
		}
		else if (fwbwCheckQueue.size() >= 1)
		{
			matchkey = fwbwCheckQueue.back();
			fwbwCheckQueue.pop_back();
		}
		else
		{
			assert(nextUnqueued >= firstFwBwMatch);
			assert(nextUnqueued < totalMatchCount);
			matchkey = nextUnqueued;
			nextUnqueued += 1;
		}
		inQueue[matchkey] = false;
		std::tie(groupi, posi) = getLinearTo2dIndices(matchkey, kmermatchToGroup, numKmerMatchesBeforeGroupStart);
		if (matches[groupi].rightFw)
		{
			addedAny = extendBreakpointsFwFw(readLengths, breakpoints, matches[groupi].leftRead, matches[groupi].matches[posi].leftStart, matches[groupi].matches[posi].leftStart+matches[groupi].matches[posi].length, matches[groupi].rightRead, matches[groupi].matches[posi].rightStart, matches[groupi].matches[posi].rightStart+matches[groupi].matches[posi].length);
		}
		else
		{
			addedAny = extendBreakpointsFwBw(readLengths, breakpoints, matches[groupi].leftRead, matches[groupi].matches[posi].leftStart, matches[groupi].matches[posi].leftStart+matches[groupi].matches[posi].length, matches[groupi].rightRead, matches[groupi].matches[posi].rightStart, matches[groupi].matches[posi].rightStart+matches[groupi].matches[posi].length);
		}
		if (addedAny.first)
		{
			auto span = getMatchSpan(readLengths, matches[groupi], matches[groupi].matches[posi], matches[groupi].leftRead);
			iterateIntervalTreeAlnMatches(readLengths, relevantMatches[matches[groupi].leftRead], matches, kmermatchToGroup, numKmerMatchesBeforeGroupStart, relevantMatchTreeMaxRight[matches[groupi].leftRead], relevantMatchTreeMinLeft[matches[groupi].leftRead], matches[groupi].leftRead, span.first, span.second, [firstFwBwMatch, &fwfwCheckQueue, &fwbwCheckQueue, &inQueue, matchkey, &kmermatchToGroup, &numKmerMatchesBeforeGroupStart](size_t aln)
			{
				size_t rawAln = aln & maskUint64_t;
				if (inQueue[rawAln]) return;
				if (rawAln == matchkey) return;
				if (rawAln < firstFwBwMatch)
				{
					fwfwCheckQueue.push_back(rawAln);
				}
				else
				{
					fwbwCheckQueue.push_back(rawAln);
				}
				inQueue[rawAln] = true;
			});
		}
		if (addedAny.second)
		{
			auto span = getMatchSpan(readLengths, matches[groupi], matches[groupi].matches[posi], matches[groupi].rightRead);
			iterateIntervalTreeAlnMatches(readLengths, relevantMatches[matches[groupi].rightRead], matches, kmermatchToGroup, numKmerMatchesBeforeGroupStart, relevantMatchTreeMaxRight[matches[groupi].rightRead], relevantMatchTreeMinLeft[matches[groupi].rightRead], matches[groupi].rightRead, span.first, span.second, [firstFwBwMatch, &fwfwCheckQueue, &fwbwCheckQueue, &inQueue, matchkey, &kmermatchToGroup, &numKmerMatchesBeforeGroupStart](size_t aln)
			{
				size_t rawAln = aln & maskUint64_t;
				if (inQueue[rawAln]) return;
				if (rawAln == matchkey) return;
				if (rawAln < firstFwBwMatch)
				{
					fwfwCheckQueue.push_back(rawAln);
				}
				else
				{
					fwbwCheckQueue.push_back(rawAln);
				}
				inQueue[rawAln] = true;
			});
		}
	}
	assert(nextUnqueued == totalMatchCount);
	for (size_t i = 0; i < breakpoints.size(); i++) breakpoints[i].buildRanks();
	return breakpoints;
}

std::vector<size_t> getNodeCoverage(const std::vector<uint64_t>& segments, const phmap::flat_hash_map<uint64_t, size_t>& segmentToNode)
{
	std::vector<size_t> result;
	result.resize(segmentToNode.size(), 0);
	for (size_t i = 0; i < segments.size(); i++)
	{
		uint64_t node = segmentToNode.at(segments[i] & maskUint64_t);
		assert(node < result.size());
		result[node] += 1;
	}
	for (size_t i = 0; i < result.size(); i++)
	{
		assert(result[i] != 0);
	}
	return result;
}

std::vector<size_t> getNodeLengths(const std::vector<uint64_t>& segments, const phmap::flat_hash_map<uint64_t, size_t>& segmentToNode, const std::vector<RankBitvector>& breakpoints)
{
	std::vector<size_t> result;
	result.resize(segmentToNode.size(), std::numeric_limits<size_t>::max());
	size_t i = 0;
	for (size_t readi = 0; readi < breakpoints.size(); readi++)
	{
		size_t lastBreakpoint = 0;
		assert(breakpoints[readi].get(0));
		for (size_t pos = 1; pos < breakpoints[readi].size(); pos++)
		{
			if (!breakpoints[readi].get(pos)) continue;
			const size_t node = segmentToNode.at(segments[i] & maskUint64_t);
			const size_t length = pos - lastBreakpoint;
			assert(node < result.size());
			assert(length >= 1);
			if (!(result[node] == std::numeric_limits<size_t>::max() || result[node] == length))
			{
				std::cerr << length << " " << result[node] << std::endl;
				std::cerr << i << " " << segments[i] << " " << node << std::endl;
			}
			assert(result[node] == std::numeric_limits<size_t>::max() || result[node] == length);
			result[node] = length;
			i += 1;
			lastBreakpoint = pos;
		}
	}
	assert(i == segments.size());
	for (size_t i = 0; i < result.size(); i++)
	{
		assert(result[i] != std::numeric_limits<size_t>::max());
	}
	return result;
}

void makeGraph(const std::vector<size_t>& readLengths, const std::vector<MatchGroup>& matches, const size_t minCoverage, const std::string& outputFileName, const size_t k)
{
	std::vector<RankBitvector> breakpoints = extendBreakpoints(readLengths, matches);
	size_t countBreakpoints = 0;
	for (size_t i = 0; i < breakpoints.size(); i++)
	{
		for (size_t j = 0; j < breakpoints[i].size(); j++)
		{
			if (breakpoints[i].get(j)) countBreakpoints += 1;
		}
	}
	std::cerr << countBreakpoints << " breakpoints" << std::endl;
	std::vector<uint64_t> segments = mergeSegments(readLengths, matches, breakpoints, countBreakpoints);
	phmap::flat_hash_map<uint64_t, size_t> segmentToNode = getSegmentToNode(segments);
	std::vector<size_t> nodeCoverage = getNodeCoverage(segments, segmentToNode);
	std::vector<size_t> nodeLength = getNodeLengths(segments, segmentToNode, breakpoints);
	phmap::flat_hash_map<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>, size_t> edgeCoverages = getEdgeCoverages(readLengths, segmentToNode, segments, breakpoints);
	writeGraph(outputFileName, nodeCoverage, nodeLength, edgeCoverages, minCoverage, k);
}

template <typename F>
void iterateKmerMatchPositions(const uint64_t kmer, const phmap::flat_hash_map<uint64_t, uint32_t>& firstPositions, const phmap::flat_hash_map<uint64_t, std::vector<uint32_t>>& extraPositions, F callback)
{
	auto found = firstPositions.find(kmer);
	if (found == firstPositions.end()) return;
	callback(found->second);
	auto found2 = extraPositions.find(kmer);
	if (found2 == extraPositions.end()) return;
	for (auto pos : found2->second) callback(pos);
}

template <typename F1, typename F2>
void iterateSyncmers(const std::vector<TwobitString>& readSequences, const size_t k, const size_t w, const size_t maxLen, F1 callback, F2 getchar)
{
	static std::vector<std::tuple<size_t, uint64_t>> smerOrder;
	assert(w < k);
	assert(w >= 3);
	const uint64_t mask = (1ull << (2ull*k)) - 1;
	const size_t s = k-w+1;
	const uint64_t smask = (1ull << (2ull*s)) - 1;
	uint64_t kmer = 0;
	uint64_t smer = 0;
	for (size_t i = 0; i < s; i++)
	{
		uint64_t c = getchar(i);
		smer <<= 2;
		smer += c;
	}
	kmer = smer;
	smerOrder.emplace_back(0, smer);
	for (size_t i = s; i < k; i++)
	{
		uint64_t c = getchar(i);
		kmer <<= 2;
		kmer += c;
		smer <<= 2;
		smer += c;
		smer &= smask;
		while (smerOrder.size() > 0 && std::get<1>(smerOrder.back()) > smer) smerOrder.pop_back();
		smerOrder.emplace_back(i-s+1, smer);
	}
	if ((std::get<0>(smerOrder.front()) == 0) || (std::get<1>(smerOrder.back()) == std::get<1>(smerOrder.front()) && std::get<0>(smerOrder.back()) == w-1))
	{
		callback(kmer, 0);
	}
	for (size_t i = k; i < maxLen; i++)
	{
		uint64_t c = getchar(i);
		kmer <<= 2;
		kmer += c;
		kmer &= mask;
		smer <<= 2;
		smer += c;
		smer &= smask;
		while (smerOrder.size() > 0 && std::get<1>(smerOrder.back()) > smer) smerOrder.pop_back();
		while (smerOrder.size() > 0 && std::get<0>(smerOrder.front()) <= i-s+1-w) smerOrder.erase(smerOrder.begin());
		smerOrder.emplace_back(i-s+1, smer);
		if ((std::get<0>(smerOrder.front()) == i-s+2-w) || (std::get<1>(smerOrder.back()) == std::get<1>(smerOrder.front()) && std::get<0>(smerOrder.back()) == i-s+1))
		{
			callback(kmer, i-k+1);
		}
	}
	smerOrder.clear();
}

template <typename F>
void iterateSyncmers(const std::vector<TwobitString>& readSequences, const size_t k, const size_t w, const size_t read, const size_t readStart, const size_t readEnd, const bool fw, F callback)
{
	if (fw)
	{
		iterateSyncmers(readSequences, k, w, readEnd-readStart, callback, [&readSequences, read, readStart, readEnd](size_t index){ return readSequences[read].get(readStart+index); });
	}
	else
	{
		iterateSyncmers(readSequences, k, w, readEnd-readStart, callback, [&readSequences, read, readStart, readEnd](size_t index){ return 3-readSequences[read].get(readSequences[read].size() - 1 - (readStart+index)); });
	}
}

void getKmerMatches(const std::vector<TwobitString>& readSequences, MatchGroup& mappingMatch, const size_t k, const size_t w)
{
	assert(k <= 31);
	assert(k % 2 == 1);
	size_t leftStart = mappingMatch.leftStart;
	size_t leftEnd = mappingMatch.leftEnd;
	size_t rightStart = mappingMatch.rightStart;
	size_t rightEnd = mappingMatch.rightEnd;
	bool rightFw = mappingMatch.rightFw;
	size_t left = mappingMatch.leftRead;
	size_t right = mappingMatch.rightRead;
	phmap::flat_hash_map<uint64_t, uint32_t> firstKmerPositionInLeft;
	phmap::flat_hash_map<uint64_t, std::vector<uint32_t>> extraKmerPositionsInLeft;
	iterateSyncmers(readSequences, k, 20, left, leftStart, leftEnd, true, [&firstKmerPositionInLeft, &extraKmerPositionsInLeft](const size_t kmer, const size_t pos)
	{
		if (firstKmerPositionInLeft.count(kmer) == 0)
		{
			firstKmerPositionInLeft[kmer] = pos;
		}
		else
		{
			extraKmerPositionsInLeft[kmer].push_back(pos);
		}
	});
	std::vector<std::pair<size_t, size_t>> currentMatchesPerDiagonal;
	size_t diagonalCount = 2*w+1;
	size_t zeroDiagonal = w;
	if (leftEnd-leftStart > rightEnd-rightStart)
	{
		diagonalCount += (leftEnd-leftStart)-(rightEnd-rightStart);
		zeroDiagonal = w + (leftEnd-leftStart)-(rightEnd-rightStart);
	}
	if (rightEnd-rightStart > leftEnd-leftStart) diagonalCount += (rightEnd-rightStart)-(leftEnd-leftStart);
	currentMatchesPerDiagonal.resize(diagonalCount, std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max()));
	iterateSyncmers(readSequences, k, 20, right, rightStart, rightEnd, rightFw, [&firstKmerPositionInLeft, &extraKmerPositionsInLeft, &currentMatchesPerDiagonal, diagonalCount, zeroDiagonal, rightFw, left, right, leftStart, rightStart, leftEnd, rightEnd, w, k, &mappingMatch](const size_t kmer, const size_t rightPos)
	{
		size_t interpolatedLeftPos = (double)(rightPos) / (double)(rightEnd-rightStart) * (double)(leftEnd-leftStart);
		assert(rightPos + zeroDiagonal >= interpolatedLeftPos + w);
		assert(rightPos + zeroDiagonal + w >= interpolatedLeftPos);
		size_t minDiagonal = rightPos + zeroDiagonal - interpolatedLeftPos - w;
		size_t maxDiagonal = rightPos + zeroDiagonal - interpolatedLeftPos + w;
		assert(maxDiagonal <= diagonalCount);
		iterateKmerMatchPositions(kmer, firstKmerPositionInLeft, extraKmerPositionsInLeft, [zeroDiagonal, rightPos, minDiagonal, maxDiagonal, &currentMatchesPerDiagonal, &mappingMatch, rightFw, left, right, leftStart, rightStart, diagonalCount, k](const size_t leftPos)
		{
			if (leftPos > zeroDiagonal + rightPos) return;
			if (zeroDiagonal + rightPos - leftPos >= diagonalCount) return;
			size_t diagonal = zeroDiagonal + rightPos - leftPos;
			if (diagonal < minDiagonal || diagonal > maxDiagonal) return;
			if (currentMatchesPerDiagonal[diagonal].first != std::numeric_limits<size_t>::max() && currentMatchesPerDiagonal[diagonal].second + k > rightPos)
			{
				currentMatchesPerDiagonal[diagonal].second = rightPos+1;
				return;
			}
			if (currentMatchesPerDiagonal[diagonal].first != std::numeric_limits<size_t>::max())
			{
				assert(currentMatchesPerDiagonal[diagonal].second != std::numeric_limits<size_t>::max());
				assert(currentMatchesPerDiagonal[diagonal].second > currentMatchesPerDiagonal[diagonal].first);
				assert(currentMatchesPerDiagonal[diagonal].first + zeroDiagonal >= diagonal);
				size_t leftMatchStart = leftStart + currentMatchesPerDiagonal[diagonal].first + zeroDiagonal - diagonal;
				size_t rightMatchStart = rightStart + currentMatchesPerDiagonal[diagonal].first;
				size_t length = currentMatchesPerDiagonal[diagonal].second - currentMatchesPerDiagonal[diagonal].first;
				mappingMatch.matches.emplace_back();
				mappingMatch.matches.back().leftStart = leftMatchStart;
				mappingMatch.matches.back().rightStart = rightMatchStart;
				mappingMatch.matches.back().length = length;
			}
			currentMatchesPerDiagonal[diagonal].first = rightPos;
			currentMatchesPerDiagonal[diagonal].second = rightPos+1;
		});
	});
	for (size_t diagonal = 0; diagonal < diagonalCount; diagonal++)
	{
		if (currentMatchesPerDiagonal[diagonal].first == std::numeric_limits<size_t>::max()) continue;
		assert(currentMatchesPerDiagonal[diagonal].second != std::numeric_limits<size_t>::max());
		assert(currentMatchesPerDiagonal[diagonal].second > currentMatchesPerDiagonal[diagonal].first);
		assert(currentMatchesPerDiagonal[diagonal].first + zeroDiagonal >= diagonal);
		size_t leftMatchStart = leftStart + currentMatchesPerDiagonal[diagonal].first + zeroDiagonal - diagonal;
		size_t rightMatchStart = rightStart + currentMatchesPerDiagonal[diagonal].first;
		size_t length = currentMatchesPerDiagonal[diagonal].second - currentMatchesPerDiagonal[diagonal].first;
		mappingMatch.matches.emplace_back();
		mappingMatch.matches.back().leftStart = leftMatchStart;
		mappingMatch.matches.back().rightStart = rightMatchStart;
		mappingMatch.matches.back().length = length;
	}
}

int main(int argc, char** argv)
{
	size_t numThreads = std::stoi(argv[1]);
	size_t k = std::stoull(argv[2]);
	size_t numWindows = std::stoull(argv[3]);
	size_t windowSize = std::stoull(argv[4]);
	size_t minAlignmentLength = std::stoull(argv[5]);
	const size_t graphk = 31;
	const size_t minCoverage = 2;
	std::vector<std::string> readFiles;
	for (size_t i = 6; i < argc; i++)
	{
		readFiles.emplace_back(argv[i]);
	}
	MatchIndex matchIndex { k, numWindows, windowSize };
	ReadStorage storage;
	std::vector<TwobitString> readSequences;
	for (auto file : readFiles)
	{
		std::mutex indexMutex;
		storage.iterateReadsFromFile(file, numThreads, false, [&matchIndex, &indexMutex, &readSequences](size_t readName, const std::string& sequence)
		{
			matchIndex.addMatchesFromRead(readName, indexMutex, sequence);
			readSequences.emplace_back(sequence);
		});
	}
	size_t numWindowChunks = matchIndex.numWindowChunks();
	size_t numUniqueChunks = matchIndex.numUniqueChunks();
	matchIndex.clearConstructionVariablesAndCompact();
	const std::vector<size_t>& readBasepairLengths = storage.getRawReadLengths();
	std::vector<size_t> readKmerLengths;
	readKmerLengths.resize(readBasepairLengths.size());
	for (size_t i = 0; i < readBasepairLengths.size(); i++)
	{
		readKmerLengths[i] = readBasepairLengths[i] + 1 - graphk;
	}
	const std::vector<std::string>& readNames = storage.getNames();
	std::mutex printMutex;
	std::vector<MatchGroup> matches;
	auto result = matchIndex.iterateMatchChains(numThreads, storage.getRawReadLengths(), [&printMutex, &matches, minAlignmentLength](const size_t left, const size_t leftstart, const size_t leftend, const bool leftFw, const size_t right, const size_t rightstart, const size_t rightend, const bool rightFw)
	{
		if (leftend+1-leftstart < minAlignmentLength) return;
		if (rightend+1-rightstart < minAlignmentLength) return;
		std::lock_guard<std::mutex> lock { printMutex };
		assert(leftFw);
		matches.emplace_back();
		matches.back().leftRead = left;
		matches.back().rightRead = right;
		matches.back().rightFw = rightFw;
		matches.back().leftStart = leftstart;
		matches.back().leftEnd = leftend+1;
		matches.back().rightStart = rightstart;
		matches.back().rightEnd = rightend+1;
	});
	// fw-fw matches first, fw-bw matches later
	std::stable_sort(matches.begin(), matches.end(), [](const auto& left, const auto& right){ return left.rightFw && !right.rightFw; });
	std::cerr << matches.size() << " mapping matches" << std::endl;
	size_t kmerMatchCount = 0;
	for (auto& t : matches)
	{
		getKmerMatches(readSequences, t, graphk, 50);
		kmerMatchCount += t.matches.size();
		for (auto pos : t.matches)
		{
			assert(pos.leftStart >= t.leftStart);
			assert(pos.rightStart >= t.rightStart);
			assert(pos.leftStart + pos.length <= t.leftEnd);
			assert(pos.rightStart + pos.length <= t.rightEnd);
			// std::cout << readNames[std::get<0>(match)] << "\t" << readKmerLengths[std::get<0>(match)] << "\t" << std::get<1>(match) << "\t" << std::get<2>(match) << "\t" << readNames[std::get<3>(match)] << "\t" << readKmerLengths[std::get<3>(match)] << "\t" << std::get<4>(match) << "\t" << std::get<5>(match) << "\t" << (std::get<6>(match) ? "fw" : "bw") << std::endl;
		}
	}
	std::cerr << kmerMatchCount << " kmer matches" << std::endl;
	makeGraph(readKmerLengths, matches, minCoverage, "graph.gfa", graphk);
	std::cerr << result.numberReads << " reads" << std::endl;
	std::cerr << numWindowChunks << " distinct windowchunks" << std::endl;
	std::cerr << result.totalReadChunkMatches << " read-windowchunk matches (except unique)" << std::endl;
	std::cerr << numUniqueChunks << " windowchunks have only one read" << std::endl;
	std::cerr << result.readsWithMatch << " reads with a match" << std::endl;
	std::cerr << result.readPairMatches << " read-read matches" << std::endl;
	std::cerr << result.readChainMatches << " chain matches" << std::endl;
	std::cerr << result.totalMatches << " window matches" << std::endl;
	std::cerr << result.maxPerChunk << " max windowchunk size" << std::endl;
}
