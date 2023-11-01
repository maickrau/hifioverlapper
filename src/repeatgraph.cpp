#include <mutex>
#include <iostream>
#include <vector>
#include <string>
#include "CommonUtils.h"
#include "ReadStorage.h"
#include "MatchIndex.h"

size_t popcount(uint64_t x);

// add if missing
std::pair<size_t, bool> getBreakpoint(std::vector<std::tuple<size_t, size_t, size_t, bool>>& readBreakpoints, size_t readIndex, size_t wantedPos)
{
	for (auto pos : readBreakpoints)
	{
		if (std::get<0>(pos) != wantedPos) continue;
		return std::make_pair(std::get<0>(pos), false);
	}
	readBreakpoints.emplace_back(wantedPos, readIndex, wantedPos, true);
	return std::make_pair(wantedPos, true);
}

// for index finding and never add
size_t findBreakpointIndex(const std::vector<std::vector<std::tuple<size_t, size_t, size_t, bool>>>& breakpoints, size_t readIndex, size_t wantedPos)
{
	size_t found = std::numeric_limits<size_t>::max();
	for (size_t i = 0; i < breakpoints[readIndex].size(); i++)
	{
		if (std::get<0>(breakpoints[readIndex][i]) != wantedPos) continue;
		assert(found == std::numeric_limits<size_t>::max());
		found = i;
	}
	assert(found != std::numeric_limits<size_t>::max());
	return found;
}

std::tuple<size_t, size_t, bool> find(std::vector<std::vector<std::tuple<size_t, size_t, size_t, bool>>>& breakpoints, size_t readIndex, size_t wantedPos)
{
	size_t found = findBreakpointIndex(breakpoints, readIndex, wantedPos);
	if (std::get<1>(breakpoints[readIndex][found]) == readIndex && std::get<2>(breakpoints[readIndex][found]) == wantedPos)
	{
		assert(std::get<3>(breakpoints[readIndex][found]));
		return std::make_tuple(readIndex, wantedPos, true);
	}
	auto parent = find(breakpoints, std::get<1>(breakpoints[readIndex][found]), std::get<2>(breakpoints[readIndex][found]));
	if (!std::get<3>(breakpoints[readIndex][found])) std::get<2>(parent) = !std::get<2>(parent);
	std::get<1>(breakpoints[readIndex][found]) = std::get<0>(parent);
	std::get<2>(breakpoints[readIndex][found]) = std::get<1>(parent);
	std::get<3>(breakpoints[readIndex][found]) = std::get<2>(parent);
	return parent;
}

bool merge(std::vector<std::vector<std::tuple<size_t, size_t, size_t, bool>>>& breakpoints, size_t leftRead, size_t leftPos, size_t rightRead, size_t rightPos, bool forward)
{
	auto leftParent = find(breakpoints, leftRead, leftPos);
	auto rightParent = find(breakpoints, rightRead, rightPos);
	size_t leftIndex = findBreakpointIndex(breakpoints, std::get<0>(leftParent), std::get<1>(leftParent));
	size_t rightIndex = findBreakpointIndex(breakpoints, std::get<0>(rightParent), std::get<1>(rightParent));
	if (std::get<0>(leftParent) == std::get<0>(rightParent) && std::get<1>(leftParent) == std::get<1>(rightParent))
	{
		assert(std::get<2>(leftParent) ^ std::get<2>(rightParent) ^ forward);
		return false;
	}
	assert(std::get<1>(breakpoints[std::get<0>(leftParent)][leftIndex]) == std::get<0>(leftParent));
	assert(std::get<2>(breakpoints[std::get<0>(leftParent)][leftIndex]) == std::get<1>(leftParent));
	assert(std::get<1>(breakpoints[std::get<0>(rightParent)][rightIndex]) == std::get<0>(rightParent));
	assert(std::get<2>(breakpoints[std::get<0>(rightParent)][rightIndex]) == std::get<1>(rightParent));
	std::get<1>(breakpoints[std::get<0>(rightParent)][rightIndex]) = std::get<0>(leftParent);
	std::get<2>(breakpoints[std::get<0>(rightParent)][rightIndex]) = std::get<1>(leftParent);
	std::get<3>(breakpoints[std::get<0>(rightParent)][rightIndex]) = std::get<2>(leftParent) ^ std::get<2>(rightParent) ^ forward;
	return true;
}

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

std::tuple<uint32_t, uint32_t, bool> find(std::vector<std::vector<std::tuple<uint32_t, uint32_t, bool>>>& result, const size_t read, const size_t index)
{
	if (std::get<0>(result[read][index]) == read && std::get<1>(result[read][index]) == index)
	{
		assert(std::get<2>(result[read][index]));
		return result[read][index];
	}
	std::tuple<size_t, size_t, bool> foundPos = result[read][index];
	while (std::get<0>(result[std::get<0>(foundPos)][std::get<1>(foundPos)]) != std::get<0>(foundPos) || std::get<1>(result[std::get<0>(foundPos)][std::get<1>(foundPos)]) != std::get<1>(foundPos))
	{
		auto nextFoundPos = result[std::get<0>(foundPos)][std::get<1>(foundPos)];
		if (!std::get<2>(foundPos)) std::get<2>(nextFoundPos) = !std::get<2>(nextFoundPos);
		foundPos = nextFoundPos;
	}
	assert(std::get<0>(result[std::get<0>(foundPos)][std::get<1>(foundPos)]) == std::get<0>(foundPos));
	assert(std::get<1>(result[std::get<0>(foundPos)][std::get<1>(foundPos)]) == std::get<1>(foundPos));
	assert(std::get<2>(result[std::get<0>(foundPos)][std::get<1>(foundPos)]));
	std::tuple<size_t, size_t, bool> finalPos = foundPos;
	foundPos = result[read][index];
	result[read][index] = finalPos;
	while (std::get<0>(result[std::get<0>(foundPos)][std::get<1>(foundPos)]) != std::get<0>(foundPos) || std::get<1>(result[std::get<0>(foundPos)][std::get<1>(foundPos)]) != std::get<1>(foundPos))
	{
		auto nextFoundPos = result[std::get<0>(foundPos)][std::get<1>(foundPos)];
		if (!std::get<2>(foundPos))
		{
			std::get<2>(finalPos) = !std::get<2>(finalPos);
		}
		result[std::get<0>(foundPos)][std::get<1>(foundPos)] = finalPos;
		foundPos = nextFoundPos;
	}
	return result[read][index];
}

void mergeSegments(std::vector<std::vector<std::tuple<uint32_t, uint32_t, bool>>>& result, const size_t leftRead, const size_t leftIndex, const size_t rightRead, const size_t rightIndex, const bool fw)
{
	auto leftParent = find(result, leftRead, leftIndex);
	auto rightParent = find(result, rightRead, rightIndex);
	assert(std::get<0>(leftParent) == std::get<0>(result[std::get<0>(leftParent)][std::get<1>(leftParent)]));
	assert(std::get<1>(leftParent) == std::get<1>(result[std::get<0>(leftParent)][std::get<1>(leftParent)]));
	assert(std::get<0>(rightParent) == std::get<0>(result[std::get<0>(rightParent)][std::get<1>(rightParent)]));
	assert(std::get<1>(rightParent) == std::get<1>(result[std::get<0>(rightParent)][std::get<1>(rightParent)]));
	if (std::get<0>(leftParent) == std::get<0>(rightParent) && std::get<1>(leftParent) == std::get<1>(rightParent))
	{
		assert(std::get<2>(rightParent) ^ std::get<2>(leftParent) ^ fw);
		return;
	}
	result[std::get<0>(rightParent)][std::get<1>(rightParent)] = std::make_tuple(std::get<0>(leftParent), std::get<1>(leftParent), std::get<2>(leftParent) ^ std::get<2>(rightParent) ^ fw);
}

void mergeSegments(const std::vector<size_t>& readLengths, std::vector<std::vector<std::tuple<uint32_t, uint32_t, bool>>>& result, const std::vector<std::vector<uint32_t>>& segmentStarts, const std::vector<RankBitvector>& breakpoints, const size_t leftRead, const bool leftFw, const size_t leftStart, const size_t leftEnd, const size_t rightRead, const bool rightFw, const size_t rightStart, const size_t rightEnd)
{
	assert(leftFw);
	assert(breakpoints[leftRead].get(leftStart));
	assert(breakpoints[leftRead].get(leftEnd));
	size_t leftFirst = breakpoints[leftRead].getRank(leftStart);
	size_t leftLast = breakpoints[leftRead].getRank(leftEnd);
	assert(leftLast > leftFirst);
	assert(segmentStarts[leftRead][leftFirst] == leftStart);
	assert(segmentStarts[leftRead][leftLast] == leftEnd);
	size_t rightFirst;
	size_t rightLast;
	if (rightFw)
	{
		assert(breakpoints[rightRead].get(rightStart));
		assert(breakpoints[rightRead].get(rightEnd));
		rightFirst = breakpoints[rightRead].getRank(rightStart);
		rightLast = breakpoints[rightRead].getRank(rightEnd);
		assert(segmentStarts[rightRead][rightFirst] == rightStart);
		assert(segmentStarts[rightRead][rightLast] == rightEnd);
	}
	else
	{
		assert(breakpoints[rightRead].get(readLengths[rightRead]-rightStart));
		assert(breakpoints[rightRead].get(readLengths[rightRead]-rightEnd));
		rightLast = breakpoints[rightRead].getRank(readLengths[rightRead]-rightStart);
		rightFirst = breakpoints[rightRead].getRank(readLengths[rightRead]-rightEnd);
		assert(segmentStarts[rightRead][rightLast] == readLengths[rightRead]-rightStart);
		assert(segmentStarts[rightRead][rightFirst] == readLengths[rightRead]-rightEnd);
	}
	assert(rightLast > rightFirst);
	assert(rightLast-rightFirst == leftLast-leftFirst);
	for (size_t i = 0; i < leftLast-leftFirst; i++)
	{
		if (rightFw)
		{
			assert(segmentStarts[leftRead][leftFirst+i+1] - segmentStarts[leftRead][leftFirst+i] == segmentStarts[rightRead][rightFirst+i+1] - segmentStarts[rightRead][rightFirst+i]);
			mergeSegments(result, leftRead, leftFirst+i, rightRead, rightFirst+i, true);
		}
		else
		{
			assert(segmentStarts[leftRead][leftFirst+i+1] - segmentStarts[leftRead][leftFirst+i] == segmentStarts[rightRead][rightLast-i] - segmentStarts[rightRead][rightLast-i-1]);
			mergeSegments(result, leftRead, leftFirst+i, rightRead, rightLast-i-1, false);
		}
	}
}

std::tuple<std::vector<std::vector<std::tuple<uint32_t, uint32_t, bool>>>, std::vector<std::vector<uint32_t>>> mergeSegments(const std::vector<size_t>& readLengths, const std::vector<std::tuple<size_t, size_t, size_t, size_t, size_t, size_t, bool>>& matches, const std::vector<RankBitvector>& breakpoints)
{
	assert(breakpoints.size() == readLengths.size());
	std::vector<std::vector<std::tuple<uint32_t, uint32_t, bool>>> result;
	std::vector<std::vector<uint32_t>> segmentStarts;
	result.resize(breakpoints.size());
	segmentStarts.resize(breakpoints.size());
	for (size_t i = 0; i < breakpoints.size(); i++)
	{
		segmentStarts[i].push_back(0);
		assert(breakpoints[i].get(0));
		assert(breakpoints[i].get(readLengths[i]));
		for (size_t j = 1; j < breakpoints[i].size(); j++)
		{
			if (!breakpoints[i].get(j)) continue;
			result[i].emplace_back(i, result[i].size(), true);
			segmentStarts[i].push_back(j);
		}
		assert(result[i].size() >= 1);
		assert(segmentStarts[i].size() == result[i].size()+1);
		assert(result[i].size() == breakpoints[i].getRank(readLengths[i]));
	}
	for (auto match : matches)
	{
		mergeSegments(readLengths, result, segmentStarts, breakpoints, std::get<0>(match), true, std::get<1>(match), std::get<2>(match), std::get<3>(match), std::get<6>(match), std::get<4>(match), std::get<5>(match));
	}
	for (size_t i = 0; i < result.size(); i++)
	{
		for (size_t j = 0; j < result[i].size(); j++)
		{
			find(result, i, j);
		}
	}
	return std::make_pair(result, segmentStarts);
}

phmap::flat_hash_map<std::pair<uint32_t, uint32_t>, size_t> getSegmentToNode(const std::vector<std::vector<std::tuple<uint32_t, uint32_t, bool>>>& segments)
{
	size_t nextIndex = 0;
	phmap::flat_hash_map<std::pair<uint32_t, uint32_t>, size_t> result;
	for (size_t i = 0; i < segments.size(); i++)
	{
		for (size_t j = 0; j < segments[i].size(); j++)
		{
			if (std::get<0>(segments[i][j]) != i) continue;
			if (std::get<1>(segments[i][j]) != j) continue;
			assert(std::get<2>(segments[i][j]));
			result[std::make_pair(std::get<0>(segments[i][j]), std::get<1>(segments[i][j]))] = nextIndex;
			nextIndex += 1;
		}
	}
	return result;
}

phmap::flat_hash_map<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>, size_t> getEdgeCoverages(const std::vector<size_t>& readLengths, const phmap::flat_hash_map<std::pair<uint32_t, uint32_t>, size_t>& segmentToNode, const std::vector<std::vector<std::tuple<uint32_t, uint32_t, bool>>>& segments)
{
	phmap::flat_hash_map<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>, size_t> edgeCoverage;
	for (size_t i = 0; i < segments.size(); i++)
	{
		for (size_t j = 1; j < segments[i].size(); j++)
		{
			std::pair<size_t, bool> edgeFrom;
			edgeFrom.first = segmentToNode.at(std::make_pair(std::get<0>(segments[i][j-1]), std::get<1>(segments[i][j-1])));
			edgeFrom.second = std::get<2>(segments[i][j-1]);
			std::pair<size_t, bool> edgeTo;
			edgeTo.first = segmentToNode.at(std::make_pair(std::get<0>(segments[i][j]), std::get<1>(segments[i][j])));
			edgeTo.second = std::get<2>(segments[i][j]);
			auto key = canon(edgeFrom, edgeTo);
			edgeCoverage[key] += 1;
		}
	}
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

std::pair<size_t, size_t> getMatchSpan(const std::vector<size_t>& readLengths, const std::tuple<size_t, size_t, size_t, size_t, size_t, size_t, bool> match, size_t readI)
{
	if (std::get<0>(match) == readI)
	{
		assert(std::get<3>(match) != readI);
		return std::make_pair(std::get<1>(match), std::get<2>(match));
	}
	else
	{
		assert(std::get<3>(match) == readI);
		if (std::get<6>(match))
		{
			return std::make_pair(std::get<4>(match), std::get<5>(match));
		}
		else
		{
			return std::make_pair(readLengths[readI] - std::get<5>(match), readLengths[readI] - std::get<4>(match));
		}
	}
}

void buildIntervalTree(const std::vector<size_t>& readLengths, const std::vector<std::tuple<size_t, size_t, size_t, size_t, size_t, size_t, bool>>& matches, std::vector<size_t>& relevantMatches, std::vector<size_t>& relevantMatchTreeMaxRight, std::vector<size_t>& relevantMatchTreeMinLeft, const size_t readi)
{
	const uint64_t firstBit = 1ull << 63ull;
	const uint64_t mask = firstBit-1;
	std::vector<size_t> sortedMatches = relevantMatches;
	std::sort(sortedMatches.begin(), sortedMatches.end(), [&matches, &readLengths, mask, firstBit, readi](size_t left, size_t right) { return getMatchSpan(readLengths, matches[left & mask], readi).first < getMatchSpan(readLengths, matches[right & mask], readi).first; });
	size_t test = sortedMatches.size();
	relevantMatches.assign(relevantMatches.size(), std::numeric_limits<size_t>::max());
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
		relevantMatchTreeMaxRight[i] = getMatchSpan(readLengths, matches[relevantMatches[i] & mask], readi).second;
		relevantMatchTreeMinLeft[i] = getMatchSpan(readLengths, matches[relevantMatches[i] & mask], readi).first;
		if (i * 2 + 1 < sortedMatches.size())
		{
			relevantMatchTreeMaxRight[i] = std::max(relevantMatchTreeMaxRight[i], relevantMatchTreeMaxRight[i*2+1]);
			relevantMatchTreeMinLeft[i] = std::min(relevantMatchTreeMinLeft[i], relevantMatchTreeMinLeft[i*2+1]);
			assert(getMatchSpan(readLengths, matches[relevantMatches[i*2+1] & mask], readi).first <= getMatchSpan(readLengths, matches[relevantMatches[i] & mask], readi).first);
		}
		if (i * 2 + 2 < sortedMatches.size())
		{
			relevantMatchTreeMaxRight[i] = std::max(relevantMatchTreeMaxRight[i], relevantMatchTreeMaxRight[i*2+2]);
			relevantMatchTreeMinLeft[i] = std::min(relevantMatchTreeMinLeft[i], relevantMatchTreeMinLeft[i*2+2]);
			assert(getMatchSpan(readLengths, matches[relevantMatches[i] & mask], readi).first <= getMatchSpan(readLengths, matches[relevantMatches[i*2+2] & mask], readi).first);
		}
	}
}

template <typename F>
void iterateIntervalTreeAlnMatches(const std::vector<size_t>& readLengths, const std::vector<std::tuple<size_t, size_t, size_t, size_t, size_t, size_t, bool>>& matches, const std::vector<size_t>& relevantMatches, const std::vector<std::pair<size_t, size_t>>& matchLeftSpan, const std::vector<std::pair<size_t, size_t>>& matchRightSpan, const std::vector<size_t>& relevantMatchTreeMaxRight, const std::vector<size_t>& relevantMatchTreeMinLeft, const size_t readi, const size_t intervalStart, const size_t intervalEnd, F callback)
{
	const uint64_t firstBit = 1ull << 63ull;
	const uint64_t mask = firstBit-1;
	if (relevantMatchTreeMaxRight[0] < intervalStart) return;
	if (relevantMatchTreeMinLeft[0] > intervalEnd) return;
	// weird stack structure so no function call recursion, weird pos / stack insert to minimize stack pushes
	std::vector<size_t> checkStack;
	size_t pos = 0;
	while (true)
	{
		size_t aln = relevantMatches[pos];
		size_t rawAln = aln & mask;
		auto interval = getMatchSpan(readLengths, matches[rawAln], readi);
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

std::vector<RankBitvector> extendBreakpoints(const std::vector<size_t>& readLengths, const std::vector<std::tuple<size_t, size_t, size_t, size_t, size_t, size_t, bool>>& matches)
{
	std::vector<RankBitvector> breakpoints;
	const uint64_t firstBit = 1ull << 63ull;
	const uint64_t mask = firstBit-1;
	breakpoints.resize(readLengths.size());
	for (size_t i = 0; i < readLengths.size(); i++)
	{
		breakpoints[i].resize(readLengths[i]+1);
		breakpoints[i].set(0, true);
		breakpoints[i].set(readLengths[i], true);
	}
	std::vector<std::vector<uint64_t>> relevantMatches;
	relevantMatches.resize(readLengths.size());
	size_t firstFwBwMatch = matches.size();
	for (size_t i = 0; i < matches.size(); i++)
	{
		if (std::get<6>(matches[i]))
		{
			assert(firstFwBwMatch == matches.size());
		}
		else
		{
			if (firstFwBwMatch == matches.size())
			{
				firstFwBwMatch = i;
			}
			assert(i >= firstFwBwMatch);
		}
		relevantMatches[std::get<0>(matches[i])].emplace_back(i);
		relevantMatches[std::get<3>(matches[i])].emplace_back(i + firstBit);
		breakpoints[std::get<0>(matches[i])].set(std::get<1>(matches[i]), true);
		breakpoints[std::get<0>(matches[i])].set(std::get<2>(matches[i]), true);
		if (std::get<6>(matches[i]))
		{
			breakpoints[std::get<3>(matches[i])].set(std::get<4>(matches[i]), true);
			breakpoints[std::get<3>(matches[i])].set(std::get<5>(matches[i]), true);
		}
		else
		{
			breakpoints[std::get<3>(matches[i])].set(readLengths[std::get<3>(matches[i])], true);
			breakpoints[std::get<3>(matches[i])].set(readLengths[std::get<3>(matches[i])], true);
		}
	}
	std::vector<std::vector<size_t>> relevantMatchTreeMaxRight;
	std::vector<std::vector<size_t>> relevantMatchTreeMinLeft;
	relevantMatchTreeMaxRight.resize(relevantMatches.size());
	relevantMatchTreeMinLeft.resize(relevantMatches.size());
	for (size_t i = 0; i < relevantMatches.size(); i++)
	{
		buildIntervalTree(readLengths, matches, relevantMatches[i], relevantMatchTreeMaxRight[i], relevantMatchTreeMinLeft[i], i);
	}
	std::vector<std::pair<size_t, size_t>> matchLeftSpan;
	std::vector<std::pair<size_t, size_t>> matchRightSpan;
	matchLeftSpan.resize(matches.size());
	matchRightSpan.resize(matches.size());
	for (size_t i = 0; i < matches.size(); i++)
	{
		matchLeftSpan[i] = getMatchSpan(readLengths, matches[i], std::get<0>(matches[i]));
		matchRightSpan[i] = getMatchSpan(readLengths, matches[i], std::get<3>(matches[i]));
	}
	std::vector<bool> inQueue;
	inQueue.resize(matches.size(), true);
	std::vector<size_t> fwfwCheckQueue;
	std::vector<size_t> fwbwCheckQueue;
	for (size_t i = 0; i < firstFwBwMatch; i++)
	{
		fwfwCheckQueue.emplace_back(i);
	}
	for (size_t i = firstFwBwMatch; i < matches.size(); i++)
	{
		fwbwCheckQueue.emplace_back(i);
	}
	while (fwfwCheckQueue.size() >= 1 || fwbwCheckQueue.size() >= 1)
	{
		size_t matchi;
		std::pair<bool, bool> addedAny;
		if (fwfwCheckQueue.size() >= 1)
		{
			matchi = fwfwCheckQueue.back();
			assert(std::get<6>(matches[matchi]));
			inQueue[matchi] = false;
			fwfwCheckQueue.pop_back();
			addedAny = extendBreakpointsFwFw(readLengths, breakpoints, std::get<0>(matches[matchi]), std::get<1>(matches[matchi]), std::get<2>(matches[matchi]), std::get<3>(matches[matchi]), std::get<4>(matches[matchi]), std::get<5>(matches[matchi]));
		}
		else
		{
			assert(fwbwCheckQueue.size() >= 1);
			matchi = fwbwCheckQueue.back();
			assert(!std::get<6>(matches[matchi]));
			inQueue[matchi] = false;
			fwbwCheckQueue.pop_back();
			addedAny = extendBreakpointsFwBw(readLengths, breakpoints, std::get<0>(matches[matchi]), std::get<1>(matches[matchi]), std::get<2>(matches[matchi]), std::get<3>(matches[matchi]), std::get<4>(matches[matchi]), std::get<5>(matches[matchi]));
		}
		if (addedAny.first)
		{
			iterateIntervalTreeAlnMatches(readLengths, matches, relevantMatches[std::get<0>(matches[matchi])], matchLeftSpan, matchRightSpan, relevantMatchTreeMaxRight[std::get<0>(matches[matchi])], relevantMatchTreeMinLeft[std::get<0>(matches[matchi])], std::get<0>(matches[matchi]), matchLeftSpan[matchi].first, matchLeftSpan[matchi].second, [firstBit, mask, firstFwBwMatch, &fwfwCheckQueue, &fwbwCheckQueue, &inQueue, &matchLeftSpan, &matchRightSpan, matchi](size_t aln)
			{
				size_t rawAln = aln & mask;
				if (inQueue[rawAln]) return;
				if (rawAln == matchi) return;
				bool isRight = (aln & firstBit) != 0;
				assert(!isRight || matchRightSpan[rawAln].first <= matchLeftSpan[matchi].second);
				assert(!isRight || matchRightSpan[rawAln].second >= matchLeftSpan[matchi].first);
				assert(isRight || matchLeftSpan[rawAln].first <= matchLeftSpan[matchi].second);
				assert(isRight || matchLeftSpan[rawAln].second >= matchLeftSpan[matchi].first);
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
			iterateIntervalTreeAlnMatches(readLengths, matches, relevantMatches[std::get<3>(matches[matchi])], matchLeftSpan, matchRightSpan, relevantMatchTreeMaxRight[std::get<3>(matches[matchi])], relevantMatchTreeMinLeft[std::get<3>(matches[matchi])], std::get<3>(matches[matchi]), matchRightSpan[matchi].first, matchRightSpan[matchi].second, [firstBit, mask, firstFwBwMatch, &fwfwCheckQueue, &fwbwCheckQueue, &inQueue, &matchLeftSpan, &matchRightSpan, matchi](size_t aln)
			{
				size_t rawAln = aln & mask;
				if (inQueue[rawAln]) return;
				if (rawAln == matchi) return;
				bool isRight = (aln & firstBit) != 0;
				assert(!isRight || matchRightSpan[rawAln].first <= matchRightSpan[matchi].second);
				assert(!isRight || matchRightSpan[rawAln].second >= matchRightSpan[matchi].first);
				assert(isRight || matchLeftSpan[rawAln].first <= matchRightSpan[matchi].second);
				assert(isRight || matchLeftSpan[rawAln].second >= matchRightSpan[matchi].first);
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
	for (size_t i = 0; i < breakpoints.size(); i++) breakpoints[i].buildRanks();
	return breakpoints;
}

std::vector<size_t> getNodeCoverage(const std::vector<std::vector<std::tuple<uint32_t, uint32_t, bool>>>& segments, const phmap::flat_hash_map<std::pair<uint32_t, uint32_t>, size_t>& segmentToNode)
{
	std::vector<size_t> result;
	result.resize(segmentToNode.size(), 0);
	for (size_t i = 0; i < segments.size(); i++)
	{
		for (size_t j = 0; j < segments[i].size(); j++)
		{
			size_t node = segmentToNode.at(std::make_pair(std::get<0>(segments[i][j]), std::get<1>(segments[i][j])));
			result[node] += 1;
		}
	}
	return result;
}

std::vector<size_t> getNodeLengths(const std::vector<std::vector<std::tuple<uint32_t, uint32_t, bool>>>& segments, const phmap::flat_hash_map<std::pair<uint32_t, uint32_t>, size_t>& segmentToNode, const std::vector<std::vector<uint32_t>>& segmentStarts)
{
	std::vector<size_t> result;
	result.resize(segmentToNode.size(), std::numeric_limits<size_t>::max());
	for (size_t i = 0; i < segments.size(); i++)
	{
		for (size_t j = 0; j < segments[i].size(); j++)
		{
			size_t node = segmentToNode.at(std::make_pair(std::get<0>(segments[i][j]), std::get<1>(segments[i][j])));
			size_t length = segmentStarts[i][j+1] - segmentStarts[i][j];
			if (!(result[node] == std::numeric_limits<size_t>::max() || result[node] == length))
			{
				std::cerr << result[node] << " " << length << std::endl;
			}
			assert(result[node] == std::numeric_limits<size_t>::max() || result[node] == length);
			result[node] = length;
		}
	}
	return result;
}

void makeGraph(const std::vector<size_t>& readLengths, const std::vector<std::tuple<size_t, size_t, size_t, size_t, size_t, size_t, bool>>& matches, const size_t minCoverage, const std::string& outputFileName, const size_t k)
{
	std::vector<RankBitvector> breakpoints = extendBreakpoints(readLengths, matches);
	std::vector<std::vector<std::tuple<uint32_t, uint32_t, bool>>> segments;
	std::vector<std::vector<uint32_t>> segmentStarts;
	std::tie(segments, segmentStarts) = mergeSegments(readLengths, matches, breakpoints);
	phmap::flat_hash_map<std::pair<uint32_t, uint32_t>, size_t> segmentToNode = getSegmentToNode(segments);
	std::vector<size_t> nodeCoverage = getNodeCoverage(segments, segmentToNode);
	std::vector<size_t> nodeLength = getNodeLengths(segments, segmentToNode, segmentStarts);
	phmap::flat_hash_map<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>, size_t> edgeCoverages = getEdgeCoverages(readLengths, segmentToNode, segments);
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

std::vector<std::tuple<size_t, size_t, size_t, size_t, size_t, size_t, bool>> getKmerMatches(const std::vector<TwobitString>& readSequences, const size_t left, const size_t leftStart, const size_t leftEnd, const size_t right, const size_t rightStart, const size_t rightEnd, const bool rightFw, const size_t k, const size_t w)
{
	assert(k <= 31);
	assert(k % 2 == 1);
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
	std::vector<std::tuple<size_t, size_t, size_t, size_t, size_t, size_t, bool>> result;
	iterateSyncmers(readSequences, k, 20, right, rightStart, rightEnd, rightFw, [&firstKmerPositionInLeft, &extraKmerPositionsInLeft, &currentMatchesPerDiagonal, diagonalCount, zeroDiagonal, rightFw, left, right, leftStart, rightStart, leftEnd, rightEnd, w, k, &result](const size_t kmer, const size_t rightPos)
	{
		size_t interpolatedLeftPos = (double)(rightPos) / (double)(rightEnd-rightStart) * (double)(leftEnd-leftStart);
		assert(rightPos + zeroDiagonal >= interpolatedLeftPos + w);
		assert(rightPos + zeroDiagonal + w >= interpolatedLeftPos);
		size_t minDiagonal = rightPos + zeroDiagonal - interpolatedLeftPos - w;
		size_t maxDiagonal = rightPos + zeroDiagonal - interpolatedLeftPos + w;
		assert(maxDiagonal <= diagonalCount);
		iterateKmerMatchPositions(kmer, firstKmerPositionInLeft, extraKmerPositionsInLeft, [zeroDiagonal, rightPos, minDiagonal, maxDiagonal, &currentMatchesPerDiagonal, &result, rightFw, left, right, leftStart, rightStart, diagonalCount, k](const size_t leftPos)
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
				result.emplace_back(left, leftMatchStart, leftMatchStart + length, right, rightMatchStart, rightMatchStart + length, rightFw);
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
		result.emplace_back(left, leftMatchStart, leftMatchStart + length, right, rightMatchStart, rightMatchStart + length, rightFw);
	}
	return result;
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
	std::vector<std::tuple<size_t, size_t, size_t, size_t, size_t, size_t, bool>> mappingmatches;
	auto result = matchIndex.iterateMatchChains(numThreads, storage.getRawReadLengths(), [&printMutex, &mappingmatches, minAlignmentLength](const size_t left, const size_t leftstart, const size_t leftend, const bool leftFw, const size_t right, const size_t rightstart, const size_t rightend, const bool rightFw)
	{
		if (leftend+1-leftstart < minAlignmentLength) return;
		if (rightend+1-rightstart < minAlignmentLength) return;
		std::lock_guard<std::mutex> lock { printMutex };
		assert(leftFw);
		mappingmatches.emplace_back(left, leftstart, leftend+1, right, rightstart, rightend+1, rightFw); // param coordinates are inclusive, switch to exclusive
	});
	// fw-fw matches first, fw-bw matches later
	std::stable_sort(mappingmatches.begin(), mappingmatches.end(), [](auto left, auto right){ return std::get<6>(left) && !std::get<6>(right); });
	std::vector<std::tuple<size_t, size_t, size_t, size_t, size_t, size_t, bool>> kmermatches;
	std::cerr << mappingmatches.size() << " mapping matches" << std::endl;
	for (auto t : mappingmatches)
	{
		auto matches = getKmerMatches(readSequences, std::get<0>(t), std::get<1>(t), std::get<2>(t), std::get<3>(t), std::get<4>(t), std::get<5>(t), std::get<6>(t), graphk, 50);
		for (auto match : matches)
		{
			assert(std::get<2>(match) > std::get<1>(match));
			assert(std::get<1>(match) >= std::get<1>(t));
			assert(std::get<2>(match) <= std::get<2>(t));
			assert(std::get<5>(match) > std::get<4>(match));
			assert(std::get<4>(match) >= std::get<4>(t));
			assert(std::get<5>(match) <= std::get<5>(t));
			// std::cout << readNames[std::get<0>(match)] << "\t" << readKmerLengths[std::get<0>(match)] << "\t" << std::get<1>(match) << "\t" << std::get<2>(match) << "\t" << readNames[std::get<3>(match)] << "\t" << readKmerLengths[std::get<3>(match)] << "\t" << std::get<4>(match) << "\t" << std::get<5>(match) << "\t" << (std::get<6>(match) ? "fw" : "bw") << std::endl;
		}
		kmermatches.insert(kmermatches.end(), matches.begin(), matches.end());
	}
	std::cerr << kmermatches.size() << " kmer matches" << std::endl;
	makeGraph(readKmerLengths, kmermatches, minCoverage, "graph.gfa", graphk);
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
