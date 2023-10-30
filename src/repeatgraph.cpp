#include <mutex>
#include <iostream>
#include <vector>
#include <string>
#include "CommonUtils.h"
#include "ReadStorage.h"
#include "MatchIndex.h"

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

phmap::flat_hash_map<std::pair<size_t, size_t>, size_t> getBreakpointToNode(const std::vector<std::vector<std::tuple<size_t, size_t, size_t, bool>>>& breakpoints)
{
	size_t nextIndex = 0;
	phmap::flat_hash_map<std::pair<size_t, size_t>, size_t> result;
	for (size_t i = 0; i < breakpoints.size(); i++)
	{
		for (size_t j = 0; j < breakpoints[i].size(); j++)
		{
			if (std::get<1>(breakpoints[i][j]) != i) continue;
			if (std::get<2>(breakpoints[i][j]) != std::get<0>(breakpoints[i][j])) continue;
			result[std::make_pair(i, std::get<0>(breakpoints[i][j]))] = nextIndex;
			nextIndex += 1;
		}
	}
	return result;
}

std::vector<size_t> getBreakpointCoverage(const phmap::flat_hash_map<std::pair<size_t, size_t>, size_t>& breakpointToNode, const std::vector<std::vector<std::tuple<size_t, size_t, size_t, bool>>>& breakpoints)
{
	std::vector<size_t> result;
	result.resize(breakpointToNode.size(), 0);
	for (size_t i = 0; i < breakpoints.size(); i++)
	{
		for (size_t j = 0; j < breakpoints[i].size(); j++)
		{
			size_t index = breakpointToNode.at(std::make_pair(std::get<1>(breakpoints[i][j]), std::get<2>(breakpoints[i][j])));
			assert(index < result.size());
			result[index] += 1;
		}
	}
	return result;
}

std::tuple<size_t, size_t, bool> find(std::vector<std::vector<std::tuple<size_t, size_t, bool>>>& result, const size_t read, const size_t index)
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

void mergeSegments(std::vector<std::vector<std::tuple<size_t, size_t, bool>>>& result, const size_t leftRead, const size_t leftIndex, const size_t rightRead, const size_t rightIndex, const bool fw)
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

void mergeSegments(const std::vector<size_t>& readLengths, std::vector<std::vector<std::tuple<size_t, size_t, bool>>>& result, const std::vector<std::vector<size_t>>& segmentStarts, const size_t leftRead, const bool leftFw, const size_t leftStart, const size_t leftEnd, const size_t rightRead, const bool rightFw, const size_t rightStart, const size_t rightEnd)
{
	assert(leftFw);
	size_t leftFirst = std::numeric_limits<size_t>::max();
	size_t leftLast = std::numeric_limits<size_t>::max();
	for (size_t i = 0; i < segmentStarts[leftRead].size(); i++)
	{
		if (segmentStarts[leftRead][i] == leftStart)
		{
			assert(leftFirst == std::numeric_limits<size_t>::max());
			leftFirst = i;
		}
		if (segmentStarts[leftRead][i] == leftEnd)
		{
			assert(leftLast == std::numeric_limits<size_t>::max());
			leftLast = i;
		}
	}
	assert(leftFirst != std::numeric_limits<size_t>::max());
	assert(leftLast != std::numeric_limits<size_t>::max());
	assert(leftLast > leftFirst);
	size_t rightFirst = std::numeric_limits<size_t>::max();
	size_t rightLast = std::numeric_limits<size_t>::max();
	for (size_t i = 0; i < segmentStarts[rightRead].size(); i++)
	{
		if (rightFw && segmentStarts[rightRead][i] == rightStart)
		{
			assert(rightFirst == std::numeric_limits<size_t>::max());
			rightFirst = i;
		}
		if (rightFw && segmentStarts[rightRead][i] == rightEnd)
		{
			assert(rightLast == std::numeric_limits<size_t>::max());
			rightLast = i;
		}
		if (!rightFw && segmentStarts[rightRead][i] == readLengths[rightRead]-rightEnd)
		{
			assert(rightFirst == std::numeric_limits<size_t>::max());
			rightFirst = i;
		}
		if (!rightFw && segmentStarts[rightRead][i] == readLengths[rightRead]-rightStart)
		{
			assert(rightLast == std::numeric_limits<size_t>::max());
			rightLast = i;
		}
	}
	assert(rightFirst != std::numeric_limits<size_t>::max());
	assert(rightLast != std::numeric_limits<size_t>::max());
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

std::tuple<std::vector<std::vector<std::tuple<size_t, size_t, bool>>>, std::vector<std::vector<size_t>>> mergeSegments(const std::vector<size_t>& readLengths, const std::vector<std::tuple<size_t, size_t, size_t, size_t, size_t, size_t, bool>>& matches, const std::vector<std::vector<bool>>& breakpoints)
{
	assert(breakpoints.size() == readLengths.size());
	std::vector<std::vector<std::tuple<size_t, size_t, bool>>> result;
	std::vector<std::vector<size_t>> segmentStarts;
	result.resize(breakpoints.size());
	segmentStarts.resize(breakpoints.size());
	for (size_t i = 0; i < breakpoints.size(); i++)
	{
		segmentStarts[i].push_back(0);
		assert(breakpoints[i][0]);
		for (size_t j = 1; j < breakpoints[i].size(); j++)
		{
			if (!breakpoints[i][j]) continue;
			result[i].emplace_back(i, result[i].size(), true);
			segmentStarts[i].push_back(j);
		}
		assert(segmentStarts[i].size() == result[i].size()+1);
		assert(segmentStarts[i][0] == 0);
		assert(segmentStarts[i].back() == readLengths[i]);
		assert(result[i].size() >= 1);
	}
	for (auto match : matches)
	{
		mergeSegments(readLengths, result, segmentStarts, std::get<0>(match), true, std::get<1>(match), std::get<2>(match), std::get<3>(match), std::get<6>(match), std::get<4>(match), std::get<5>(match));
	}
	for (size_t i = 0; i < result.size(); i++)
	{
		for (size_t j = 0; j < result[i].size(); j++)
		{
			find(result, i, j);
		}
	}
	return std::make_tuple(result, segmentStarts);
}

phmap::flat_hash_map<std::pair<size_t, size_t>, size_t> getSegmentToNode(const std::vector<std::vector<std::tuple<size_t, size_t, bool>>>& segments)
{
	size_t nextIndex = 0;
	phmap::flat_hash_map<std::pair<size_t, size_t>, size_t> result;
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

phmap::flat_hash_map<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>, size_t> getEdgeCoverages(const std::vector<size_t>& readLengths, const phmap::flat_hash_map<std::pair<size_t, size_t>, size_t>& segmentToNode, const std::vector<std::vector<std::tuple<size_t, size_t, bool>>>& segments)
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

std::pair<bool, bool> extendBreakpoints(const std::vector<size_t>& readLengths, std::vector<std::vector<bool>>& breakpoints, size_t leftRead, bool leftFw, size_t leftStart, size_t leftEnd, size_t rightRead, bool rightFw, size_t rightStart, size_t rightEnd)
{
	assert(leftFw);
	assert(leftRead != rightRead);
	assert(leftEnd - leftStart == rightEnd - rightStart);
	bool addedAnyLeft = false;
	bool addedAnyRight = false;
	size_t leftIndex = leftStart;
	size_t rightIndex = rightStart;
	if (!rightFw) rightIndex = readLengths[rightRead] - rightStart;
	for (size_t i = 0; i <= leftEnd-leftStart; i++)
	{
		bool leftbit = breakpoints[leftRead][leftIndex];
		bool rightbit = breakpoints[rightRead][rightIndex];
		if (leftbit && !rightbit)
		{
			breakpoints[rightRead][rightIndex] = true;
			addedAnyRight = true;
		}
		else if (!leftbit && rightbit)
		{
			breakpoints[leftRead][leftIndex] = true;
			addedAnyLeft = true;
		}
		leftIndex += 1;
		rightIndex += (rightFw ? 1 : -1);
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

std::vector<std::vector<bool>> extendBreakpoints(const std::vector<size_t>& readLengths, const std::vector<std::tuple<size_t, size_t, size_t, size_t, size_t, size_t, bool>>& matches)
{
	std::vector<std::vector<bool>> breakpoints;
	uint64_t firstBit = 1ull << 63ull;
	uint64_t mask = firstBit-1;
	breakpoints.resize(readLengths.size());
	for (size_t i = 0; i < readLengths.size(); i++)
	{
		breakpoints[i].resize(readLengths[i]+1, false);
		breakpoints[i][0] = true;
		breakpoints[i].back() = true;
	}
	std::vector<std::vector<uint64_t>> relevantMatches;
	std::vector<std::pair<size_t, size_t>> matchLeftSpan;
	std::vector<std::pair<size_t, size_t>> matchRightSpan;
	std::vector<bool> inQueue;
	inQueue.resize(matches.size(), true);
	matchLeftSpan.resize(matches.size());
	matchRightSpan.resize(matches.size());
	relevantMatches.resize(readLengths.size());
	for (size_t i = 0; i < matches.size(); i++)
	{
		relevantMatches[std::get<0>(matches[i])].emplace_back(i);
		relevantMatches[std::get<3>(matches[i])].emplace_back(i + firstBit);
		matchLeftSpan[i] = getMatchSpan(readLengths, matches[i], std::get<0>(matches[i]));
		matchRightSpan[i] = getMatchSpan(readLengths, matches[i], std::get<3>(matches[i]));
		breakpoints[std::get<0>(matches[i])][std::get<1>(matches[i])] = true;
		breakpoints[std::get<0>(matches[i])][std::get<2>(matches[i])] = true;
		if (std::get<6>(matches[i]))
		{
			breakpoints[std::get<3>(matches[i])][std::get<4>(matches[i])] = true;
			breakpoints[std::get<3>(matches[i])][std::get<5>(matches[i])] = true;
		}
		else
		{
			breakpoints[std::get<3>(matches[i])][readLengths[std::get<3>(matches[i])]] = true;
			breakpoints[std::get<3>(matches[i])][readLengths[std::get<3>(matches[i])]] = true;
		}
	}
	std::vector<size_t> checkQueue;
	for (size_t i = 0; i < matches.size(); i++)
	{
		checkQueue.emplace_back(i);
	}
	while (checkQueue.size() >= 1)
	{
		size_t matchi = checkQueue.back();
		inQueue[matchi] = false;
		checkQueue.pop_back();
		std::pair<bool, bool> addedAny = extendBreakpoints(readLengths, breakpoints, std::get<0>(matches[matchi]), true, std::get<1>(matches[matchi]), std::get<2>(matches[matchi]), std::get<3>(matches[matchi]), std::get<6>(matches[matchi]), std::get<4>(matches[matchi]), std::get<5>(matches[matchi]));
		if (addedAny.first)
		{
			for (auto aln : relevantMatches[std::get<0>(matches[matchi])])
			{
				size_t rawAln = aln & mask;
				if (inQueue[rawAln]) continue;
				if (rawAln == matchi) continue;
				bool isRight = (aln & firstBit) != 0;
				if (isRight && (matchRightSpan[rawAln].first > matchLeftSpan[matchi].second || matchRightSpan[rawAln].second < matchLeftSpan[matchi].first)) continue;
				if (!isRight && (matchLeftSpan[rawAln].first > matchLeftSpan[matchi].second || matchLeftSpan[rawAln].second < matchLeftSpan[matchi].first)) continue;
				checkQueue.push_back(rawAln);
				inQueue[rawAln] = true;
			}
		}
		if (addedAny.second)
		{
			for (auto aln : relevantMatches[std::get<3>(matches[matchi])])
			{
				size_t rawAln = aln & mask;
				if (inQueue[rawAln]) continue;
				if (rawAln == matchi) continue;
				bool isRight = (aln & firstBit) != 0;
				if (isRight && (matchRightSpan[rawAln].first > matchRightSpan[matchi].second || matchRightSpan[rawAln].second < matchRightSpan[matchi].first)) continue;
				if (!isRight && (matchLeftSpan[rawAln].first > matchRightSpan[matchi].second || matchLeftSpan[rawAln].second < matchRightSpan[matchi].first)) continue;
				checkQueue.push_back(rawAln);
				inQueue[rawAln] = true;
			}
		}
	}
	return breakpoints;
}

std::vector<size_t> getNodeCoverage(const std::vector<std::vector<std::tuple<size_t, size_t, bool>>>& segments, const phmap::flat_hash_map<std::pair<size_t, size_t>, size_t>& segmentToNode)
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

std::vector<size_t> getNodeLengths(const std::vector<std::vector<std::tuple<size_t, size_t, bool>>>& segments, const phmap::flat_hash_map<std::pair<size_t, size_t>, size_t>& segmentToNode, const std::vector<std::vector<size_t>>& segmentStarts)
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
	std::vector<std::vector<bool>> breakpoints = extendBreakpoints(readLengths, matches);
	std::vector<std::vector<std::tuple<size_t, size_t, bool>>> segments;
	std::vector<std::vector<size_t>> segmentStarts;
	std::tie(segments, segmentStarts) = mergeSegments(readLengths, matches, breakpoints);
	phmap::flat_hash_map<std::pair<size_t, size_t>, size_t> segmentToNode = getSegmentToNode(segments);
	std::vector<size_t> nodeCoverage = getNodeCoverage(segments, segmentToNode);
	std::vector<size_t> nodeLength = getNodeLengths(segments, segmentToNode, segmentStarts);
	phmap::flat_hash_map<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>, size_t> edgeCoverages = getEdgeCoverages(readLengths, segmentToNode, segments);
	writeGraph(outputFileName, nodeCoverage, nodeLength, edgeCoverages, minCoverage, k);
}

std::vector<std::tuple<size_t, size_t, size_t, size_t, size_t, size_t, bool>> getKmerMatches(const std::vector<TwobitString>& readSequences, const size_t left, const size_t leftStart, const size_t leftEnd, const size_t right, const size_t rightStart, const size_t rightEnd, const bool rightFw, const size_t k, const size_t w)
{
	assert(k <= 31);
	assert(k % 2 == 1);
	phmap::flat_hash_map<uint64_t, std::vector<size_t>> kmerPositionsInLeft;
	uint64_t mask = (1ull << (2ull*k)) - 1;
	uint64_t leftKmer = 0;
	for (size_t i = 0; i < k; i++)
	{
		leftKmer <<= 2;
		leftKmer += readSequences[left].get(leftStart+i);
	}
	kmerPositionsInLeft[leftKmer].push_back(0);
	std::vector<std::tuple<size_t, size_t, size_t, size_t, size_t, size_t, bool>> result;
	for (size_t i = k; i < leftEnd-leftStart; i++)
	{
		leftKmer <<= 2;
		leftKmer &= mask;
		leftKmer += readSequences[left].get(leftStart+i);
		kmerPositionsInLeft[leftKmer].push_back(i-k+1);
	}
	uint64_t rightKmer = 0;
	for (size_t i = 0; i < k; i++)
	{
		rightKmer <<= 2;
		if (rightFw)
		{
			rightKmer += readSequences[right].get(rightStart+i);
		}
		else
		{
			rightKmer += 3-readSequences[right].get(readSequences[right].size()-1-(rightStart+i));
		}
	}
	std::vector<std::pair<size_t, size_t>> currentMatches;
	for (auto pos : kmerPositionsInLeft[rightKmer])
	{
		size_t interpolatedPos = (double)pos / (double)(leftEnd-leftStart) * (double)(rightEnd-rightStart);
		if (interpolatedPos > w) continue;
		currentMatches.emplace_back(pos, 1);
	}
	for (size_t rightpos = k; rightpos < rightEnd-rightStart; rightpos++)
	{
		rightKmer <<= 2;
		rightKmer &= mask;
		if (rightFw)
		{
			rightKmer += readSequences[right].get(rightStart+rightpos);
		}
		else
		{
			rightKmer += 3-readSequences[right].get(readSequences[right].size()-1-(rightStart+rightpos));
		}
		phmap::flat_hash_set<size_t> poses;
		for (auto pos : kmerPositionsInLeft[rightKmer])
		{
			size_t interpolatedPos = (double)pos / (double)(leftEnd-leftStart) * (double)(rightEnd-rightStart);
			if (interpolatedPos > rightpos-k+1+w) continue;
			if (rightpos-k+1 > interpolatedPos+w) continue;
			poses.emplace(pos);
		}
		phmap::flat_hash_set<size_t> foundPoses;
		for (size_t i = currentMatches.size()-1; i < currentMatches.size(); i--)
		{
			if (poses.count(currentMatches[i].first+currentMatches[i].second) == 0)
			{
				result.emplace_back(left, leftStart + currentMatches[i].first, leftStart + currentMatches[i].first + currentMatches[i].second, right, rightStart + rightpos + 1 - k - currentMatches[i].second, rightStart + rightpos + 1 - k, rightFw);
				std::swap(currentMatches.back(), currentMatches[i]);
				currentMatches.pop_back();
			}
			else
			{
				foundPoses.insert(currentMatches[i].first+currentMatches[i].second);
				currentMatches[i].second += 1;
			}
		}
		for (auto pos : poses)
		{
			if (foundPoses.count(pos) == 1) continue;
			currentMatches.emplace_back(pos, 1);
		}
	}
	for (auto match : currentMatches)
	{
		result.emplace_back(left, leftStart + match.first, leftStart + match.first + match.second, right, rightStart + (rightEnd-rightStart) + 1 - k - match.second, rightStart + (rightEnd-rightStart) + 1 - k, rightFw);
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
	const std::vector<size_t>& readLengths = storage.getRawReadLengths();
	const std::vector<std::string>& readNames = storage.getNames();
	std::mutex printMutex;
	std::vector<std::tuple<size_t, size_t, size_t, size_t, size_t, size_t, bool>> mappingmatches;
	auto result = matchIndex.iterateMatchChains(numThreads, storage.getRawReadLengths(), [&printMutex, &mappingmatches, minAlignmentLength, &readLengths](const size_t left, const size_t leftstart, const size_t leftend, const bool leftFw, const size_t right, const size_t rightstart, const size_t rightend, const bool rightFw)
	{
		if (leftend+1-leftstart < minAlignmentLength) return;
		if (rightend+1-rightstart < minAlignmentLength) return;
		std::lock_guard<std::mutex> lock { printMutex };
		assert(leftFw);
		mappingmatches.emplace_back(left, leftstart, leftend+1, right, rightstart, rightend+1, rightFw); // param coordinates are inclusive, switch to exclusive
	});
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
		}
		kmermatches.insert(kmermatches.end(), matches.begin(), matches.end());
	}
	std::cerr << kmermatches.size() << " kmer matches" << std::endl;
	makeGraph(readLengths, kmermatches, minCoverage, "graph.gfa", graphk);
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
