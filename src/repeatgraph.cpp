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

bool extendAndMergeBreakpoints(const std::vector<size_t>& readLengths, std::vector<std::vector<std::tuple<size_t, size_t, size_t, bool>>>& breakpoints, size_t leftRead, bool leftFw, size_t leftStart, size_t leftEnd, size_t rightRead, bool rightFw, size_t rightStart, size_t rightEnd)
{
	assert(leftFw);
	assert(leftRead != rightRead);
	assert(leftEnd - leftStart == rightEnd - rightStart);
	bool addedAny = false;
	size_t countLeftBreakpoints = 0;
	size_t countRightBreakpoints = 0;
	for (auto breakpoint : breakpoints[leftRead])
	{
		if (std::get<0>(breakpoint) < leftStart) continue;
		if (std::get<0>(breakpoint) > leftEnd) continue;
		countLeftBreakpoints += 1;
		size_t positionWithinMatch = std::get<0>(breakpoint) - leftStart;
		assert(positionWithinMatch <= leftEnd - leftStart);
		size_t positionInRight = rightStart + positionWithinMatch;
		assert(positionInRight <= readLengths[rightRead]);
		if (rightFw)
		{
			auto gotMatch = getBreakpoint(breakpoints[rightRead], rightRead, positionInRight);
			if (gotMatch.second) addedAny = true;
			bool merged = merge(breakpoints, leftRead, std::get<0>(breakpoint), rightRead, gotMatch.first, true);
			if (merged) addedAny = true;
		}
		else
		{
			positionInRight = readLengths[rightRead] - positionInRight;
			assert(positionInRight <= readLengths[rightRead]);
			auto gotMatch = getBreakpoint(breakpoints[rightRead], rightRead, positionInRight);
			if (gotMatch.second) addedAny = true;
			bool merged = merge(breakpoints, leftRead, std::get<0>(breakpoint), rightRead, gotMatch.first, false);
			if (merged) addedAny = true;
		}
	}
	for (auto breakpoint : breakpoints[rightRead])
	{
		size_t positionWithinMatch;
		if (rightFw)
		{
			if (std::get<0>(breakpoint) < rightStart) continue;
			if (std::get<0>(breakpoint) > rightEnd) continue;
			positionWithinMatch = std::get<0>(breakpoint) - rightStart;
		}
		else
		{
			if (std::get<0>(breakpoint) < readLengths[rightRead]-rightEnd) continue;
			if (std::get<0>(breakpoint) > readLengths[rightRead]-rightStart) continue;
			positionWithinMatch = (readLengths[rightRead]-rightStart) - std::get<0>(breakpoint);
		}
		assert(positionWithinMatch <= rightEnd - rightStart);
		countRightBreakpoints += 1;
		size_t positionInLeft = leftStart + positionWithinMatch;
		assert(positionInLeft <= readLengths[leftRead]);
		auto gotMatch = getBreakpoint(breakpoints[leftRead], leftRead, positionInLeft);
		if (gotMatch.second) addedAny = true;
		bool merged = merge(breakpoints, leftRead, gotMatch.first, rightRead, std::get<0>(breakpoint), rightFw);
		if (merged) addedAny = true;
	}
	if (!(addedAny || countLeftBreakpoints == countRightBreakpoints))
	{
		std::cerr << leftStart << " " << leftEnd << " " << readLengths[leftRead] << std::endl;
		for (auto breakpoint : breakpoints[leftRead])
		{
			std::cerr << std::get<0>(breakpoint) << " ";
		}
		std::cerr << std::endl;
		std::cerr << rightStart << " " << rightEnd << " " << readLengths[rightRead] << std::endl;
		for (auto breakpoint : breakpoints[rightRead])
		{
			std::cerr << std::get<0>(breakpoint) << " ";
		}
		std::cerr << std::endl;
		std::cerr << (rightFw ? "fw" : "bw") << std::endl;
		std::cerr << countLeftBreakpoints << " " << countRightBreakpoints << std::endl;
	}
	assert(addedAny || countLeftBreakpoints == countRightBreakpoints);
	return addedAny;
}

std::vector<std::vector<std::tuple<size_t, size_t, size_t, bool>>> extendAndMergeBreakpoints(const std::vector<size_t>& readLengths, const std::vector<std::tuple<size_t, size_t, size_t, size_t, size_t, size_t, bool>>& matches)
{
	std::vector<std::vector<size_t>> relevantMatches;
	std::vector<std::vector<std::tuple<size_t, size_t, size_t, bool>>> breakpoints;
	relevantMatches.resize(readLengths.size());
	breakpoints.resize(readLengths.size());
	for (size_t i = 0; i < matches.size(); i++)
	{
		relevantMatches[std::get<0>(matches[i])].emplace_back(i);
		relevantMatches[std::get<3>(matches[i])].emplace_back(i);
		breakpoints[std::get<0>(matches[i])].emplace_back(std::get<1>(matches[i]), std::get<0>(matches[i]), std::get<1>(matches[i]), true);
		breakpoints[std::get<0>(matches[i])].emplace_back(std::get<2>(matches[i]), std::get<0>(matches[i]), std::get<2>(matches[i]), true);
		if (std::get<6>(matches[i]))
		{
			breakpoints[std::get<3>(matches[i])].emplace_back(std::get<4>(matches[i]), std::get<3>(matches[i]), std::get<4>(matches[i]), true);
			breakpoints[std::get<3>(matches[i])].emplace_back(std::get<5>(matches[i]), std::get<3>(matches[i]), std::get<5>(matches[i]), true);
		}
		else
		{
			breakpoints[std::get<3>(matches[i])].emplace_back(readLengths[std::get<3>(matches[i])]-std::get<5>(matches[i]), std::get<3>(matches[i]), readLengths[std::get<3>(matches[i])]-std::get<5>(matches[i]), true);
			breakpoints[std::get<3>(matches[i])].emplace_back(readLengths[std::get<3>(matches[i])]-std::get<4>(matches[i]), std::get<3>(matches[i]), readLengths[std::get<3>(matches[i])]-std::get<4>(matches[i]), true);
		}
	}
	std::vector<size_t> checkQueue;
	for (size_t i = 0; i < readLengths.size(); i++)
	{
		checkQueue.emplace_back(i);
		breakpoints[i].emplace_back(0, i, 0, true);
		breakpoints[i].emplace_back(readLengths[i], i, readLengths[i], true);
		std::sort(breakpoints[i].begin(), breakpoints[i].end(), [](auto left, auto right) { return std::get<0>(left) < std::get<0>(right); });
		if (breakpoints[i].size() >= 1)
		{
			for (size_t j = breakpoints[i].size()-1; j > 0; j--)
			{
				if (breakpoints[i][j] != breakpoints[i][j-1]) continue;
				std::swap(breakpoints[i][j], breakpoints[i].back());
				breakpoints[i].pop_back();
			}
		}
		std::sort(breakpoints[i].begin(), breakpoints[i].end(), [](auto left, auto right) { return std::get<0>(left) < std::get<0>(right); });
	}
	while (checkQueue.size() >= 1)
	{
		size_t topRead = checkQueue.back();
		checkQueue.pop_back();
		phmap::flat_hash_set<size_t> newCheckables;
		for (auto matchi : relevantMatches[topRead])
		{
			bool addedAny = extendAndMergeBreakpoints(readLengths, breakpoints, std::get<0>(matches[matchi]), true, std::get<1>(matches[matchi]), std::get<2>(matches[matchi]), std::get<3>(matches[matchi]), std::get<6>(matches[matchi]), std::get<4>(matches[matchi]), std::get<5>(matches[matchi]));
			if (addedAny)
			{
				newCheckables.emplace(std::get<0>(matches[matchi]));
				newCheckables.emplace(std::get<3>(matches[matchi]));
			}
		}
		checkQueue.insert(checkQueue.end(), newCheckables.begin(), newCheckables.end());
	}
	for (size_t i = 0; i < readLengths.size(); i++)
	{
		std::sort(breakpoints[i].begin(), breakpoints[i].end(), [](auto left, auto right) { return std::get<0>(left) < std::get<0>(right); });
		for (size_t j = 0; j < breakpoints[i].size(); j++)
		{
			find(breakpoints, i, std::get<0>(breakpoints[i][j]));
		}
		for (size_t j = 1; j < breakpoints[i].size(); j++)
		{
			assert(std::get<0>(breakpoints[i][j]) > std::get<0>(breakpoints[i][j-1]));
		}
	}
	return breakpoints;
}

void writeGraph(std::string outputFileName, const std::vector<size_t>& breakpointCoverage, const std::vector<size_t>& edgeCoverage, const std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>& edges, const std::vector<size_t>& edgeLength, const size_t minCoverage)
{
	std::ofstream graph { outputFileName };
	size_t countBreakpoints = 0;
	size_t countSegments = 0;
	for (size_t i = 0; i < breakpointCoverage.size(); i++)
	{
		if (breakpointCoverage[i] < minCoverage) continue;
		graph << "S\tbreakpoint_" << i << "\t*\tLN:i:0\tll:f:" << breakpointCoverage[i] << "\tFC:i:" << breakpointCoverage[i] << std::endl;
		countBreakpoints += 1;
	}
	for (size_t i = 0; i < edgeCoverage.size(); i++)
	{
		if (edgeCoverage[i] < minCoverage) continue;
		assert(breakpointCoverage[edges[i].first.first] >= minCoverage);
		assert(breakpointCoverage[edges[i].second.first] >= minCoverage);
		graph << "S\tsegment_" << i << "\t*\tLN:i:" << edgeLength[i] << "\tll:f:" << edgeCoverage[i] << "\tFC:i:" << edgeCoverage[i] * edgeLength[i] << std::endl;
		graph << "L\tbreakpoint_" << edges[i].first.first << "\t" << (edges[i].first.second ? "+" : "-") << "\tsegment_" << i << "\t+\t0M\tec:i:" << edgeCoverage[i] << std::endl;
		graph << "L\tsegment_" << i << "\t+\tbreakpoint_" << edges[i].second.first << "\t" << (edges[i].second.second ? "+" : "-") << "\t0M\tec:i:" << edgeCoverage[i] << std::endl;
		countSegments += 1;
	}
	std::cerr << countBreakpoints << " graph breakpoints" << std::endl;
	std::cerr << countSegments << " graph segments" << std::endl;
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

void mergeSegments(const std::vector<size_t>& readLengths, std::vector<std::vector<std::tuple<size_t, size_t, bool>>>& result, const std::vector<std::vector<std::pair<size_t, size_t>>>& segmentIntervals, const std::vector<std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>>& segmentEndpoints, const size_t leftRead, const bool leftFw, const size_t leftStart, const size_t leftEnd, const size_t rightRead, const bool rightFw, const size_t rightStart, const size_t rightEnd)
{
	assert(leftFw);
	size_t leftFirst = std::numeric_limits<size_t>::max();
	size_t leftLast = std::numeric_limits<size_t>::max();
	for (size_t i = 0; i < segmentIntervals[leftRead].size(); i++)
	{
		if (segmentIntervals[leftRead][i].first == leftStart)
		{
			assert(leftFirst == std::numeric_limits<size_t>::max());
			leftFirst = i;
		}
		if (segmentIntervals[leftRead][i].second == leftEnd)
		{
			assert(leftLast == std::numeric_limits<size_t>::max());
			leftLast = i;
		}
	}
	assert(leftFirst != std::numeric_limits<size_t>::max());
	assert(leftLast != std::numeric_limits<size_t>::max());
	assert(leftLast >= leftFirst);
	size_t rightFirst = std::numeric_limits<size_t>::max();
	size_t rightLast = std::numeric_limits<size_t>::max();
	for (size_t i = 0; i < segmentIntervals[rightRead].size(); i++)
	{
		if (rightFw && segmentIntervals[rightRead][i].first == rightStart)
		{
			assert(rightFirst == std::numeric_limits<size_t>::max());
			rightFirst = i;
		}
		if (rightFw && segmentIntervals[rightRead][i].second == rightEnd)
		{
			assert(rightLast == std::numeric_limits<size_t>::max());
			rightLast = i;
		}
		if (!rightFw && segmentIntervals[rightRead][i].first == readLengths[rightRead]-rightEnd)
		{
			assert(rightFirst == std::numeric_limits<size_t>::max());
			rightFirst = i;
		}
		if (!rightFw && segmentIntervals[rightRead][i].second == readLengths[rightRead]-rightStart)
		{
			assert(rightLast == std::numeric_limits<size_t>::max());
			rightLast = i;
		}
	}
	assert(rightFirst != std::numeric_limits<size_t>::max());
	assert(rightLast != std::numeric_limits<size_t>::max());
	assert(rightLast >= rightFirst);
	if (!(rightLast-rightFirst == leftLast-leftFirst))
	{
		std::cerr << leftStart << " " << leftEnd << std::endl;
		std::cerr << rightStart << " " << rightEnd << std::endl;
		std::cerr << (rightFw ? "fw" : "bw") << std::endl;
		std::cerr << rightFirst << " " << rightLast << std::endl;
		std::cerr << leftFirst << " " << leftLast << std::endl;
	}
	assert(rightLast-rightFirst == leftLast-leftFirst);
	for (size_t i = 0; i <= leftLast-leftFirst; i++)
	{
		if (rightFw)
		{
			assert(segmentIntervals[leftRead][leftFirst+i].second - segmentIntervals[leftRead][leftFirst+i].first == segmentIntervals[rightRead][rightFirst+i].second - segmentIntervals[rightRead][rightFirst+i].first);
			mergeSegments(result, leftRead, leftFirst+i, rightRead, rightFirst+i, true);
		}
		else
		{
			assert(segmentIntervals[leftRead][leftFirst+i].second - segmentIntervals[leftRead][leftFirst+i].first == segmentIntervals[rightRead][rightLast-i].second - segmentIntervals[rightRead][rightLast-i].first);
			mergeSegments(result, leftRead, leftFirst+i, rightRead, rightLast-i, false);
		}
	}
}

std::tuple<std::vector<std::vector<std::tuple<size_t, size_t, bool>>>, std::vector<std::vector<std::pair<size_t, size_t>>>, std::vector<std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>>> mergeSegments(const std::vector<size_t>& readLengths, const std::vector<std::tuple<size_t, size_t, size_t, size_t, size_t, size_t, bool>>& matches, const std::vector<std::vector<std::tuple<size_t, size_t, size_t, bool>>>& breakpoints, const phmap::flat_hash_map<std::pair<size_t, size_t>, size_t>& breakpointToNode)
{
	assert(breakpoints.size() == readLengths.size());
	std::vector<std::vector<std::tuple<size_t, size_t, bool>>> result;
	std::vector<std::vector<std::pair<size_t, size_t>>> segmentIntervals;
	std::vector<std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>> segmentEndpoints;
	result.resize(readLengths.size());
	segmentIntervals.resize(readLengths.size());
	segmentEndpoints.resize(readLengths.size());
	for (size_t i = 0; i < breakpoints.size(); i++)
	{
		size_t countSegments = breakpoints[i].size()-1;
		result[i].resize(countSegments);
		segmentIntervals[i].resize(countSegments);
		segmentEndpoints[i].resize(countSegments);
		for (size_t j = 0; j < countSegments; j++)
		{
			std::get<0>(result[i][j]) = i;
			std::get<1>(result[i][j]) = j;
			std::get<2>(result[i][j]) = true;
		}
		size_t lastPos = std::get<0>(breakpoints[i][0]);
		std::pair<size_t, bool> lastEndpoint { breakpointToNode.at(std::make_pair(std::get<1>(breakpoints[i][0]), std::get<2>(breakpoints[i][0]))), std::get<3>(breakpoints[i][0]) };
		assert(std::get<0>(breakpoints[i][0]) == 0);
		assert(std::get<0>(breakpoints[i].back()) == readLengths[i]);
		for (size_t j = 1; j < breakpoints[i].size(); j++)
		{
			size_t nextPos = std::get<0>(breakpoints[i][j]);
			assert(nextPos != 0);
			assert(j == breakpoints[i].size()-1 || nextPos != readLengths[i]);
			std::pair<size_t, bool> endpoint;
			endpoint.first = breakpointToNode.at(std::make_pair(std::get<1>(breakpoints[i][j]), std::get<2>(breakpoints[i][j])));
			endpoint.second = std::get<3>(breakpoints[i][j]);
			segmentIntervals[i][j-1].first = lastPos;
			segmentIntervals[i][j-1].second = nextPos;
			segmentEndpoints[i][j-1].first = lastEndpoint;
			segmentEndpoints[i][j-1].second = endpoint;
			lastPos = nextPos;
			lastEndpoint = endpoint;
		}
	}
	for (auto match : matches)
	{
		mergeSegments(readLengths, result, segmentIntervals, segmentEndpoints, std::get<0>(match), true, std::get<1>(match), std::get<2>(match), std::get<3>(match), std::get<6>(match), std::get<4>(match), std::get<5>(match));
	}
	for (size_t i = 0; i < result.size(); i++)
	{
		for (size_t j = 0; j < result[i].size(); j++)
		{
			find(result, i, j);
		}
	}
	return std::make_tuple(result, segmentIntervals, segmentEndpoints);
}

phmap::flat_hash_map<std::pair<size_t, size_t>, size_t> getSegmentToEdge(const std::vector<std::vector<std::tuple<size_t, size_t, bool>>>& segments)
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

std::tuple<std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>, std::vector<size_t>, std::vector<size_t>> getEdgeInformation(const std::vector<size_t>& readLengths, const phmap::flat_hash_map<std::pair<size_t, size_t>, size_t>& segmentToEdge, const std::vector<std::vector<std::tuple<size_t, size_t, bool>>>& segments, const std::vector<std::vector<std::pair<size_t, size_t>>>& segmentIntervals, const std::vector<std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>>& segmentEndpoints)
{
	std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>> edges;
	std::vector<size_t> edgeCoverage;
	std::vector<size_t> edgeLength;
	edges.resize(segmentToEdge.size());
	edgeCoverage.resize(segmentToEdge.size(), 0);
	edgeLength.resize(segmentToEdge.size(), std::numeric_limits<size_t>::max());
	for (size_t i = 0; i < segments.size(); i++)
	{
		for (size_t j = 0; j < segments[i].size(); j++)
		{
			std::pair<size_t, bool> edge;
			edge.first = segmentToEdge.at(std::make_pair(std::get<0>(segments[i][j]), std::get<1>(segments[i][j])));
			edge.second = std::get<2>(segments[i][j]);
			if (edge.second)
			{
				edges[edge.first] = segmentEndpoints[i][j];
			}
			else
			{
				edges[edge.first].first = reverse(segmentEndpoints[i][j].second);
				edges[edge.first].second = reverse(segmentEndpoints[i][j].first);
			}
			edgeCoverage[edge.first] += 1;
			assert(edgeLength[edge.first] == std::numeric_limits<size_t>::max() || edgeLength[edge.first] == (segmentIntervals[i][j].second - segmentIntervals[i][j].first));
			edgeLength[edge.first] = (segmentIntervals[i][j].second - segmentIntervals[i][j].first);
		}
	}
	return std::make_tuple(edges, edgeCoverage, edgeLength);
}

void makeGraph(const std::vector<size_t>& readLengths, const std::vector<std::tuple<size_t, size_t, size_t, size_t, size_t, size_t, bool>>& matches, size_t minCoverage, std::string outputFileName)
{
	std::vector<std::vector<std::tuple<size_t, size_t, size_t, bool>>> breakpoints = extendAndMergeBreakpoints(readLengths, matches);
	phmap::flat_hash_map<std::pair<size_t, size_t>, size_t> breakpointToNode = getBreakpointToNode(breakpoints);
	std::vector<size_t> breakpointCoverage = getBreakpointCoverage(breakpointToNode, breakpoints);
	std::vector<std::vector<std::tuple<size_t, size_t, bool>>> segments;
	std::vector<std::vector<std::pair<size_t, size_t>>> segmentIntervals;
	std::vector<std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>> segmentEndpoints;
	std::tie(segments, segmentIntervals, segmentEndpoints) = mergeSegments(readLengths, matches, breakpoints, breakpointToNode);
	phmap::flat_hash_map<std::pair<size_t, size_t>, size_t> segmentToEdge = getSegmentToEdge(segments);
	std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>> edges;
	std::vector<size_t> edgeCoverage;
	std::vector<size_t> edgeLength;
	std::tie(edges, edgeCoverage, edgeLength) = getEdgeInformation(readLengths, segmentToEdge, segments, segmentIntervals, segmentEndpoints);
	writeGraph(outputFileName, breakpointCoverage, edgeCoverage, edges, edgeLength, minCoverage);
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
		auto matches = getKmerMatches(readSequences, std::get<0>(t), std::get<1>(t), std::get<2>(t), std::get<3>(t), std::get<4>(t), std::get<5>(t), std::get<6>(t), 21, 50);
		for (auto match : matches)
		{
			assert(std::get<2>(match) > std::get<1>(match));
			assert(std::get<1>(match) >= std::get<1>(t));
			assert(std::get<2>(match) <= std::get<2>(t));
			assert(std::get<5>(match) > std::get<4>(match));
			assert(std::get<4>(match) >= std::get<4>(t));
			assert(std::get<5>(match) <= std::get<5>(t));
			std::cout << readNames[std::get<0>(match)] << "\tfw\t" << readLengths[std::get<0>(match)] << "\t" << std::get<1>(match) << "\t" << std::get<2>(match) + 20 << "\t" << readNames[std::get<3>(match)] << "\t" << (std::get<6>(match) ? "fw" : "bw") << "\t" << readLengths[std::get<3>(match)] << "\t" << std::get<4>(match) << "\t" << std::get<5>(match) + 20 << std::endl;
		}
		kmermatches.insert(kmermatches.end(), matches.begin(), matches.end());
	}
	std::cerr << kmermatches.size() << " kmer matches" << std::endl;
	makeGraph(readLengths, kmermatches, 2, "graph.gfa");
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
