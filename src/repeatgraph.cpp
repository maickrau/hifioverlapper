#include <mutex>
#include <iostream>
#include <vector>
#include <string>
#include "ReadStorage.h"
#include "MatchIndex.h"

// for approx pos finding and add if missing
std::pair<double, bool> getBreakpoint(std::vector<std::tuple<double, size_t, double, bool>>& readBreakpoints, size_t readIndex, double wantedPos)
{
	for (auto pos : readBreakpoints)
	{
		if (std::get<0>(pos) < wantedPos - 0.01) continue;
		if (std::get<0>(pos) > wantedPos + 0.01) continue;
		return std::make_pair(std::get<0>(pos), false);
	}
	readBreakpoints.emplace_back(wantedPos, readIndex, wantedPos, true);
	return std::make_pair(wantedPos, true);
}

// for exact pos finding and never add
size_t findBreakpointIndex(const std::vector<std::vector<std::tuple<double, size_t, double, bool>>>& breakpoints, size_t readIndex, double wantedPos)
{
	std::numeric_limits<size_t>::max();
	size_t found = std::numeric_limits<size_t>::max();
	for (size_t i = 0; i < breakpoints[readIndex].size(); i++)
	{
		if (std::get<0>(breakpoints[readIndex][i]) != wantedPos) continue; // equal without epsilon is safe bc pos is always exact
		assert(found == std::numeric_limits<size_t>::max());
		found = i;
	}
	assert(found != std::numeric_limits<size_t>::max());
	return found;
}

std::tuple<size_t, double, bool> find(std::vector<std::vector<std::tuple<double, size_t, double, bool>>>& breakpoints, size_t readIndex, double wantedPos)
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

void merge(std::vector<std::vector<std::tuple<double, size_t, double, bool>>>& breakpoints, size_t leftRead, double leftPos, size_t rightRead, double rightPos, bool forward)
{
	auto leftParent = find(breakpoints, leftRead, leftPos);
	auto rightParent = find(breakpoints, rightRead, rightPos);
	size_t leftIndex = findBreakpointIndex(breakpoints, std::get<0>(leftParent), std::get<1>(leftParent));
	size_t rightIndex = findBreakpointIndex(breakpoints, std::get<0>(rightParent), std::get<1>(rightParent));
	if (std::get<0>(leftParent) == std::get<0>(rightParent) && std::get<1>(leftParent) == std::get<1>(rightParent))
	{
		assert(std::get<2>(leftParent) ^ std::get<2>(rightParent) ^ forward);
		return;
	}
	assert(std::get<1>(breakpoints[std::get<0>(leftParent)][leftIndex]) == std::get<0>(leftParent));
	assert(std::get<2>(breakpoints[std::get<0>(leftParent)][leftIndex]) == std::get<1>(leftParent));
	assert(std::get<1>(breakpoints[std::get<0>(rightParent)][rightIndex]) == std::get<0>(rightParent));
	assert(std::get<2>(breakpoints[std::get<0>(rightParent)][rightIndex]) == std::get<1>(rightParent));
	std::get<1>(breakpoints[std::get<0>(rightParent)][rightIndex]) = std::get<0>(leftParent);
	std::get<2>(breakpoints[std::get<0>(rightParent)][rightIndex]) = std::get<1>(leftParent);
	std::get<3>(breakpoints[std::get<0>(rightParent)][rightIndex]) = std::get<2>(leftParent) ^ std::get<2>(rightParent) ^ forward;
}

bool extendAndMergeBreakpoints(const std::vector<size_t>& readLengths, std::vector<std::vector<std::tuple<double, size_t, double, bool>>>& breakpoints, size_t leftRead, bool leftFw, size_t leftStart, size_t leftEnd, size_t rightRead, bool rightFw, size_t rightStart, size_t rightEnd)
{
	assert(leftFw);
	assert(leftRead != rightRead);
	bool addedAny = false;
	double leftStartFraction = (double)leftStart / (double)readLengths[leftRead];
	double leftEndFraction = (double)leftEnd / (double)readLengths[leftRead];
	double rightStartFraction = (double)rightStart / (double)readLengths[rightRead];
	double rightEndFraction = (double)rightEnd / (double)readLengths[rightRead];
	for (auto breakpoint : breakpoints[leftRead])
	{
		if (std::get<0>(breakpoint) < leftStartFraction) continue; // comparison without epsilon is safe but only bc fraction is computed same way in matches and here
		if (std::get<0>(breakpoint) > leftEndFraction) continue;
		double positionWithinMatch = (std::get<0>(breakpoint) - leftStartFraction) / (leftEndFraction - leftStartFraction);
		double positionInRight = rightStartFraction + positionWithinMatch * (rightEndFraction - rightStartFraction);
		if (rightFw)
		{
			auto gotMatch = getBreakpoint(breakpoints[rightRead], rightRead, positionInRight);
			if (gotMatch.second) addedAny = true;
			merge(breakpoints, leftRead, std::get<0>(breakpoint), rightRead, gotMatch.first, true);
		}
		else
		{
			positionInRight = 1.0 - positionInRight;
			auto gotMatch = getBreakpoint(breakpoints[rightRead], rightRead, positionInRight);
			if (gotMatch.second) addedAny = true;
			merge(breakpoints, leftRead, std::get<0>(breakpoint), rightRead, gotMatch.first, false);
		}
	}
	for (auto breakpoint : breakpoints[rightRead])
	{
		double positionWithinMatch;
		if (rightFw)
		{
			if (std::get<0>(breakpoint) < rightStartFraction) continue; // comparison without epsilon is safe but only bc fraction is computed same way in matches and here
			if (std::get<0>(breakpoint) > rightEndFraction) continue;
			positionWithinMatch = (std::get<0>(breakpoint) - rightStartFraction) / (rightEndFraction - rightStartFraction);
		}
		else
		{
			if (std::get<0>(breakpoint) < 1.0-rightStartFraction) continue; // comparison without epsilon is safe but only bc fraction is computed same way in matches and here
			if (std::get<0>(breakpoint) > 1.0-rightEndFraction) continue;
			positionWithinMatch = (rightEndFraction - std::get<0>(breakpoint)) / (rightEndFraction - rightStartFraction);
		}
		double positionInLeft = leftStartFraction + positionWithinMatch * (leftEndFraction - leftStartFraction);
		auto gotMatch = getBreakpoint(breakpoints[leftRead], leftRead, positionInLeft);
		if (gotMatch.second) addedAny = true;
		merge(breakpoints, leftRead, gotMatch.first, rightRead, std::get<0>(breakpoint), rightFw);
	}
	return addedAny;
}

std::vector<std::vector<std::tuple<double, size_t, double, bool>>> extendAndMergeBreakpoints(const std::vector<size_t>& readLengths, const std::vector<std::tuple<size_t, size_t, size_t, size_t, size_t, size_t, bool>>& matches)
{
	std::vector<std::vector<size_t>> relevantMatches;
	std::vector<std::vector<std::tuple<double, size_t, double, bool>>> breakpoints;
	relevantMatches.resize(readLengths.size());
	breakpoints.resize(readLengths.size());
	for (size_t i = 0; i < matches.size(); i++)
	{
		relevantMatches[std::get<0>(matches[i])].emplace_back(i);
		relevantMatches[std::get<3>(matches[i])].emplace_back(i);
		breakpoints[std::get<0>(matches[i])].emplace_back((double)std::get<1>(matches[i])/(double)readLengths[std::get<0>(matches[i])], std::get<0>(matches[i]), (double)std::get<1>(matches[i])/(double)readLengths[std::get<0>(matches[i])], true);
		breakpoints[std::get<0>(matches[i])].emplace_back((double)std::get<2>(matches[i])/(double)readLengths[std::get<0>(matches[i])], std::get<0>(matches[i]), (double)std::get<2>(matches[i])/(double)readLengths[std::get<0>(matches[i])], true);
		if (std::get<6>(matches[i]))
		{
			breakpoints[std::get<3>(matches[i])].emplace_back((double)std::get<4>(matches[i])/(double)readLengths[std::get<3>(matches[i])], std::get<3>(matches[i]), (double)std::get<4>(matches[i])/(double)readLengths[std::get<3>(matches[i])], true);
			breakpoints[std::get<3>(matches[i])].emplace_back((double)std::get<5>(matches[i])/(double)readLengths[std::get<3>(matches[i])], std::get<3>(matches[i]), (double)std::get<5>(matches[i])/(double)readLengths[std::get<3>(matches[i])], true);
		}
		else
		{
			breakpoints[std::get<3>(matches[i])].emplace_back(1.0-(double)std::get<5>(matches[i])/(double)readLengths[std::get<3>(matches[i])], std::get<3>(matches[i]), 1.0-(double)std::get<5>(matches[i])/(double)readLengths[std::get<3>(matches[i])], true);
			breakpoints[std::get<3>(matches[i])].emplace_back(1.0-(double)std::get<4>(matches[i])/(double)readLengths[std::get<3>(matches[i])], std::get<3>(matches[i]), 1.0-(double)std::get<4>(matches[i])/(double)readLengths[std::get<3>(matches[i])], true);
		}
	}
	std::vector<size_t> checkQueue;
	for (size_t i = 0; i < readLengths.size(); i++)
	{
		checkQueue.emplace_back(i);
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

phmap::flat_hash_map<std::pair<size_t, double>, size_t> getBreakpointToNode(const std::vector<std::vector<std::tuple<double, size_t, double, bool>>>& breakpoints)
{
	size_t nextIndex = 0;
	phmap::flat_hash_map<std::pair<size_t, double>, size_t> result;
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

std::vector<size_t> getBreakpointCoverage(const phmap::flat_hash_map<std::pair<size_t, double>, size_t>& breakpointToNode, const std::vector<std::vector<std::tuple<double, size_t, double, bool>>>& breakpoints)
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

void mergeSegments(std::vector<std::vector<std::tuple<size_t, size_t, bool>>>& result, const std::vector<std::vector<std::pair<double, double>>>& segmentIntervals, const std::vector<std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>>& segmentEndpoints, const size_t leftRead, const bool leftFw, const double leftStart, const double leftEnd, const size_t rightRead, const bool rightFw, const double rightStart, const double rightEnd)
{
	assert(leftFw);
	size_t leftIndex = std::numeric_limits<size_t>::max();
	for (size_t i = 0; i < segmentIntervals[leftRead].size(); i++)
	{
		if (segmentIntervals[leftRead][i].second > leftStart)
		{
			leftIndex = i;
			break;
		}
	}
	assert(leftIndex != std::numeric_limits<size_t>::max());
	if (rightFw)
	{
		size_t rightIndex = std::numeric_limits<size_t>::max();
		for (size_t i = 0; i < segmentIntervals[rightRead].size(); i++)
		{
			if (segmentIntervals[rightRead][i].second > rightStart)
			{
				rightIndex = i;
				break;
			}
		}
		assert(rightIndex != std::numeric_limits<size_t>::max());
		while (leftIndex < segmentIntervals[leftRead].size() && rightIndex < segmentIntervals[rightRead].size() && segmentIntervals[leftRead][leftIndex].first <= leftEnd && segmentIntervals[rightRead][rightIndex].first <= rightEnd)
		{
			double leftStartWithinAlignment = (segmentIntervals[leftRead][leftIndex].first - leftStart) / (leftEnd - leftStart);
			double leftEndWithinAlignment = (segmentIntervals[leftRead][leftIndex].second - leftStart) / (leftEnd - leftStart);
			double rightStartWithinAlignment = (segmentIntervals[rightRead][rightIndex].first - rightStart) / (rightEnd - rightStart);
			double rightEndWithinAlignment = (segmentIntervals[rightRead][rightIndex].second - rightStart) / (rightEnd - rightStart);
			if (leftEndWithinAlignment > rightStartWithinAlignment && leftStartWithinAlignment < rightEndWithinAlignment)
			{
				if (segmentEndpoints[leftRead][leftIndex] == segmentEndpoints[rightRead][rightIndex] && segmentEndpoints[leftRead][leftIndex].first.first != std::numeric_limits<size_t>::max() && segmentEndpoints[leftRead][leftIndex].second.first != std::numeric_limits<size_t>::max())
				{
					mergeSegments(result, leftRead, leftIndex, rightRead, rightIndex, true);
				}
			}
			if (segmentIntervals[leftRead][leftIndex].second > segmentIntervals[rightRead][rightIndex].second)
			{
				rightIndex += 1;
			}
			else if (segmentIntervals[rightRead][rightIndex].second > segmentIntervals[leftRead][leftIndex].second)
			{
				leftIndex += 1;
			}
			else
			{
				leftIndex += 1;
				rightIndex += 1;
			}
		}
	}
	else
	{
		size_t rightIndex = std::numeric_limits<size_t>::max();
		for (size_t i = segmentIntervals[rightRead].size()-1; i < segmentIntervals[rightRead].size(); i--)
		{
			if (segmentIntervals[rightRead][i].first < 1.0-rightStart)
			{
				rightIndex = i;
				break;
			}
		}
		assert(rightIndex != std::numeric_limits<size_t>::max());
		while (leftIndex < segmentIntervals[leftRead].size() && rightIndex < segmentIntervals[rightRead].size() && segmentIntervals[leftRead][leftIndex].first <= leftEnd && segmentIntervals[rightRead][rightIndex].second >= 1.0-rightEnd)
		{
			double leftStartWithinAlignment = (segmentIntervals[leftRead][leftIndex].first - leftStart) / (leftEnd - leftStart);
			double leftEndWithinAlignment = (segmentIntervals[leftRead][leftIndex].second - leftStart) / (leftEnd - leftStart);
			double rightStartWithinAlignment = ((1.0-rightStart) - segmentIntervals[rightRead][rightIndex].second) / (rightEnd - rightStart);
			double rightEndWithinAlignment = ((1.0-rightStart) - segmentIntervals[rightRead][rightIndex].first) / (rightEnd - rightStart);
			if (leftEndWithinAlignment > rightStartWithinAlignment && leftStartWithinAlignment < rightEndWithinAlignment)
			{
				if (segmentEndpoints[leftRead][leftIndex].first == reverse(segmentEndpoints[rightRead][rightIndex].second) && segmentEndpoints[leftRead][leftIndex].second == reverse(segmentEndpoints[rightRead][rightIndex].first) && segmentEndpoints[leftRead][leftIndex].first.first != std::numeric_limits<size_t>::max() && segmentEndpoints[leftRead][leftIndex].second.first != std::numeric_limits<size_t>::max())
				{
					mergeSegments(result, leftRead, leftIndex, rightRead, rightIndex, false);
				}
			}
			if (segmentIntervals[leftRead][leftIndex].second > 1.0-segmentIntervals[rightRead][rightIndex].first)
			{
				rightIndex -= 1;
			}
			else if (1.0-segmentIntervals[rightRead][rightIndex].first > segmentIntervals[leftRead][leftIndex].second)
			{
				leftIndex += 1;
			}
			else
			{
				leftIndex += 1;
				rightIndex -= 1;
			}
		}
	}
}

std::tuple<std::vector<std::vector<std::tuple<size_t, size_t, bool>>>, std::vector<std::vector<std::pair<double, double>>>, std::vector<std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>>> mergeSegments(const std::vector<size_t>& readLengths, const std::vector<std::tuple<size_t, size_t, size_t, size_t, size_t, size_t, bool>>& matches, const std::vector<std::vector<std::tuple<double, size_t, double, bool>>>& breakpoints, const phmap::flat_hash_map<std::pair<size_t, double>, size_t>& breakpointToNode)
{
	assert(breakpoints.size() == readLengths.size());
	std::vector<std::vector<std::tuple<size_t, size_t, bool>>> result;
	std::vector<std::vector<std::pair<double, double>>> segmentIntervals;
	std::vector<std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>> segmentEndpoints;
	result.resize(readLengths.size());
	segmentIntervals.resize(readLengths.size());
	segmentEndpoints.resize(readLengths.size());
	for (size_t i = 0; i < breakpoints.size(); i++)
	{
		size_t countSegments = breakpoints[i].size()+1;
		result[i].resize(countSegments);
		segmentIntervals[i].resize(countSegments);
		segmentEndpoints[i].resize(countSegments);
		for (size_t j = 0; j < countSegments; j++)
		{
			std::get<0>(result[i][j]) = i;
			std::get<1>(result[i][j]) = j;
			std::get<2>(result[i][j]) = true;
		}
		double lastPos = 0;
		std::pair<size_t, bool> lastEndpoint { std::numeric_limits<size_t>::max(), true };
		for (size_t j = 0; j < breakpoints[i].size(); j++)
		{
			double nextPos = std::get<0>(breakpoints[i][j]);
			std::pair<size_t, bool> endpoint;
			endpoint.first = breakpointToNode.at(std::make_pair(std::get<1>(breakpoints[i][j]), std::get<2>(breakpoints[i][j])));
			endpoint.second = std::get<3>(breakpoints[i][j]);
			segmentIntervals[i][j].first = lastPos;
			segmentIntervals[i][j].second = nextPos;
			segmentEndpoints[i][j].first = lastEndpoint;
			segmentEndpoints[i][j].second = endpoint;
			lastPos = nextPos;
			lastEndpoint = endpoint;
		}
		segmentIntervals[i].back().first = lastPos;
		segmentIntervals[i].back().second = 1.0;
		segmentEndpoints[i].back().first = lastEndpoint;
		segmentEndpoints[i].back().second = std::make_pair(std::numeric_limits<size_t>::max(), true);
	}
	for (auto match : matches)
	{
		double leftStartFraction = (double)std::get<1>(match)/(double)readLengths[std::get<0>(match)];
		double leftEndFraction = (double)std::get<2>(match)/(double)readLengths[std::get<0>(match)];
		double rightStartFraction = (double)std::get<4>(match)/(double)readLengths[std::get<3>(match)];
		double rightEndFraction = (double)std::get<5>(match)/(double)readLengths[std::get<3>(match)];
		mergeSegments(result, segmentIntervals, segmentEndpoints, std::get<0>(match), true, leftStartFraction, leftEndFraction, std::get<3>(match), std::get<6>(match), rightStartFraction, rightEndFraction);
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

std::tuple<std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>, std::vector<size_t>, std::vector<size_t>> getEdgeInformation(const std::vector<size_t>& readLengths, const phmap::flat_hash_map<std::pair<size_t, size_t>, size_t>& segmentToEdge, const std::vector<std::vector<std::tuple<size_t, size_t, bool>>>& segments, const std::vector<std::vector<std::pair<double, double>>>& segmentIntervals, const std::vector<std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>>& segmentEndpoints)
{
	std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>> edges;
	std::vector<size_t> edgeCoverage;
	std::vector<size_t> edgeLength;
	edges.resize(segmentToEdge.size());
	edgeCoverage.resize(segmentToEdge.size(), 0);
	edgeLength.resize(segmentToEdge.size(), 0);
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
			edgeLength[edge.first] += (segmentIntervals[i][j].second - segmentIntervals[i][j].first) * readLengths[i];
		}
	}
	for (size_t i = 0; i < edgeLength.size(); i++)
	{
		edgeLength[i] /= edgeCoverage[i];
	}
	return std::make_tuple(edges, edgeCoverage, edgeLength);
}

void makeGraph(const std::vector<size_t>& readLengths, const std::vector<std::tuple<size_t, size_t, size_t, size_t, size_t, size_t, bool>>& matches, size_t minCoverage, std::string outputFileName)
{
	std::vector<std::vector<std::tuple<double, size_t, double, bool>>> breakpoints = extendAndMergeBreakpoints(readLengths, matches);
	phmap::flat_hash_map<std::pair<size_t, double>, size_t> breakpointToNode = getBreakpointToNode(breakpoints);
	std::vector<size_t> breakpointCoverage = getBreakpointCoverage(breakpointToNode, breakpoints);
	std::vector<std::vector<std::tuple<size_t, size_t, bool>>> segments;
	std::vector<std::vector<std::pair<double, double>>> segmentIntervals;
	std::vector<std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>> segmentEndpoints;
	std::tie(segments, segmentIntervals, segmentEndpoints) = mergeSegments(readLengths, matches, breakpoints, breakpointToNode);
	phmap::flat_hash_map<std::pair<size_t, size_t>, size_t> segmentToEdge = getSegmentToEdge(segments);
	std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>> edges;
	std::vector<size_t> edgeCoverage;
	std::vector<size_t> edgeLength;
	std::tie(edges, edgeCoverage, edgeLength) = getEdgeInformation(readLengths, segmentToEdge, segments, segmentIntervals, segmentEndpoints);
	writeGraph(outputFileName, breakpointCoverage, edgeCoverage, edges, edgeLength, minCoverage);
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
	for (auto file : readFiles)
	{
		std::mutex indexMutex;
		storage.iterateReadsFromFile(file, numThreads, false, [&matchIndex, &indexMutex](size_t readName, const std::string& sequence)
		{
			matchIndex.addMatchesFromRead(readName, indexMutex, sequence);
		});
	}
	size_t numWindowChunks = matchIndex.numWindowChunks();
	size_t numUniqueChunks = matchIndex.numUniqueChunks();
	matchIndex.clearConstructionVariablesAndCompact();
	const std::vector<size_t>& readLengths = storage.getRawReadLengths();
	std::mutex printMutex;
	std::vector<std::tuple<size_t, size_t, size_t, size_t, size_t, size_t, bool>> matches;
	auto result = matchIndex.iterateMatchChains(numThreads, storage.getRawReadLengths(), [&printMutex, &matches, minAlignmentLength, &readLengths](const size_t left, const size_t leftstart, const size_t leftend, const bool leftFw, const size_t right, const size_t rightstart, const size_t rightend, const bool rightFw)
	{
		if (leftend+1-leftstart < minAlignmentLength) return;
		if (rightend+1-rightstart < minAlignmentLength) return;
		std::lock_guard<std::mutex> lock { printMutex };
		assert(leftFw);
		matches.emplace_back(left, leftstart, leftend+1, right, rightstart, rightend+1, rightFw); // param coordinates are inclusive, switch to exclusive
		//matches.emplace_back(left, readLengths[left]-leftend-1, readLengths[left]-leftstart, right, readLengths[right]-rightend-1, readLengths[right]-rightstart, !rightFw); // param coordinates are inclusive, switch to exclusive
	});
	std::sort(matches.begin(), matches.end(), [](auto left, auto right) { return std::get<2>(left)-std::get<1>(left) + std::get<5>(left)-std::get<4>(left) > std::get<2>(right)-std::get<1>(right) + std::get<5>(right)-std::get<4>(right); });
	makeGraph(readLengths, matches, 2, "graph.gfa");
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
