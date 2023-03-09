#include <mutex>
#include <iostream>
#include <vector>
#include <string>
#include <tuple>
#include "ReadStorage.h"
#include "MatchIndex.h"
#include "UnitigResolver.h"
#include "UnitigGraph.h"

// defined in MBG.cpp but not exposed, so declare them here too
void loadReadsAsHashesMultithread(HashList& result, const size_t kmerSize, const ReadpartIterator& partIterator, const size_t numThreads, std::ostream& log);
UnitigGraph getUnitigGraph(HashList& hashlist, const size_t minCoverage, const double minUnitigCoverage, const bool keepGaps, const bool oneCovHeuristic);
std::vector<ReadPath> getReadPaths(const UnitigGraph& graph, const HashList& hashlist, const size_t numThreads, const ReadpartIterator& partIterator, const size_t kmerSize);

std::pair<size_t, bool> findBubble(const UnitigGraph& unitigs, std::pair<size_t, bool> start)
{
	if (unitigs.edges[start].size() != 2) return std::make_pair(std::numeric_limits<size_t>::max(), false);
	std::pair<size_t, bool> allele1;
	std::pair<size_t, bool> allele2;
	allele1 = *unitigs.edges[start].begin();
	allele2 = *(++unitigs.edges[start].begin());
	if (allele1.first == allele2.first) return std::make_pair(std::numeric_limits<size_t>::max(), false);
	if (unitigs.edges[allele1].size() != 1) return std::make_pair(std::numeric_limits<size_t>::max(), false);
	if (unitigs.edges[allele2].size() != 1) return std::make_pair(std::numeric_limits<size_t>::max(), false);
	if (unitigs.edges[reverse(allele1)].size() != 1) return std::make_pair(std::numeric_limits<size_t>::max(), false);
	if (unitigs.edges[reverse(allele2)].size() != 1) return std::make_pair(std::numeric_limits<size_t>::max(), false);
	std::pair<size_t, bool> end = *(unitigs.edges[allele1].begin());
	if (end != *unitigs.edges[allele2].begin()) return std::make_pair(std::numeric_limits<size_t>::max(), false);
	if (unitigs.edges[reverse(end)].size() != 2) return std::make_pair(std::numeric_limits<size_t>::max(), false);
	return end;
}

std::vector<std::pair<std::vector<size_t>, std::vector<size_t>>> getBubbleInfos(const UnitigGraph& unitigs, const std::vector<ReadPath>& readPaths, const std::unordered_map<std::string, size_t>& readNameToId)
{
	std::vector<std::pair<std::vector<size_t>, std::vector<size_t>>> result;
	std::unordered_set<std::pair<size_t, bool>> processed;
	for (size_t i = 0; i < unitigs.unitigs.size(); i++)
	{
		std::pair<size_t, bool> fw { i, true };
		if (processed.count(fw) == 0)
		{
			auto foundBubble = findBubble(unitigs, fw);
			if (foundBubble.first != std::numeric_limits<size_t>::max())
			{
				processed.insert(reverse(foundBubble));
				std::pair<size_t, bool> allele1;
				std::pair<size_t, bool> allele2;
				assert(unitigs.edges[fw].size() == 2);
				allele1 = *unitigs.edges[fw].begin();
				allele2 = *(++unitigs.edges[fw].begin());
				assert(allele1.first != allele2.first);
				result.emplace_back();
				for (const auto& path : readPaths)
				{
					for (const auto& node : path.path)
					{
						if (node.id() == allele1.first)
						{
							result.back().first.push_back(readNameToId.at(path.readName.first));
						}
						else if (node.id() == allele2.first)
						{
							result.back().second.push_back(readNameToId.at(path.readName.first));
						}
					}
				}
				std::sort(result.back().first.begin(), result.back().first.end());
				std::sort(result.back().second.begin(), result.back().second.end());
			}
		}
		std::pair<size_t, bool> bw { i, false };
		if (processed.count(bw) == 0)
		{
			auto foundBubble = findBubble(unitigs, bw);
			if (foundBubble.first != std::numeric_limits<size_t>::max())
			{
				processed.insert(reverse(foundBubble));
				std::pair<size_t, bool> allele1;
				std::pair<size_t, bool> allele2;
				assert(unitigs.edges[bw].size() == 2);
				allele1 = *unitigs.edges[bw].begin();
				allele2 = *(++unitigs.edges[bw].begin());
				assert(allele1.first != allele2.first);
				result.emplace_back();
				for (const auto& path : readPaths)
				{
					for (const auto& node : path.path)
					{
						if (node.id() == allele1.first)
						{
							result.back().first.push_back(readNameToId.at(path.readName.first));
						}
						else if (node.id() == allele2.first)
						{
							result.back().second.push_back(readNameToId.at(path.readName.first));
						}
					}
				}
				std::sort(result.back().first.begin(), result.back().first.end());
				std::sort(result.back().second.begin(), result.back().second.end());
			}
		}
	}
	return result;
}

size_t getNumConnections(const std::vector<size_t>& left, const std::vector<size_t>& right)
{
	size_t lefti = 0;
	size_t righti = 0;
	size_t result = 0;
	while (lefti < left.size() && righti < right.size())
	{
		if (left[lefti] == right[righti])
		{
			result += 1;
			lefti += 1;
			righti += 1;
		}
		else if (left[lefti] > right[righti])
		{
			righti += 1;
		}
		else if (right[righti] > left[lefti])
		{
			lefti += 1;
		}
	}
	return result;
}

bool isValidPartition(const std::pair<std::vector<size_t>, std::vector<size_t>>& left, const std::pair<std::vector<size_t>, std::vector<size_t>>& right, size_t minCoverage)
{
	size_t cisConnections = getNumConnections(left.first, right.first) + getNumConnections(left.second, right.second);
	size_t transConnections = getNumConnections(left.first, right.second) + getNumConnections(left.second, right.first);
	if (transConnections == 0 && cisConnections >= minCoverage) return true;
	if (cisConnections == 0 && transConnections >= minCoverage) return true;
	return false;
}

size_t canonnum(size_t left, size_t right)
{
	if (left >= right) std::swap(left, right);
	assert(left < std::numeric_limits<uint32_t>::max());
	assert(right < std::numeric_limits<uint32_t>::max());
	return (left << 32ull) + right;
}

std::unordered_set<uint64_t> forbidDifferentHaplos(const MatchIndex& matchIndex, const ReadStorage& storage, const std::unordered_set<uint64_t>& initialForbidden, const ErrorMasking errorMasking)
{
	std::unordered_set<uint64_t> result;
	size_t kmerSize = 501;
	size_t windowSize = 250;
	ReadpartIterator partIterator { kmerSize, windowSize, errorMasking, 1, std::vector<std::string>{}, false, "" };
	std::unordered_map<std::string, size_t> readNameToId;
	storage.iterateReadsFromStorage([&partIterator, &storage, &readNameToId](size_t readid, const std::string& readSequence)
	{
		std::pair<std::string, std::string> nameAndSequence;
		nameAndSequence.first = storage.getNames()[readid];
		readNameToId[nameAndSequence.first] = readid;
		nameAndSequence.second = readSequence;
		partIterator.addMemoryRead(nameAndSequence);
	});
	matchIndex.iterateMatchReadPairs([&initialForbidden, &result, &partIterator, &readNameToId, kmerSize, windowSize](size_t read, const std::unordered_set<size_t>& matches)
	{
		std::vector<size_t> useThese;
		useThese.push_back(read);
		for (auto match : matches)
		{
			if (initialForbidden.count(canonnum(read, match)) == 1) continue;
			useThese.push_back(match);
		}
		if (useThese.size() <= 5) return;
		partIterator.setMemoryReadIterables(useThese);
		HashList reads { kmerSize };
		loadReadsAsHashesMultithread(reads, kmerSize, partIterator, 1, std::cerr);
		auto unitigs = getUnitigGraph(reads, 1, 1, false, false);
		auto readPaths = getReadPaths(unitigs, reads, 1, partIterator, kmerSize);
		std::tie(unitigs, readPaths) = resolveUnitigs(unitigs, reads, readPaths, partIterator, 1, kmerSize, 15000, 0, false, false, false, true, false, std::cerr);
		auto bubbleInfo = getBubbleInfos(unitigs, readPaths, readNameToId);
		for (size_t i = 0; i < bubbleInfo.size(); i++)
		{
			for (size_t j = i+1; j < bubbleInfo.size(); j++)
			{
				if (isValidPartition(bubbleInfo[i], bubbleInfo[j], 5))
				{
					for (size_t leftread : bubbleInfo[i].first)
					{
						for (size_t rightread : bubbleInfo[i].second)
						{
							result.insert(canonnum(leftread, rightread));
						}
					}
					for (size_t leftread : bubbleInfo[j].first)
					{
						for (size_t rightread : bubbleInfo[j].second)
						{
							result.insert(canonnum(leftread, rightread));
						}
					}
				}
			}
		}
	});
	return result;
}

int main(int argc, char** argv)
{
	size_t numThreads = std::stoi(argv[1]);
	size_t k = std::stoull(argv[2]);
	size_t numWindows = std::stoull(argv[3]);
	size_t windowSize = std::stoull(argv[4]);
	std::vector<std::string> readFiles;
	for (size_t i = 5; i < argc; i++)
	{
		readFiles.emplace_back(argv[i]);
	}
	MatchIndex matchIndex { k, numWindows, windowSize };
	ReadStorage storage;
	for (auto file : readFiles)
	{
		std::mutex indexMutex;
		storage.iterateReadsFromFile(file, numThreads, true, [&matchIndex, &indexMutex](size_t readName, const std::string& sequence)
		{
			matchIndex.addMatchesFromRead(readName, indexMutex, sequence);
		});
	}
	std::unordered_set<uint64_t> haploForbiddenMatches = forbidDifferentHaplos(matchIndex, storage, std::unordered_set<uint64_t>{}, ErrorMasking::Microsatellite);
	std::unordered_set<uint64_t> haploForbiddenMatchesHpc = forbidDifferentHaplos(matchIndex, storage, haploForbiddenMatches, ErrorMasking::Hpc);
	haploForbiddenMatches.insert(haploForbiddenMatchesHpc.begin(), haploForbiddenMatchesHpc.end());
	std::unordered_set<uint64_t> haploForbiddenMatchesNomask = forbidDifferentHaplos(matchIndex, storage, haploForbiddenMatches, ErrorMasking::No);
	haploForbiddenMatches.insert(haploForbiddenMatchesNomask.begin(), haploForbiddenMatchesNomask.end());
	matchIndex.iterateMatchReadPairs([&haploForbiddenMatches, &storage](size_t read, const std::unordered_set<size_t>& matches)
	{
		for (size_t match : matches)
		{
			if (haploForbiddenMatches.count(canonnum(read, match)) == 1) continue;
			std::cout << storage.getNames()[read] << "\t" << storage.getNames()[match] << std::endl;
		}
	});
	std::cerr << haploForbiddenMatches.size() << " read-pairs forbidden by haplotype analysis" << std::endl;
}
