#include <mutex>
#include <iostream>
#include <vector>
#include <string>
#include <tuple>
#include <set>
#include "ReadStorage.h"
#include "MatchIndex.h"
#include "UnitigResolver.h"
#include "UnitigGraph.h"
#include "UnitigHelper.h"
#include "HPCConsensus.h"

// defined in MBG.cpp but not exposed, so declare them here too
void loadReadsAsHashesMultithread(HashList& result, const size_t kmerSize, const ReadpartIterator& partIterator, const size_t numThreads, std::ostream& log);
UnitigGraph getUnitigGraph(HashList& hashlist, const size_t minCoverage, const double minUnitigCoverage, const bool keepGaps, const bool oneCovHeuristic);
std::vector<ReadPath> getReadPaths(const UnitigGraph& graph, const HashList& hashlist, const size_t numThreads, const ReadpartIterator& partIterator, const size_t kmerSize);

std::vector<std::pair<size_t, bool>> getUniqueReplacementPath(const UnitigGraph& graph, const std::pair<size_t, bool> start, const std::pair<size_t, bool> end, const std::unordered_set<size_t>& pathNodes, const double minCoverage)
{
	std::set<std::pair<size_t, bool>> fwReachable;
	std::vector<std::pair<size_t, bool>> checkStack;
	checkStack.push_back(start);
	while (checkStack.size() > 0)
	{
		auto top = checkStack.back();
		checkStack.pop_back();
		if (fwReachable.count(top) == 1) continue;
		fwReachable.insert(top);
		for (auto edge : graph.edges[top])
		{
			if (fwReachable.count(edge) == 1) continue;
			if (pathNodes.count(edge.first) == 1) continue;
			if (graph.averageCoverage(edge.first) < minCoverage) continue;
			checkStack.push_back(edge);
		}
	}
	std::vector<std::pair<size_t, bool>> result;
	result.push_back(reverse(end));
	std::set<std::pair<size_t, bool>> visited;
	while (result.back() != reverse(start))
	{
		std::pair<size_t, bool> uniqueOption { graph.unitigs.size(), false };
		for (auto edge : graph.edges[result.back()])
		{
			if (fwReachable.count(reverse(edge)) == 0) continue;
			if (graph.averageCoverage(edge.first) < minCoverage) continue;
			if (pathNodes.count(edge.first) == 1 && edge != reverse(start)) continue;
			if (uniqueOption.first == graph.unitigs.size())
			{
				uniqueOption = edge;
			}
			else
			{
				return std::vector<std::pair<size_t, bool>>{};
			}
		}
		if (uniqueOption.first == graph.unitigs.size()) return std::vector<std::pair<size_t, bool>>{};
		assert(visited.count(uniqueOption) == 0);
		visited.insert(uniqueOption);
		result.push_back(uniqueOption);
	}
	assert(result.back() == reverse(start));
	result.pop_back();
	std::reverse(result.begin(), result.end());
	for (auto& node : result)
	{
		node = reverse(node);
	}
	return result;
}

std::vector<Node> getCorrectedPath(const UnitigGraph& graph, const std::vector<Node>& originalPath, const double minSolidCoverage, const double minAmbiguousCoverage)
{
	size_t lastValid = std::numeric_limits<size_t>::max();
	std::vector<Node> result;
	std::unordered_set<size_t> pathNodes;
	for (auto node : originalPath)
	{
		pathNodes.insert(node.id());
	}
	for (size_t i = 0; i < originalPath.size(); i++)
	{
		if (graph.averageCoverage(i) < minSolidCoverage) continue;
		if (lastValid == std::numeric_limits<size_t>::max())
		{
			assert(result.size() == 0);
			result.insert(result.end(), originalPath.begin(), originalPath.begin()+i+1);
			lastValid = i;
			continue;
		}
		if (i == lastValid+1)
		{
			result.push_back(originalPath[i]);
			lastValid = i;
			continue;
		}
		auto uniqueReplacementPath = getUniqueReplacementPath(graph, originalPath[lastValid], originalPath[i], pathNodes, minAmbiguousCoverage);
		if (uniqueReplacementPath.size() > 0)
		{
			bool valid = true;
			for (auto node : uniqueReplacementPath)
			{
				if (graph.averageCoverage(node.first) < minSolidCoverage) valid = false;
			}
			if (valid)
			{
				result.insert(result.end(), uniqueReplacementPath.begin(), uniqueReplacementPath.end());
				lastValid = i;
				continue;
			}
		}
		result.insert(result.end(), originalPath.begin()+lastValid+1, originalPath.begin()+i+1);
		lastValid = i;
		continue;
	}
	if (lastValid < originalPath.size())
	{
		result.insert(result.end(), originalPath.begin()+lastValid+1, originalPath.end());
	}
	return result;
}

bool found(const std::vector<std::pair<size_t, bool>>& vec, Node val)
{
	for (Node v : vec) if (v == val) return true;
	return false;
}

std::string getSequence(const HashList& hashlist, const size_t kmerSize, const std::vector<CompressedSequenceType>& unitigSequences, const StringIndex& stringIndex, const UnitigGraph& graph, const std::vector<Node>& path)
{
	std::string result;
	for (size_t i = 0; i < path.size(); i++)
	{
		std::string add = unitigSequences[path[i].id()].getExpandedSequence(stringIndex);
		if (!path[i].forward())
		{
			add = revCompRaw(add);
		}
		if (i > 0)
		{
			assert(found(graph.edges[path[i-1]], path[i]));
			result += add.substr(getUnitigOverlap(hashlist, kmerSize, graph, path[i-1], path[i]));
		}
		else
		{
			assert(result.size() == 0);
			result = add;
		}
	}
	return result;
}

template <typename F>
void iterateCorrectedSequences(const MatchIndex& matchIndex, const ReadStorage& storage, F callback)
{
	size_t kmerSize = 101;
	size_t windowSize = 50;
	ReadpartIterator partIterator { kmerSize, windowSize, ErrorMasking::No, 1, std::vector<std::string>{}, false, "" };
	std::unordered_map<std::string, size_t> readNameToId;
	storage.iterateReadsFromStorage([&partIterator, &storage, &readNameToId](size_t readid, const std::string& readSequence)
	{
		std::pair<std::string, std::string> nameAndSequence;
		nameAndSequence.first = storage.getNames()[readid];
		readNameToId[nameAndSequence.first] = readid;
		nameAndSequence.second = readSequence;
		partIterator.addMemoryRead(nameAndSequence);
	});
	matchIndex.iterateMatchReadPairs([&partIterator, &storage, kmerSize, callback](size_t read, const std::unordered_set<size_t>& matches)
	{
		std::vector<size_t> useThese;
		useThese.push_back(read);
		for (auto match : matches)
		{
			useThese.push_back(match);
		}
		if (useThese.size() <= 5) return;
		std::string thisReadName = storage.getNames()[read];
		partIterator.setMemoryReadIterables(useThese);
		HashList reads { kmerSize };
		loadReadsAsHashesMultithread(reads, kmerSize, partIterator, 1, std::cerr);
		auto unitigs = getUnitigGraph(reads, 1, 1, false, false);
		auto readPaths = getReadPaths(unitigs, reads, 1, partIterator, kmerSize);
		std::tie(unitigs, readPaths) = resolveUnitigs(unitigs, reads, readPaths, partIterator, 1, kmerSize, 15000, 0, false, false, false, true, false, std::cerr);
		std::vector<Node> path;
		for (auto read : readPaths)
		{
			if (read.readName.first != thisReadName) continue;
			path = read.path;
		}
		if (path.size() == 0) return;
		path = getCorrectedPath(unitigs, path, 5, 2);
		std::vector<CompressedSequenceType> unitigSequences;
		StringIndex stringIndex;
		std::tie(unitigSequences, stringIndex) = getHPCUnitigSequences(reads, unitigs, readPaths, kmerSize, partIterator, 1);
		std::string result = getSequence(reads, kmerSize, unitigSequences, stringIndex, unitigs, path);
		callback(read, result);
	});
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
	std::ofstream correctedOut { "corrected.fa" };
	std::vector<bool> wasCorrected;
	wasCorrected.resize(storage.size(), false);
	iterateCorrectedSequences(matchIndex, storage, [&matchIndex, &storage, &correctedOut, &wasCorrected](size_t readId, const std::string& readSeq)
	{
		correctedOut << ">" << storage.getNames()[readId] << std::endl;
		correctedOut << readSeq << std::endl;
		assert(!wasCorrected[readId]);
		wasCorrected[readId] = true;
	});
	std::ofstream uncorrectedOut { "uncorrected.fa" };
	for (size_t i = 0; i < storage.size(); i++)
	{
		if (wasCorrected[i]) continue;
		auto uncorrected = storage.getRead(i);
		uncorrectedOut << ">" << uncorrected.first << std::endl;
		uncorrectedOut << uncorrected.second << std::endl;
	}
}
