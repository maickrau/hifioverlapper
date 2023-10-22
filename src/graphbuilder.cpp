#include <mutex>
#include <iostream>
#include <vector>
#include <string>
#include "ReadStorage.h"
#include "MatchIndex.h"

const size_t minCoverage = 3;

void getUnitig(std::vector<size_t>& belongsToUnitig, const std::vector<std::vector<size_t>>& outEdges, const std::vector<std::vector<size_t>>& inEdges, size_t start, std::vector<size_t>& unitigEnd)
{
	size_t realStart = start;
	assert(belongsToUnitig[start] == std::numeric_limits<size_t>::max());
	while (inEdges[start].size() == 1)
	{
		size_t next = inEdges[start][0];
		if (outEdges[next].size() != 1) break;
		if (next == start) break;
		start = next;
		if (start == realStart) break;
		assert(belongsToUnitig[start] == std::numeric_limits<size_t>::max());
	}
	size_t unitigStart = start;
	assert(belongsToUnitig[start] == std::numeric_limits<size_t>::max());
	belongsToUnitig[start] = unitigEnd.size();
	while (outEdges[start].size() == 1)
	{
		size_t next = outEdges[start][0];
		if (inEdges[next].size() != 1) break;
		if (next == start) break;
		start = next;
		if (start == unitigStart) break;
		assert(belongsToUnitig[start] == std::numeric_limits<size_t>::max());
		belongsToUnitig[start] = unitigEnd.size();
	}
	unitigEnd.push_back(start);
}

void removeArtefactEdges(std::vector<std::vector<size_t>>& outEdges, std::vector<std::vector<size_t>>& inEdges)
{
	std::vector<std::pair<size_t, size_t>> removeEdges;
	for (size_t i = 0; i < outEdges.size(); i++)
	{
		if (outEdges[i].size() != 2) continue;
		if (outEdges[outEdges[i][0]].size() == 1 && outEdges[outEdges[i][0]][0] == outEdges[i][1])
		{
			removeEdges.emplace_back(i, outEdges[i][0]);
		}
		if (outEdges[outEdges[i][1]].size() == 1 && outEdges[outEdges[i][1]][0] == outEdges[i][0])
		{
			removeEdges.emplace_back(i, outEdges[i][1]);
		}
	}
	if (removeEdges.size() == 0) return;
	for (auto edge : removeEdges)
	{
		for (size_t i = outEdges[edge.first].size()-1; i < outEdges[edge.first].size(); i--)
		{
			if (outEdges[edge.first][i] == edge.second)
			{
				std::swap(outEdges[edge.first][i], outEdges[edge.first].back());
				outEdges[edge.first].pop_back();
			}
		}
	}
	inEdges.clear();
	inEdges.resize(outEdges.size());
	for (size_t i = 0; i < outEdges.size(); i++)
	{
		for (size_t j = 0; j < outEdges[i].size(); j++)
		{
			inEdges[outEdges[i][j]].push_back(i);
		}
	}
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
		storage.iterateReadsFromFile(file, numThreads, false, [&matchIndex, &indexMutex](size_t readName, const std::string& sequence)
		{
			matchIndex.addMatchesFromRead(readName, indexMutex, sequence);
		});
	}
	std::vector<size_t> merCoverage;
	merCoverage.resize(matchIndex.multinumberSize());
	for (size_t i = 0; i < merCoverage.size(); i++)
	{
		merCoverage[i] = matchIndex.multinumberReadCount(i);
		if (merCoverage[i] < minCoverage) merCoverage[i] = 0;
	}
	std::vector<std::vector<size_t>> readMerMatches;
	readMerMatches.resize(storage.getNames().size());
	matchIndex.iterateReadMultinumberMatches(minCoverage, [&readMerMatches](const size_t mer, const size_t read, const size_t readStart, const size_t readEnd)
	{
		assert(read < readMerMatches.size());
		assert(mer < (size_t)std::numeric_limits<uint32_t>::max());
		readMerMatches[read].emplace_back((readStart << 32) + mer);
	});
	std::vector<size_t> edges;
	for (size_t i = 0; i < readMerMatches.size(); i++)
	{
		std::sort(readMerMatches[i].begin(), readMerMatches[i].end());
		for (size_t j = 1; j < readMerMatches[i].size(); j++)
		{
			size_t prevPos = readMerMatches[i][j-1] >> 32;
			size_t thisPos = readMerMatches[i][j] >> 32;
			if ((prevPos & 0x80000000) != (thisPos & 0x80000000)) continue;
			size_t prevMer = readMerMatches[i][j-1] & 0xFFFFFFFF;
			size_t thisMer = readMerMatches[i][j] & 0xFFFFFFFF;
			edges.push_back((prevMer << 32) + thisMer);
		}
	}
	std::sort(edges.begin(), edges.end());
	std::vector<std::vector<size_t>> outEdges;
	outEdges.resize(merCoverage.size());
	std::vector<std::vector<size_t>> inEdges;
	inEdges.resize(merCoverage.size());
	size_t edgeCount = 0;
	for (size_t i = minCoverage-1; i < edges.size(); i++)
	{
		if (i >= minCoverage && edges[i] == edges[i-minCoverage]) continue;
		if (edges[i+1-minCoverage] != edges[i]) continue;
		size_t prevMer = edges[i] >> 32;
		size_t thisMer = edges[i] & 0xFFFFFFFF;
		outEdges[prevMer].push_back(thisMer);
		inEdges[thisMer].push_back(prevMer);
		edgeCount += 1;
	}
	removeArtefactEdges(outEdges, inEdges);
	std::vector<size_t> belongsToUnitig;
	belongsToUnitig.resize(merCoverage.size(), std::numeric_limits<size_t>::max());
	std::vector<size_t> unitigEnd;
	for (size_t i = 0; i < belongsToUnitig.size(); i++)
	{
		if (belongsToUnitig[i] != std::numeric_limits<size_t>::max()) continue;
		getUnitig(belongsToUnitig, outEdges, inEdges, i, unitigEnd);
		if (!(belongsToUnitig[i] != std::numeric_limits<size_t>::max()))
		{
			std::cerr << i << std::endl;
		}
		assert(belongsToUnitig[i] != std::numeric_limits<size_t>::max());
	}
	std::vector<size_t> unitigCoverageSum;
	std::vector<size_t> unitigCoverageDiv;
	unitigCoverageDiv.resize(unitigEnd.size(), 0);
	unitigCoverageSum.resize(unitigEnd.size(), 0);
	for (size_t i = 0; i < belongsToUnitig.size(); i++)
	{
		assert(belongsToUnitig[i] != std::numeric_limits<size_t>::max());
		assert(belongsToUnitig[i] < unitigEnd.size());
		unitigCoverageSum[belongsToUnitig[i]] += merCoverage[i];
		unitigCoverageDiv[belongsToUnitig[i]] += 1;
	}
	size_t unitigEdges = 0;
	size_t unitigCount = 0;
	std::ofstream graphfile { "graph.gfa" };
	for (size_t i = 0; i < unitigEnd.size(); i++)
	{
		if (unitigCoverageDiv[i] == 0) continue;
		size_t size = (windowSize * (numWindows-0.5)) + unitigCoverageDiv[i]*windowSize/2;
		size_t coverage = unitigCoverageSum[i] / unitigCoverageDiv[i];
		if (coverage < minCoverage) continue;
		unitigCount += 1;
		graphfile << "S\t" << i << "\t*\tLN:i:" << size << "\tll:f:" << coverage << "\tFC:i:" << (coverage*size) << std::endl;
		for (auto edge : outEdges[unitigEnd[i]])
		{
			graphfile << "L\t" << i << "\t+\t" << belongsToUnitig[edge] << "\t+\t0M" << std::endl;
			unitigEdges += 1;
		}
	}
	std::cerr << merCoverage.size() << " mers" << std::endl;
	std::cerr << edgeCount << " mer-edges" << std::endl;
	std::cerr << unitigCount << " unitigs" << std::endl;
	std::cerr << unitigEdges << " unitig-edges" << std::endl;
}
