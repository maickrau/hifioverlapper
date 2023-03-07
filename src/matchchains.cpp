#include <mutex>
#include <iostream>
#include <vector>
#include <string>
#include "ReadStorage.h"
#include "MatchIndex.h"

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
		storage.iterateReadsFromFile(file, numThreads, [&matchIndex, &indexMutex](size_t readName, const std::string& sequence)
		{
			matchIndex.addMatchesFromRead(readName, indexMutex, sequence);
		});
	}
	auto result = matchIndex.iterateMatchNames(storage.getNames(), storage.getRawReadLengths(), [](const std::string& left, const size_t leftstart, const size_t leftend, const bool leftFw, const std::string& right, const size_t rightstart, const size_t rightend, const bool rightFw)
	{
		std::cout << left << "\t" << leftstart << "\t" << leftend << "\t" << (leftFw ? "fw" : "bw") << "\t" << right << "\t" << rightstart << "\t" << rightend << "\t" << (rightFw ? "fw" : "bw") << std::endl;
	});
	std::cerr << result.numberReads << " reads" << std::endl;
	std::cerr << result.numerWindowChunks << " distinct windowchunks" << std::endl;
	std::cerr << result.totalReadChunkMatches << " read-windowchunk matches (except unique)" << std::endl;
	std::cerr << result.uniqueChunks << " windowchunks have only one read" << std::endl;
	std::cerr << result.readsWithMatch << " reads with a match" << std::endl;
	std::cerr << result.readPairMatches << " read-read matches" << std::endl;
	std::cerr << result.readChainMatches << " chain matches" << std::endl;
	std::cerr << result.totalMatches << " window matches" << std::endl;
	std::cerr << result.maxPerChunk << " max windowchunk size" << std::endl;
}
