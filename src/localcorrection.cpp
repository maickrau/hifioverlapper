#include <mutex>
#include <iostream>
#include <vector>
#include <string>
#include <tuple>
#include <set>
#include "KmerCorrector.h"
#include "ReadStorage.h"
#include "MatchIndex.h"

template <typename F>
void iterateCorrectedSequences(size_t kmerSize, const MatchIndex& matchIndex, const ReadStorage& storage, F callback)
{
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
		if (useThese.size() <= 5)
		{
			callback(read, storage.getRead(read).second, false);
			return;
		}
		partIterator.setMemoryReadIterables(useThese);
		KmerCorrector corrector { kmerSize, 5, 2 };
		corrector.buildGraph(partIterator, 1);
		partIterator.setMemoryReadIterables(std::vector<size_t> { read });
		partIterator.iterateHashes([&corrector, callback, read](const ReadInfo& read, const SequenceCharType& seq, const SequenceLengthType& poses, const std::string& rawSeq, const std::vector<size_t>& positions, const std::vector<HashType>& hashes)
		{
			std::pair<std::string, bool> corrected = corrector.getCorrectedSequence(rawSeq, positions, hashes);
			callback(read, corrected.first, corrected.second);
		});
	});
}


int main(int argc, char** argv)
{
	size_t numThreads = std::stoi(argv[1]);
	size_t searchk = std::stoull(argv[2]);
	size_t numWindows = std::stoull(argv[3]);
	size_t windowSize = std::stoull(argv[4]);
	size_t correctk = std::stoull(argv[5]);
	std::vector<std::string> readFiles;
	for (size_t i = 6; i < argc; i++)
	{
		readFiles.emplace_back(argv[i]);
	}
	MatchIndex matchIndex { searchk, numWindows, windowSize };
	ReadStorage storage;
	std::cerr << "loading reads" << std::endl;
	for (auto file : readFiles)
	{
		std::mutex indexMutex;
		storage.iterateReadsFromFile(file, numThreads, true, [&matchIndex, &indexMutex](size_t readName, const std::string& sequence)
		{
			matchIndex.addMatchesFromRead(readName, indexMutex, sequence);
		});
	}
	std::mutex writeMutex;
	size_t numCorrected = 0;
	size_t numNotCorrected = 0;
	std::cerr << "correcting" << std::endl;
	iterateCorrectedSequences(correctk, matchIndex, storage, [&matchIndex, &storage, &writeMutex, &numCorrected, &numNotCorrected](size_t readId, const std::string& readSeq, bool corrected)
	{
		std::lock_guard<std::mutex> lock { writeMutex };
		std::cout << ">" << storage.getNames()[readId] << std::endl;
		std::cout << readSeq << std::endl;
		if (corrected)
		{
			numCorrected += 1;
		}
		else
		{
			numNotCorrected += 1;
		}
	});
	std::cerr << numCorrected << " corrected" << std::endl;
	std::cerr << numNotCorrected << " not corrected" << std::endl;
}
