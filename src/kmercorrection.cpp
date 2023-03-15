#include <mutex>
#include <iostream>
#include <vector>
#include <string>
#include <tuple>
#include <set>
#include "ReadHelper.h"
#include "KmerCorrector.h"

template <typename F>
void iterateCorrectedSequences(const ReadpartIterator& partIterator, const KmerCorrector& corrector, F callback)
{
	partIterator.iterateHashes([&corrector, callback](const ReadInfo& read, const SequenceCharType& seq, const SequenceLengthType& poses, const std::string& rawSeq, const std::vector<size_t>& positions, const std::vector<HashType>& hashes)
	{
		std::pair<std::string, bool> corrected = corrector.getCorrectedSequence(rawSeq, positions, hashes);
		callback(read.readName.first, corrected.first, corrected.second);
	});
}

int main(int argc, char** argv)
{
	size_t numThreads = std::stoi(argv[1]);
	size_t kmerSize = std::stoull(argv[2]);
	std::vector<std::string> readFiles;
	for (size_t i = 3; i < argc; i++)
	{
		readFiles.emplace_back(argv[i]);
	}
	size_t windowSize = kmerSize/2;
	size_t solidCoverage = 5;
	size_t ambiguousCoverage = 2;
	ReadpartIterator partIterator { kmerSize, windowSize, ErrorMasking::No, numThreads, readFiles, false, "" };
	KmerCorrector corrector { kmerSize, solidCoverage, ambiguousCoverage };
	corrector.buildGraph(partIterator, numThreads);
	std::cerr << "correcting reads" << std::endl;
	std::mutex writeMutex;
	size_t countCorrected = 0;
	size_t countNotCorrected = 0;
	iterateCorrectedSequences(partIterator, corrector, [&writeMutex, &countCorrected, &countNotCorrected](const std::string& readname, const std::string& readseq, const bool corrected)
	{
		std::lock_guard<std::mutex> lock { writeMutex };
		if (corrected)
		{
			countCorrected += 1;
		}
		else
		{
			countNotCorrected += 1;
		}
		std::cout << ">" << readname << std::endl;
		std::cout << readseq << std::endl;
	});
	std::cerr << countCorrected << " reads corrected" << std::endl;
	std::cerr << countNotCorrected << " reads not corrected" << std::endl;
}
