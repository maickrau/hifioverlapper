#include "UnitigKmerCorrector.h"

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
	size_t ambiguousCoverage = 3;
	ReadpartIterator partIterator { kmerSize, windowSize, ErrorMasking::No, numThreads, readFiles, false, "" };
	UnitigKmerCorrector corrector { kmerSize };
	corrector.build(partIterator);
	std::cerr << "correcting reads" << std::endl;
	std::mutex writeMutex;
	corrector.iterateCorrectedSequences([&writeMutex](size_t index, const std::string& readname, const std::string& readseq)
	{
		std::lock_guard<std::mutex> lock { writeMutex };
		std::cout << ">" << readname << std::endl;
		std::cout << readseq << std::endl;
	});
}
