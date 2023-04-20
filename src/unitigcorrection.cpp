#include <mutex>
#include <atomic>
#include <iostream>
#include <concurrentqueue.h>
#include <cxxopts.hpp>
#include "MatchIndex.h"
#include "UnitigKmerCorrector.h"

template <typename F>
void iterateCorrectedReads(const UnitigKmerCorrector& corrector, const MatchIndex& matchIndex, size_t numThreads, size_t minAmbiguousCoverage, size_t minSafeCoverage, F callback)
{
	std::atomic<bool> readDone;
	readDone = false;
	moodycamel::ConcurrentQueue<std::shared_ptr<std::pair<size_t, std::vector<size_t>>>> sequenceQueue;
	std::vector<std::thread> threads;
	std::vector<bool> processed;
	processed.resize(corrector.numReads(), false);
	std::mutex processedMutex;
	for (size_t i = 0; i < numThreads; i++)
	{
		threads.emplace_back([&readDone, &sequenceQueue, &processed, &processedMutex, &corrector, minAmbiguousCoverage, minSafeCoverage, callback]()
		{
			while (true)
			{
				std::shared_ptr<std::pair<size_t, std::vector<size_t>>> got;
				if (!sequenceQueue.try_dequeue(got))
				{
					bool tryBreaking = readDone;
					if (!sequenceQueue.try_dequeue(got))
					{
						if (tryBreaking) return;
						std::this_thread::sleep_for(std::chrono::milliseconds(10));
						continue;
					}
				}
				{
					std::lock_guard lock { processedMutex };
					processed[got->first] = true;
				}
				if (got->second.size() <= 5)
				{
					callback(got->first, corrector.getName(got->first), corrector.getRawSequence(got->first), false);
					continue;
				}
				got->second.push_back(got->first);
				auto fixedContext = corrector.filterDifferentHaplotypesOut(got->first, got->second, minAmbiguousCoverage, minSafeCoverage);
				std::string corrected;
				bool changed;
				std::tie(corrected, changed) = corrector.getCorrectedSequence(got->first, fixedContext, minAmbiguousCoverage, minSafeCoverage);
				callback(got->first, corrector.getName(got->first), corrected, changed);
			}
		});
	}
	matchIndex.iterateMatchReadPairs([&sequenceQueue](size_t read, const std::unordered_set<size_t>& matches)
	{
		std::shared_ptr<std::pair<size_t, std::vector<size_t>>> ptr = std::make_shared<std::pair<size_t, std::vector<size_t>>>();
		ptr->second.insert(ptr->second.end(), matches.begin(), matches.end());
		ptr->first = read;
		bool queued = sequenceQueue.try_enqueue(ptr);
		if (queued) return;
		size_t triedSleeping = 0;
		while (triedSleeping < 100)
		{
			std::this_thread::sleep_for(std::chrono::milliseconds(10));
			queued = sequenceQueue.try_enqueue(ptr);
			if (queued) return;
			triedSleeping += 1;
		}
		sequenceQueue.enqueue(ptr);
	});
	readDone = true;
	for (size_t i = 0; i < numThreads; i++)
	{
		threads[i].join();
	}
	assert(sequenceQueue.size_approx() == 0);
	for (size_t i = 0; i < processed.size(); i++)
	{
		if (processed[i]) continue;
		callback(i, corrector.getName(i), corrector.getRawSequence(i), false);
	}
}

MatchIndex buildMatchIndex(const UnitigKmerCorrector& corrector, size_t searchk, size_t numWindows, size_t windowSize)
{
	MatchIndex matchIndex { searchk, numWindows, windowSize };
	std::mutex indexMutex;
	corrector.iterateRawSequences([&matchIndex, &indexMutex](size_t readId, const std::string& readName, const std::string& sequence)
	{
		matchIndex.addMatchesFromRead(readId, indexMutex, sequence);
	});
	return matchIndex;
}

int main(int argc, char** argv)
{
	cxxopts::Options options { "ribotin-ref" };
	options.add_options()
		("h,help", "Print help")
		("v,version", "Print version")
		("i,in", "Input reads. Multiple files can be input with -i file1.fa -i file2.fa etc (required)", cxxopts::value<std::vector<std::string>>())
		("hpc", "Collapse homopolymers in input reads")
		("searchk", "Read matching k", cxxopts::value<size_t>()->default_value("201"))
		("searchw", "Read matching window size", cxxopts::value<size_t>()->default_value("500"))
		("searchn", "Read matching window count", cxxopts::value<size_t>()->default_value("4"))
		("correctk", "Read correction k", cxxopts::value<size_t>()->default_value("101"))
		("correctw", "Read correction window size (default correctk/2)", cxxopts::value<size_t>())
		("solidcov", "Solid k-mer coverage threshold", cxxopts::value<size_t>()->default_value("5"))
		("ambiguouscov", "Ambiguous k-mer coverage threshold", cxxopts::value<size_t>()->default_value("3"))
		("t", "Number of threads (default 1)", cxxopts::value<size_t>()->default_value("1"))
	;
	auto params = options.parse(argc, argv);
	if (params.count("v") == 1)
	{
		std::cerr << "Version: " << VERSION << std::endl;
		std::exit(0);
	}
	if (params.count("h") == 1)
	{
		std::cerr << options.help() << std::endl;
		std::exit(0);
	}
	bool paramError = false;
	if (params.count("i") == 0)
	{
		std::cerr << "Input reads (-i) are required" << std::endl;
		paramError = true;
	}
	if (params["searchk"].as<size_t>() % 2 == 0)
	{
		std::cerr << "Read matching k must be odd" << std::endl;
		paramError = true;
	}
	if (params["correctk"].as<size_t>() % 2 == 0)
	{
		std::cerr << "Read matching k must be odd" << std::endl;
		paramError = true;
	}
	if (params["ambiguouscov"].as<size_t>() > params["solidcov"].as<size_t>())
	{
		std::cerr << "Ambiguous coverage threshold can't be higher than solid coverage threshold" << std::endl;
		paramError = true;
	}
	if (paramError)
	{
		std::abort();
	}
	size_t numThreads = params["t"].as<size_t>();
	size_t searchk = params["searchk"].as<size_t>();
	size_t numWindows = params["searchn"].as<size_t>();
	size_t searchwindowSize = params["searchw"].as<size_t>();
	size_t correctk = params["correctk"].as<size_t>();
	size_t correctw;
	if (params.count("correctw") == 1)
	{
		correctw = params["correctw"].as<size_t>();
	}
	else
	{
		correctw = correctk/2;
	}
	std::vector<std::string> readFiles = params["i"].as<std::vector<std::string>>();
	size_t minSafeCoverage = params["solidcov"].as<size_t>();
	size_t minAmbiguousCoverage = params["ambiguouscov"].as<size_t>();
	ErrorMasking masking = ErrorMasking::No;
	if (params.count("hpc") == 1) masking = ErrorMasking::Collapse;
	ReadpartIterator partIterator { correctk, correctw, masking, numThreads, readFiles, false, "" };
	UnitigKmerCorrector corrector { correctk };
	corrector.build(partIterator);
	std::cerr << "build match index" << std::endl;
	MatchIndex matchIndex = buildMatchIndex(corrector, searchk, numWindows, searchwindowSize);
	std::cerr << "correcting reads" << std::endl;
	std::mutex writeMutex;
	size_t numCorrected = 0;
	size_t numNotCorrected = 0;
	iterateCorrectedReads(corrector, matchIndex, numThreads, minAmbiguousCoverage, minSafeCoverage, [&writeMutex, &numCorrected, &numNotCorrected](size_t index, const std::string& readname, const std::string& readseq, bool corrected)
	{
		std::lock_guard<std::mutex> lock { writeMutex };
		std::cout << ">" << readname << std::endl;
		std::cout << readseq << std::endl;
		if (corrected)
		{
			numCorrected += 1;
		}
		else
		{
			numNotCorrected += 1;
		}
	});
	std::cerr << numCorrected << " reads corrected" << std::endl;
	std::cerr << numNotCorrected << " reads not corrected" << std::endl;
}
