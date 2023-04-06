#include <mutex>
#include <atomic>
#include <iostream>
#include <concurrentqueue.h>
#include "MatchIndex.h"
#include "UnitigKmerCorrector.h"
#include "/usr/local/include/valgrind/callgrind.h"

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
	CALLGRIND_START_INSTRUMENTATION;
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
					callback(got->first, corrector.getName(got->first), corrector.getRawSequence(got->first));
					continue;
				}
				got->second.push_back(got->first);
				std::string corrected = corrector.getCorrectedSequence(got->first, got->second, minAmbiguousCoverage, minSafeCoverage);
				callback(got->first, corrector.getName(got->first), corrected);
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
	CALLGRIND_STOP_INSTRUMENTATION;
	assert(sequenceQueue.size_approx() == 0);
	for (size_t i = 0; i < processed.size(); i++)
	{
		if (processed[i]) continue;
		callback(i, corrector.getName(i), corrector.getRawSequence(i));
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
	size_t numThreads = std::stoi(argv[1]);
	size_t searchk = std::stoull(argv[2]);
	size_t numWindows = std::stoull(argv[3]);
	size_t searchwindowSize = std::stoull(argv[4]);
	size_t correctk = std::stoull(argv[5]);
	std::vector<std::string> readFiles;
	for (size_t i = 6; i < argc; i++)
	{
		readFiles.emplace_back(argv[i]);
	}
	size_t correctwindowSize = correctk/2;
	size_t minSafeCoverage = 5;
	size_t minAmbiguousCoverage = 3;
	ReadpartIterator partIterator { correctk, correctwindowSize, ErrorMasking::No, numThreads, readFiles, false, "" };
	UnitigKmerCorrector corrector { correctk };
	corrector.build(partIterator);
	std::cerr << "build match index" << std::endl;
	MatchIndex matchIndex = buildMatchIndex(corrector, searchk, numWindows, searchwindowSize);
	std::cerr << "correcting reads" << std::endl;
	std::mutex writeMutex;
	iterateCorrectedReads(corrector, matchIndex, numThreads, minAmbiguousCoverage, minSafeCoverage, [&writeMutex](size_t index, const std::string& readname, const std::string& readseq)
	{
		std::lock_guard<std::mutex> lock { writeMutex };
		std::cout << ">" << readname << std::endl;
		std::cout << readseq << std::endl;
	});
}
