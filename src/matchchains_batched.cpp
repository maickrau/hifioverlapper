#include <mutex>
#include <iostream>
#include <vector>
#include <string>
#include <concurrentqueue.h> //https://github.com/cameron314/concurrentqueue
#include "ReadStorage.h"
#include "MatchIndex.h"

int main(int argc, char** argv)
{
	size_t numThreads = std::stoi(argv[1]);
	size_t k = std::stoull(argv[2]);
	size_t numWindows = std::stoull(argv[3]);
	size_t windowSize = std::stoull(argv[4]);
	std::string indexedReadFile = argv[5];
	std::vector<std::string> matchedReadFiles;
	for (size_t i = 6; i < argc; i++)
	{
		matchedReadFiles.emplace_back(argv[i]);
	}
	MatchIndex matchIndex { k, numWindows, windowSize };
	ReadStorage storage;
	std::mutex indexMutex;
	std::cerr << "indexing" << std::endl;
	storage.iterateReadsFromFile(indexedReadFile, numThreads, false, [&matchIndex, &indexMutex](size_t readName, const std::string& sequence)
	{
		matchIndex.addMatchesFromRead(readName, indexMutex, sequence);
	});
	size_t numWindowChunks = matchIndex.numWindowChunks();
	size_t numUniqueChunks = matchIndex.numUniqueChunks();
	// matchIndex.clearConstructionVariablesAndCompact();
	std::mutex printMutex;
	std::atomic<bool> readDone;
	readDone = false;
	std::vector<std::thread> threads;
	moodycamel::ConcurrentQueue<std::shared_ptr<FastQ>> sequenceQueue;
	size_t numMatches = 0;
	for (size_t i = 0; i < numThreads; i++)
	{
		threads.emplace_back([&readDone, &sequenceQueue, &storage, &matchIndex, &printMutex, &numMatches]()
		{
			while (true)
			{
				std::shared_ptr<FastQ> read;
				if (!sequenceQueue.try_dequeue(read))
				{
					bool tryBreaking = readDone;
					if (!sequenceQueue.try_dequeue(read))
					{
						if (tryBreaking) return;
						std::this_thread::sleep_for(std::chrono::milliseconds(10));
						continue;
					}
				}
				assert(read != nullptr);
				matchIndex.iterateMatchNamesOneRead(10000, storage.getNames(), storage.getRawReadLengths(), read->seq_id, read->sequence, [&printMutex, &numMatches](const std::string& left, const size_t leftlen, const size_t leftstart, const size_t leftend, const bool leftFw, const std::string& right, const size_t rightlen, const size_t rightstart, const size_t rightend, const bool rightFw)
				{
					std::lock_guard<std::mutex> lock { printMutex };
					std::cout << left << "\t" << leftlen << "\t" << leftstart << "\t" << leftend << "\t" << (leftFw ? "fw" : "bw") << "\t" << right << "\t" << rightlen << "\t" << rightstart << "\t" << rightend << "\t" << (rightFw ? "fw" : "bw") << std::endl;
					numMatches += 1;
				});
			}
		});
	}
	std::cerr << "matching" << std::endl;
	for (const std::string file : matchedReadFiles)
	{
		std::cerr << "Reading sequences from " << file << std::endl;
		FastQ::streamFastqFromFile(file, false, [&sequenceQueue](FastQ& read)
		{
			std::shared_ptr<FastQ> ptr = std::make_shared<FastQ>();
			std::swap(*ptr, read);
			bool queued = sequenceQueue.try_enqueue(ptr);
			if (queued) return;
			size_t triedSleeping = 0;
			while (triedSleeping < 1000)
			{
				std::this_thread::sleep_for(std::chrono::milliseconds(10));
				queued = sequenceQueue.try_enqueue(ptr);
				if (queued) return;
				triedSleeping += 1;
			}
			sequenceQueue.enqueue(ptr);
		});
	}
	readDone = true;
	for (size_t i = 0; i < threads.size(); i++)
	{
		threads[i].join();
	}
	std::cerr << numWindowChunks << " distinct windowchunks" << std::endl;
	std::cerr << numUniqueChunks << " windowchunks have only one read" << std::endl;
	std::cerr << numMatches << " chain matches" << std::endl;
}
