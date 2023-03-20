#include <mutex>
#include <iostream>
#include <vector>
#include <string>
#include <tuple>
#include <set>
#include <thread>
#include <concurrentqueue.h>
#include "KmerCorrector.h"
#include "ReadStorage.h"
#include "MatchIndex.h"

template <typename F>
void iterateCorrectedSequences(size_t kmerSize, size_t numThreads, const MatchIndex& matchIndex, const ReadStorage& storage, F callback)
{
	size_t windowSize = kmerSize/2;
	if (windowSize > 50) windowSize = 50;
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
	std::atomic<bool> readDone;
	readDone = false;
	moodycamel::ConcurrentQueue<std::shared_ptr<std::pair<size_t, std::vector<size_t>>>> sequenceQueue;
	std::vector<std::thread> threads;
	for (size_t i = 0; i < numThreads; i++)
	{
		threads.emplace_back([&readDone, &partIterator, &sequenceQueue, &storage, kmerSize, callback]()
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
				if (got->second.size() <= 5)
				{
					callback(got->first, storage.getRead(got->first).second, false);
					continue;
				}
				got->second.push_back(got->first);
				partIterator.setMemoryReadIterables(got->second);
				KmerCorrector corrector { kmerSize, 5, 3 };
				corrector.buildGraph(partIterator, 1);
				partIterator.setMemoryReadIterables(std::vector<size_t> { got->first });
				partIterator.iterateHashes([&corrector, callback, &got](const ReadInfo& read, const SequenceCharType& seq, const SequenceLengthType& poses, const std::string& rawSeq, const std::vector<size_t>& positions, const std::vector<HashType>& hashes)
				{
					std::pair<std::string, bool> corrected = corrector.getCorrectedSequence(rawSeq, positions, hashes);
					callback(got->first, corrected.first, corrected.second);
				});
			}
		});
	}
	std::vector<bool> processed;
	processed.resize(storage.size(), false);
	matchIndex.iterateMatchReadPairs([&sequenceQueue, &processed](size_t read, const std::unordered_set<size_t>& matches)
	{
		processed[read] = true;
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
		callback(i, storage.getRead(i).second, false);
	}
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
	iterateCorrectedSequences(correctk, numThreads, matchIndex, storage, [&matchIndex, &storage, &writeMutex, &numCorrected, &numNotCorrected](size_t readId, const std::string& readSeq, bool corrected)
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
