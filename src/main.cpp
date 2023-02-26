#include <mutex>
#include <iostream>
#include <vector>
#include <string>
#include <phmap.h>
#include "FastHasher.h"
#include "ReadHelper.h"

class MinimizerIterator
{
public:
	MinimizerIterator() = default;
	MinimizerIterator(size_t k) :
		k(k),
		lastHash(k)
	{
	}
	void init(const SequenceCharType& start)
	{
		assert(start.size() >= k);
		for (size_t i = 0; i < k; i++)
		{
			lastHash.addChar(start[i]);
		}
		windowSize = start.size() - k + 1;
		windowKmers.emplace_back(0, lastHash.hash());
		for (size_t i = k; i < start.size(); i++)
		{
			lastHash.addChar(start[i]);
			lastHash.removeChar(start[i-k]);
			uint64_t hash = lastHash.hash();
			while (windowKmers.size() > 0 && windowKmers.back().second > hash) windowKmers.pop_back();
			windowKmers.emplace_back(i-k+1, lastHash.hash());
		}
		pos = start.size()-k+1;
	}
	void moveChar(uint16_t added, uint16_t removed)
	{
		lastHash.addChar(added);
		lastHash.removeChar(removed);
		while (windowKmers.size() > 0 && windowKmers.front().first <= pos - windowSize) windowKmers.erase(windowKmers.begin());
		while (windowKmers.size() > 0 && windowKmers.back().second > lastHash.hash()) windowKmers.pop_back();
		windowKmers.emplace_back(pos, lastHash.hash());
		pos += 1;
	}
	size_t minimizerPosition() const
	{
		assert(windowKmers.size() >= 1);
		return windowKmers[0].first;
	}
	uint64_t minimizerHash() const
	{
		assert(windowKmers.size() >= 1);
		return windowKmers[0].second;
	}
private:
	size_t k;
	size_t windowSize;
	std::vector<std::pair<size_t, uint64_t>> windowKmers;
	FastHasher lastHash;
	size_t pos;
};

template <typename F>
void iterateWindowchunks(const SequenceCharType& seq, size_t k, size_t numWindows, size_t windowSize, F callback)
{
	if (seq.size() < numWindows * windowSize + k) return;
	std::vector<MinimizerIterator> windowIterators;
	for (size_t i = 0; i < numWindows; i++)
	{
		windowIterators.emplace_back(k);
		SequenceCharType startWindow { seq.begin() + i * windowSize, seq.begin() + (i+1) * windowSize + k };
		windowIterators[i].init(startWindow);
	}
	std::vector<uint64_t> hashesHere;
	hashesHere.resize(numWindows);
	for (size_t i = 0; i < numWindows; i++)
	{
		hashesHere[i] = windowIterators[i].minimizerHash();
	}
	callback(hashesHere);
	for (size_t i = 0; i + numWindows * windowSize + k < seq.size(); i++)
	{
		bool changed = false;
		for (size_t j = 0; j < numWindows; j++)
		{
			windowIterators[j].moveChar(seq[(j+1)*windowSize + k + i], seq[(j+1)*windowSize+i]);
			if (windowIterators[j].minimizerHash() != hashesHere[j])
			{
				changed = true;
				hashesHere[j] = windowIterators[j].minimizerHash();
			}
		}
		if (changed)
		{
			callback(hashesHere);
		}
	}
}

class ReadIdContainer
{
public:
	void addNumber(uint64_t key, uint32_t value)
	{
		assert(value < std::numeric_limits<uint32_t>::max()/2);
		auto found = firstNumberOrVectorIndex.find(key);
		if (found == firstNumberOrVectorIndex.end())
		{
			firstNumberOrVectorIndex[key] = 0x80000000 + value;
			return;
		}
		else
		{
			if (found->second & 0x80000000)
			{
				size_t oldNumber = found->second - 0x80000000;
				size_t index = numbers.size();
				numbers.emplace_back();
				firstNumberOrVectorIndex[key] = index;
				numbers[index].push_back(oldNumber);
				numbers[index].push_back(value);
			}
			else
			{
				numbers[found->second].push_back(value);
			}
		}
	}
	const std::vector<std::vector<uint32_t>>& getMultiNumbers() const
	{
		return numbers;
	}
	size_t size() const
	{
		return firstNumberOrVectorIndex.size();
	}
private:
	phmap::flat_hash_map<uint64_t, uint32_t> firstNumberOrVectorIndex;
	std::vector<std::vector<uint32_t>> numbers;
};

class IterationInfo
{
public:
	IterationInfo() = default;
	size_t numberReads;
	size_t numerWindowChunks;
	size_t totalReadChunkMatches;
	size_t uniqueChunks;
	size_t readsWithMatch;
	size_t maxPerChunk;
	size_t totalMatches;
};

template <typename F>
IterationInfo iterateMatches(size_t numThreads, size_t k, size_t numWindows, size_t windowSize, std::vector<std::string> readFiles, F callback)
{
	IterationInfo result;
	std::vector<std::string> readNames;
	ReadIdContainer readsPerWindowchunk;
	ErrorMasking errorMasking = ErrorMasking::CollapseMicrosatellite;
	ReadpartIterator partIterator { 31, 1, errorMasking, numThreads, readFiles, false, "" };
	std::mutex indexMutex;
	partIterator.iterateParts([&readNames, &readsPerWindowchunk, &indexMutex, k, numWindows, windowSize](const ReadInfo& read, const SequenceCharType& seq, const SequenceLengthType& poses, const std::string& raw)
	{
		if (seq.size() < numWindows * windowSize + k) return;
		size_t readName;
		{
			std::lock_guard<std::mutex> guard { indexMutex };
			readName = readNames.size();
			readNames.push_back(read.readName.first);
		}
		phmap::flat_hash_set<size_t> hashesHere;
		iterateWindowchunks(seq, k, numWindows, windowSize, [&hashesHere, readName](const std::vector<uint64_t>& hashes)
		{
			bool fw = true;
			for (size_t i = 0; i < hashes.size()/2; i++)
			{
				if (hashes[i] > hashes[hashes.size()-1-i]) break;
				if (hashes[i] < hashes[hashes.size()-1-i])
				{
					fw = false;
					break;
				}
			}
			size_t totalhash = 0;
			if (fw)
			{
				for (auto hash : hashes)
				{
					totalhash += hash;
					totalhash *= 3;
				}
			}
			else
			{
				for (size_t i = hashes.size()-1; i < hashes.size(); i--)
				{
					totalhash += hashes[i];
					totalhash *= 3;
				}
			}
			hashesHere.emplace(totalhash);
		});
		std::lock_guard<std::mutex> guard { indexMutex };
		for (auto hash : hashesHere)
		{
			readsPerWindowchunk.addNumber(hash, readName);
		}
	});
	result.numberReads = readNames.size();
	result.numerWindowChunks = readsPerWindowchunk.size();
	std::vector<std::vector<uint32_t>> hashesPerRead;
	hashesPerRead.resize(readNames.size());
	size_t maxPerChunk = 0;
	const auto& numbers = readsPerWindowchunk.getMultiNumbers();
	size_t uniqueChunks = readsPerWindowchunk.size() - numbers.size();
	size_t totalReadChunkMatches = 0;
	for (size_t i = 0; i < numbers.size(); i++)
	{
		totalReadChunkMatches += numbers[i].size();
		maxPerChunk = std::max(maxPerChunk, numbers[i].size());
		for (auto readname : numbers[i])
		{
			hashesPerRead[readname].push_back(i);
		}
	}
	result.totalReadChunkMatches = totalReadChunkMatches;
	result.uniqueChunks = uniqueChunks;
	size_t readsWithMatch = 0;
	size_t totalMatches = 0;
	for (size_t i = 0; i < hashesPerRead.size(); i++)
	{
		phmap::flat_hash_set<size_t> matches;
		for (auto hash : hashesPerRead[i])
		{
			for (auto read : numbers[hash])
			{
				if (read == i) continue;
				if (read < i) continue;
				matches.emplace(read);
			}
		}
		if (matches.size() > 0) readsWithMatch += 1;
		totalMatches += matches.size();
		for (auto match : matches)
		{
			callback(readNames[i], readNames[match]);
		}
	}
	result.readsWithMatch = readsWithMatch;
	result.maxPerChunk = maxPerChunk;
	result.totalMatches = totalMatches;
	return result;
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
	auto result = iterateMatches(numThreads, k, numWindows, windowSize, readFiles, [](const std::string& left, const std::string& right)
	{
		std::cout << left << "\t" << right << std::endl;
	});
	std::cerr << result.numberReads << " reads" << std::endl;
	std::cerr << result.numerWindowChunks << " distinct windowchunks" << std::endl;
	std::cerr << result.totalReadChunkMatches << " read-chunk matches (except unique)" << std::endl;
	std::cerr << result.uniqueChunks << " windowchunks have only one read" << std::endl;
	std::cerr << result.readsWithMatch << " reads with a match" << std::endl;
	std::cerr << "max per chunk: " << result.maxPerChunk << std::endl;
	std::cerr << "num matches: " << result.totalMatches << std::endl;
}
