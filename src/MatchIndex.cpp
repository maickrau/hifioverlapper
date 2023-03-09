#include <limits>
#include "fastqloader.h"
#include "MatchIndex.h"
#include "MinimizerIterator.h"
#include "ReadHelper.h"

void ReadIdContainer::addNumber(uint64_t key, uint64_t value)
{
	assert(value < std::numeric_limits<uint64_t>::max()/2);
	auto found = firstNumberOrVectorIndex.find(key);
	if (found == firstNumberOrVectorIndex.end())
	{
		firstNumberOrVectorIndex[key] = 0x8000000000000000ull + value;
		return;
	}
	else
	{
		if (found->second & 0x8000000000000000ull)
		{
			size_t oldNumber = found->second - 0x8000000000000000ull;
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

const std::vector<std::vector<uint64_t>>& ReadIdContainer::getMultiNumbers() const
{
	return numbers;
}

size_t ReadIdContainer::size() const
{
	return firstNumberOrVectorIndex.size();
}

MatchIndex::MatchIndex(size_t k, size_t numWindows, size_t windowSize) :
	k(k),
	numWindows(numWindows),
	windowSize(windowSize),
	numReads(0)
{
}

void MatchIndex::addMatch(uint64_t hash, uint64_t readKey)
{
	idContainer.addNumber(hash, readKey);
}

void MatchIndex::addMatchesFromRead(uint32_t readKey, std::mutex& indexMutex, const std::string& readSequence)
{
	if (numReads <= readKey) numReads = readKey+1;
	ErrorMasking errorMasking = ErrorMasking::CollapseMicrosatellite;
	std::vector<std::string> readFiles { };
	ReadpartIterator partIterator { 31, 1, errorMasking, 1, readFiles, false, "" };
	partIterator.iteratePartsOfRead("", readSequence, [this, &indexMutex, readKey](const ReadInfo& read, const SequenceCharType& seq, const SequenceLengthType& poses, const std::string& raw)
	{
		if (seq.size() < numWindows * windowSize + k) return;
		phmap::flat_hash_set<std::pair<size_t, uint32_t>> hashesHere;
		size_t rawReadLength = raw.size();
		size_t compressedReadLength = seq.size();
		iterateWindowchunks(seq, k, numWindows, windowSize, [this, &hashesHere, rawReadLength, compressedReadLength](const std::vector<uint64_t>& hashes, const size_t startPos)
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
			uint32_t approxPosition = (double)startPos / (double)compressedReadLength * (double)rawReadLength;
			assert(approxPosition <= rawReadLength);
			if (!fw)
			{
				assert((startPos + numWindows * windowSize) <= compressedReadLength);
				approxPosition = (double)(compressedReadLength - (startPos + numWindows * windowSize)) / (double)compressedReadLength * (double)rawReadLength;
				assert(approxPosition <= rawReadLength);
				approxPosition += 0x80000000;
			}
			hashesHere.emplace(totalhash, approxPosition);
		});
		std::lock_guard<std::mutex> guard { indexMutex };
		for (auto hash : hashesHere)
		{
			addMatch(hash.first, ((uint64_t)readKey << 32ull) + hash.second);
		}
	});
}
