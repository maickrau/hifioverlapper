#include <limits>
#include "fastqloader.h"
#include "MatchIndex.h"
#include "MinimizerIterator.h"
#include "ReadHelper.h"

void ReadIdContainer::addNumber(uint64_t key, __uint128_t value)
{
	assert(value < FirstBit);
	auto found = firstNumberOrVectorIndex.find(key);
	if (found == firstNumberOrVectorIndex.end())
	{
		firstNumberOrVectorIndex[key] = FirstBit + value;
		return;
	}
	else
	{
		if (found->second & FirstBit)
		{
			__uint128_t oldNumber = found->second - FirstBit;
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

const std::vector<std::vector<__uint128_t>>& ReadIdContainer::getMultiNumbers() const
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

void MatchIndex::addMatch(uint64_t hash, __uint128_t readKey)
{
	idContainer.addNumber(hash, readKey);
}

void MatchIndex::addMatchesFromRead(uint32_t readKey, std::mutex& indexMutex, const std::string& readSequence)
{
	addMatchesFromReadOneWay(readKey, indexMutex, readSequence, true);
	auto rev = revCompRaw(readSequence);
	addMatchesFromReadOneWay(readKey, indexMutex, rev, false);
}

void MatchIndex::addMatchesFromReadOneWay(uint32_t readKey, std::mutex& indexMutex, const std::string& readSequence, bool fw)
{
	if (numReads <= readKey) numReads = readKey+1;
	ErrorMasking errorMasking = ErrorMasking::Microsatellite;
	std::vector<std::string> readFiles { };
	ReadpartIterator partIterator { 31, 1, errorMasking, 1, readFiles, false, "" };
	partIterator.iteratePartsOfRead("", readSequence, [this, &indexMutex, readKey, fw](const ReadInfo& read, const SequenceCharType& seq, const SequenceLengthType& poses, const std::string& raw)
	{
		if (seq.size() < numWindows * windowSize + k) return;
		phmap::flat_hash_set<std::tuple<uint64_t, uint32_t, uint32_t>> hashesHere;
		size_t rawReadLength = raw.size();
		size_t compressedReadLength = seq.size();
		iterateWindowchunks(seq, k, numWindows, windowSize, [this, &hashesHere, &poses, rawReadLength, compressedReadLength, fw](const std::vector<uint64_t>& hashes, const size_t startPos, const size_t endPos)
		{
			uint64_t totalhash = 0;
			for (auto hash : hashes)
			{
				totalhash *= 3;
				totalhash += hash;
			}
			assert(endPos > startPos);
			assert(startPos < poses.size());
			assert(endPos < poses.size());
			assert(endPos+1 < poses.size());
			uint32_t realStartPos = poses[startPos];
			uint32_t realEndPos = poses[endPos+1]-1;
			assert(realEndPos > realStartPos);
			assert((realEndPos & 0x7FFFFFFF) < rawReadLength);
			assert(realStartPos < 0x80000000);
			assert(realEndPos < 0x80000000);
			if (!fw)
			{
				realStartPos += 0x80000000;
				realEndPos += 0x80000000;
			}
			hashesHere.emplace(totalhash, realStartPos, realEndPos);
		});
		std::lock_guard<std::mutex> guard { indexMutex };
		for (auto t : hashesHere)
		{
			addMatch(std::get<0>(t), ((__uint128_t)readKey << (__uint128_t)64) + ((__uint128_t)std::get<1>(t) << (__uint128_t)32) + (__uint128_t)std::get<2>(t));
		}
	});
}
