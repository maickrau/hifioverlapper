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

size_t ReadIdContainer::multinumberSize() const
{
	return numbers.size();
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
	phmap::flat_hash_set<std::tuple<uint64_t, uint32_t, uint32_t>> hashesHere;
	iterateHashesFromRead(readSequence, [readKey, &indexMutex, fw, &hashesHere](uint64_t hash, size_t start, size_t end)
	{
		assert(end < 0x80000000);
		if (!fw)
		{
			start += 0x80000000;
			end += 0x80000000;
		}
		hashesHere.emplace(hash, start, end);
	});
	std::lock_guard<std::mutex> guard { indexMutex };
	for (auto t : hashesHere)
	{
		addMatch(std::get<0>(t), ((__uint128_t)readKey << (__uint128_t)64) + ((__uint128_t)std::get<1>(t) << (__uint128_t)32) + (__uint128_t)std::get<2>(t));
	}
}

size_t MatchIndex::multinumberSize() const
{
	return idContainer.multinumberSize();
}

size_t MatchIndex::multinumberReadCount(size_t index) const
{
	return idContainer.getMultiNumbers()[index].size();
}
