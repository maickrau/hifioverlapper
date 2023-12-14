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
			numbers[index].emplace_back(oldNumber);
			numbers[index].emplace_back(value);
		}
		else
		{
			numbers[found->second].emplace_back(value);
		}
	}
}

const std::vector<ReadMatchposStorage>& ReadIdContainer::getMultiNumbers() const
{
	return numbers;
}

__uint128_t ReadIdContainer::getFirstNumberOrBucketIndex(uint64_t key) const
{
	if (firstNumberOrVectorIndex.count(key) == 0)
	{
		return 0;
	}
	else
	{
		return firstNumberOrVectorIndex.at(key);
	}
}

size_t ReadIdContainer::size() const
{
	return firstNumberOrVectorIndex.size();
}

void ReadIdContainer::clearConstructionVariablesAndCompact()
{
	{
		decltype(firstNumberOrVectorIndex) tmp;
		std::swap(tmp, firstNumberOrVectorIndex);
	}
	for (size_t i = 0; i < numbers.size(); i++)
	{
		numbers[i].compact();
	}
	std::vector<ReadMatchposStorage> compactedNumbers;
	compactedNumbers.resize(numbers.size());
	for (size_t i = 0; i < numbers.size(); i++)
	{
		std::swap(numbers[i], compactedNumbers[i]);
	}
	std::swap(numbers, compactedNumbers);
}

Match::Match(size_t leftStartPos, size_t leftEndPos, bool leftFw, size_t rightStartPos, size_t rightEndPos, bool rightFw) :
	leftStartPos(leftStartPos),
	leftEndPos(leftEndPos),
	rightStartPos(rightStartPos),
	rightEndPos(rightEndPos),
	leftFw(leftFw),
	rightFw(rightFw)
{
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
	std::vector<std::tuple<uint64_t, uint32_t, uint32_t>> windowchunks;
	iterateWindowChunksFromRead(readSequence, [&windowchunks](uint32_t startPos, uint32_t endPos, uint64_t hash)
	{
		windowchunks.emplace_back(hash, startPos, endPos);
	});
	std::lock_guard<std::mutex> guard { indexMutex };
	if (numReads <= readKey) numReads = readKey+1;
	for (auto t : windowchunks)
	{
		addMatch(std::get<0>(t), ((__uint128_t)readKey << (__uint128_t)64) + ((__uint128_t)std::get<1>(t) << (__uint128_t)32) + (__uint128_t)std::get<2>(t));
	}
}

void MatchIndex::clearConstructionVariablesAndCompact()
{
	idContainer.clearConstructionVariablesAndCompact();
}

size_t MatchIndex::numWindowChunks() const
{
	return idContainer.size();
}

size_t MatchIndex::numUniqueChunks() const
{
	return idContainer.size() - idContainer.getMultiNumbers().size();
}
