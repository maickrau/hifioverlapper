#include <limits>
#include "fastqloader.h"
#include "MatchIndex.h"
#include "MinimizerIterator.h"
#include "ReadHelper.h"

void ReadIdContainer::addNumber(uint32_t key, __uint128_t value)
{
	numbers[key].emplace_back(value);
}

void ReadIdContainer::initializeBuckets(size_t numBuckets, const std::vector<size_t>& countHashesPerBucket)
{
	assert(numBuckets == countHashesPerBucket.size());
	numbers.resize(numBuckets);
	for (size_t i = 0; i < countHashesPerBucket.size(); i++)
	{
		numbers[i].reserve(countHashesPerBucket[i]);
	}
}

const std::vector<ReadMatchposStorage>& ReadIdContainer::getMultiNumbers() const
{
	return numbers;
}

size_t ReadIdContainer::bucketsSize() const
{
	return numbers.size();
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
	windowSize(windowSize)
{
}

void MatchIndex::addHashes(uint32_t read, const std::vector<std::tuple<uint32_t, uint32_t, uint32_t>>& hashes)
{
	for (auto t : hashes)
	{
		__uint128_t val = ((__uint128_t)read << (__uint128_t)64) + ((__uint128_t)std::get<1>(t) << (__uint128_t)32) + ((__uint128_t)std::get<2>(t));
		idContainer.addNumber(std::get<0>(t), val);
	}
}

void MatchIndex::initBuckets(size_t numBuckets, const std::vector<size_t>& countHashesPerBucket)
{
	idContainer.initializeBuckets(numBuckets, countHashesPerBucket);
}
