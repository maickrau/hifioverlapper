#include <limits>
#include "fastqloader.h"
#include "MatchIndex.h"
#include "MinimizerIterator.h"
#include "ReadHelper.h"

void ReadIdContainer::addNumber(uint64_t key, __uint128_t value)
{
	auto found = nonsingletonVectorIndex.find(key);
	if (found == nonsingletonVectorIndex.end())
	{
		return;
	}
	else
	{
		numbers[found->second].emplace_back(value);
	}
}

void ReadIdContainer::clearSeen()
{
	decltype(seenOnce) tmp;
	std::swap(tmp, seenOnce);
}

void ReadIdContainer::addSeen(uint64_t hash)
{
	if (seenOnce.count(hash) == 1)
	{
		if (nonsingletonVectorIndex.count(hash) == 1) return;
		size_t id = nonsingletonVectorIndex.size();
		nonsingletonVectorIndex[hash] = id;
	}
	seenOnce.insert(hash);
}

void ReadIdContainer::initializeBuckets()
{
	numbers.resize(nonsingletonVectorIndex.size());
}

const std::vector<ReadMatchposStorage>& ReadIdContainer::getMultiNumbers() const
{
	return numbers;
}

// __uint128_t ReadIdContainer::getFirstNumberOrBucketIndex(uint64_t key) const
// {
// 	if (firstNumberOrVectorIndex.count(key) == 0)
// 	{
// 		return 0;
// 	}
// 	else
// 	{
// 		return firstNumberOrVectorIndex.at(key);
// 	}
// }

size_t ReadIdContainer::seenSize() const
{
	return seenOnce.size();
}

size_t ReadIdContainer::bucketsSize() const
{
	return std::max(numbers.size(), nonsingletonVectorIndex.size());
}

void ReadIdContainer::clearConstructionVariablesAndCompact()
{
	{
		decltype(nonsingletonVectorIndex) tmp;
		std::swap(tmp, nonsingletonVectorIndex);
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

void MatchIndex::clearSeenAndInitializeBuckets()
{
	idContainer.clearSeen();
	idContainer.initializeBuckets();
}

void MatchIndex::addMatch(uint64_t hash, __uint128_t readKey)
{
	idContainer.addNumber(hash, readKey);
}

void MatchIndex::countSeenHashes(std::mutex& indexMutex, const std::string& readSequence)
{
	phmap::flat_hash_set<uint64_t> hashesHere;
	iterateWindowChunksFromRead(readSequence, [&hashesHere](uint32_t startPos, uint32_t endPos, uint64_t hash)
	{
		hashesHere.insert(hash);
	});
	std::lock_guard<std::mutex> guard { indexMutex };
	for (auto t : hashesHere)
	{
		idContainer.addSeen(t);
	}
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
	return idContainer.bucketsSize();
}

size_t MatchIndex::numSeenChunks() const
{
	return idContainer.seenSize();
}
