#include <limits>
#include "fastqloader.h"
#include "MatchIndex.h"
#include "MinimizerIterator.h"
#include "ReadHelper.h"

void ReadIdContainer::addNumber(uint64_t key, uint32_t value)
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

const std::vector<std::vector<uint32_t>>& ReadIdContainer::getMultiNumbers() const
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
	windowSize(windowSize)
{
}

void MatchIndex::addMatchesFromFile(size_t numThreads, std::string file)
{
	ErrorMasking errorMasking = ErrorMasking::CollapseMicrosatellite;
	std::vector<std::string> readFiles { file };
	ReadpartIterator partIterator { 31, 1, errorMasking, numThreads, readFiles, false, "" };
	std::mutex indexMutex;
	partIterator.iterateParts([this, &indexMutex](const ReadInfo& read, const SequenceCharType& seq, const SequenceLengthType& poses, const std::string& raw)
	{
		if (seq.size() < numWindows * windowSize + k) return;
		size_t readName;
		{
			std::lock_guard<std::mutex> guard { indexMutex };
			readName = names.size();
			names.push_back(read.readName.first);
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
			idContainer.addNumber(hash, readName);
		}
	});

}
