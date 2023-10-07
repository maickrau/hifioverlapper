#ifndef MinimizerIterator_h
#define MinimizerIterator_h

#include <vector>
#include <tuple>
#include <cstdint>
#include "MBGCommon.h"
#include "FastHasher.h"

class MinimizerIterator
{
public:
	MinimizerIterator() = default;
	MinimizerIterator(size_t k);
	void init(const SequenceCharType& start, size_t posOffset);
	void moveChar(uint16_t added, uint16_t removed);
	size_t minimizerPosition() const;
	uint64_t minimizerHash() const;
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
		windowIterators[i].init(startWindow, i * windowSize);
	}
	std::vector<uint64_t> hashesHere;
	hashesHere.resize(numWindows);
	for (size_t i = 0; i < numWindows; i++)
	{
		hashesHere[i] = windowIterators[i].minimizerHash();
	}
	callback(hashesHere, windowIterators[0].minimizerPosition(), windowIterators.back().minimizerPosition()+k-1);
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
			callback(hashesHere, windowIterators[0].minimizerPosition(), windowIterators.back().minimizerPosition()+k-1);
		}
	}
}

#endif