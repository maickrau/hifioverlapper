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
	MinimizerIterator(size_t k, bool fw);
	void init(const MBG::SequenceCharType& start, size_t posOffset);
	void moveChar(uint16_t added, uint16_t removed);
	size_t minimizerPosition() const;
	uint64_t minimizerHash() const;
private:
	size_t k;
	size_t windowSize;
	std::vector<std::pair<size_t, uint64_t>> windowKmers;
	MBG::FastHasher lastHash;
	size_t pos;
	bool fw;
};

template <typename F>
void iterateWindowchunks(const MBG::SequenceCharType& seq, size_t k, size_t numWindows, size_t windowSize, F callback)
{
	if (seq.size() < numWindows * windowSize + k) return;
	std::vector<MinimizerIterator> windowIterators;
	for (size_t i = 0; i < numWindows; i++)
	{
		windowIterators.emplace_back(k, i < numWindows / 2);
		MBG::SequenceCharType startWindow { seq.begin() + i * windowSize, seq.begin() + (i+1) * windowSize + k };
		windowIterators[i].init(startWindow, i * windowSize);
	}
//	std::vector<uint64_t> hashesHere;
	std::vector<uint64_t> positionsHere;
//	hashesHere.resize(numWindows);
	positionsHere.resize(numWindows);
	for (size_t i = 0; i < numWindows; i++)
	{
//		hashesHere[i] = windowIterators[i].minimizerHash();
		positionsHere[i] = windowIterators[i].minimizerPosition();
	}
//	callback(hashesHere, windowIterators[0].minimizerPosition(), windowIterators.back().minimizerPosition()+k-1);
	callback(positionsHere);
	for (size_t i = 0; i + numWindows * windowSize + k < seq.size(); i++)
	{
		bool changed = false;
		for (size_t j = 0; j < numWindows; j++)
		{
			windowIterators[j].moveChar(seq[(j+1)*windowSize + k + i], seq[(j+1)*windowSize+i]);
			if (windowIterators[j].minimizerPosition() != positionsHere[j])
			{
				changed = true;
				positionsHere[j] = windowIterators[j].minimizerPosition();
			}
		}
		if (changed)
		{
//			callback(hashesHere, windowIterators[0].minimizerPosition(), windowIterators.back().minimizerPosition()+k-1);
			callback(positionsHere);
		}
	}
}

#endif