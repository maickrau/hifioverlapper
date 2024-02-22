#include "MinimizerIterator.h"

MinimizerIterator::MinimizerIterator(size_t k, bool fw) :
	k(k),
	lastHash(k),
	fw(fw)
{
}

void MinimizerIterator::init(const MBG::SequenceCharType& start, size_t posOffset)
{
	assert(start.size() >= k);
	for (size_t i = 0; i < k; i++)
	{
		lastHash.addChar(start[i]);
	}
	windowSize = start.size() - k + 1;
	windowKmers.emplace_back(posOffset, fw ? lastHash.getFwHash() : lastHash.getBwHash());
	for (size_t i = k; i < start.size(); i++)
	{
		lastHash.addChar(start[i]);
		lastHash.removeChar(start[i-k]);
		uint64_t hash = fw ? lastHash.getFwHash() : lastHash.getBwHash();
		while (windowKmers.size() > 0 && windowKmers.back().second > hash) windowKmers.pop_back();
		windowKmers.emplace_back(i-k+1+posOffset, hash);
	}
	pos = start.size()-k+1+posOffset;
}

void MinimizerIterator::moveChar(uint16_t added, uint16_t removed)
{
	lastHash.addChar(added);
	lastHash.removeChar(removed);
	// comparison rearranged, really first <= pos - windowSize but move windowSize because of underflow
	while (windowKmers.size() > 0 && windowKmers.front().first + windowSize <= pos) windowKmers.erase(windowKmers.begin());
	while (windowKmers.size() > 0 && windowKmers.back().second > (fw ? lastHash.getFwHash() : lastHash.getBwHash())) windowKmers.pop_back();
	windowKmers.emplace_back(pos, fw ? lastHash.getFwHash() : lastHash.getBwHash());
	pos += 1;
}

size_t MinimizerIterator::minimizerPosition() const
{
	assert(windowKmers.size() >= 1);
	return windowKmers[0].first;
}

uint64_t MinimizerIterator::minimizerHash() const
{
	assert(windowKmers.size() >= 1);
	return windowKmers[0].second;
}
