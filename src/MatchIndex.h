#ifndef MatchIndex_h
#define MatchIndex_h

#include <string>
#include <cstdint>
#include <vector>
#include <mutex>
#include <tuple>
#include <phmap.h>

class ReadIdContainer
{
public:
	ReadIdContainer() = default;
	void addNumber(uint64_t key, uint64_t value);
	const std::vector<std::vector<uint64_t>>& getMultiNumbers() const;
	size_t size() const;
private:
	phmap::flat_hash_map<uint64_t, uint64_t> firstNumberOrVectorIndex;
	std::vector<std::vector<uint64_t>> numbers;
};

class IterationInfo
{
public:
	IterationInfo() = default;
	size_t numberReads;
	size_t numerWindowChunks;
	size_t totalReadChunkMatches;
	size_t uniqueChunks;
	size_t readsWithMatch;
	size_t maxPerChunk;
	size_t totalMatches;
	size_t readPairMatches;
	size_t readChainMatches;
};

class MatchIndex
{
public:
	MatchIndex(size_t k, size_t numWindows, size_t windowSize);
	void addMatch(uint64_t hash, uint64_t readKey);
	void addMatchesFromRead(uint32_t readKey, std::mutex& indexMutex, const std::string& readSequence);
	template <typename F>
	IterationInfo iterateMatches(F callback) const
	{
		IterationInfo result;
		result.numberReads = numReads;
		result.numerWindowChunks = idContainer.size();
		std::vector<std::vector<std::pair<uint32_t, uint32_t>>> hashesPerRead;
		hashesPerRead.resize(numReads);
		size_t maxPerChunk = 0;
		const auto& numbers = idContainer.getMultiNumbers();
		size_t uniqueChunks = idContainer.size() - numbers.size();
		size_t totalReadChunkMatches = 0;
		for (size_t i = 0; i < numbers.size(); i++)
		{
			totalReadChunkMatches += numbers[i].size();
			maxPerChunk = std::max(maxPerChunk, numbers[i].size());
			for (uint64_t readkey : numbers[i])
			{
				uint32_t readname = readkey >> 32;
				uint32_t readPos = readkey;
				hashesPerRead[readname].emplace_back(i, readPos);
			}
		}
		result.totalReadChunkMatches = totalReadChunkMatches;
		result.uniqueChunks = uniqueChunks;
		size_t readsWithMatch = 0;
		size_t totalMatches = 0;
		size_t readPairMatches = 0;
		for (size_t i = 0; i < hashesPerRead.size(); i++)
		{
			std::vector<std::tuple<uint32_t, uint32_t, uint32_t>> matches;
			for (std::pair<uint32_t, uint32_t> pair : hashesPerRead[i])
			{
				uint32_t hash = pair.first;
				uint32_t pos = pair.second;
				for (uint64_t readkey : numbers[hash])
				{
					uint32_t read = readkey >> 32;
					uint32_t otherpos = readkey;
					if (read == i) continue;
					if (read < i) continue;
					matches.emplace_back(read, pos, otherpos);
				}
			}
			if (matches.size() > 0) readsWithMatch += 1;
			std::sort(matches.begin(), matches.end());
			totalMatches += matches.size();
			phmap::flat_hash_set<uint32_t> matchingReads;
			for (auto match : matches)
			{
				size_t otherread = std::get<0>(match);
				uint32_t thispos = std::get<1>(match);
				uint32_t otherpos = std::get<2>(match);
				bool thisbw = thispos & 0x80000000;
				bool otherbw = otherpos & 0x80000000;
				thispos &= 0x7FFFFFFF;
				otherpos &= 0x7FFFFFFF;
				callback(i, thispos, !thisbw, otherread, otherpos, !otherbw);
				matchingReads.insert(otherread);
			}
			readPairMatches += matchingReads.size();
		}
		result.readPairMatches = readPairMatches;
		result.readsWithMatch = readsWithMatch;
		result.maxPerChunk = maxPerChunk;
		result.totalMatches = totalMatches;
		return result;
	}
	template <typename F>
	IterationInfo iterateMatchChains(const std::vector<size_t>& rawReadLengths, F callback) const
	{
		size_t currentLeftRead = std::numeric_limits<size_t>::max();
		size_t currentRightRead = std::numeric_limits<size_t>::max();
		std::vector<std::pair<size_t, size_t>> currentFwMatches;
		std::vector<std::pair<size_t, size_t>> currentBwMatches;
		size_t totalChains = 0;
		auto result = iterateMatches([this, &totalChains, &currentLeftRead, &currentRightRead, &currentFwMatches, &currentBwMatches, &rawReadLengths, callback](size_t leftread, size_t leftpos, bool leftFw, size_t rightread, size_t rightpos, bool rightFw)
		{
			if (leftread != currentLeftRead || rightread != currentRightRead)
			{
				totalChains += iterateChains(currentLeftRead, currentRightRead, currentFwMatches, currentBwMatches, callback);
				currentFwMatches.clear();
				currentBwMatches.clear();
				currentLeftRead = leftread;
				currentRightRead = rightread;
			}
			if (!leftFw)
			{
				leftpos = rawReadLengths[leftread] - leftpos;
			}
			if (rightFw)
			{
				currentFwMatches.emplace_back(leftpos, rightpos);
			}
			else
			{
				currentBwMatches.emplace_back(leftpos, rightpos);
			}
		});
		totalChains += iterateChains(currentLeftRead, currentRightRead, currentFwMatches, currentBwMatches, callback);
		result.readChainMatches = totalChains;
		return result;
	}
	template <typename F>
	size_t iterateChains(size_t leftread, size_t rightread, std::vector<std::pair<size_t, size_t>>& fwMatches, std::vector<std::pair<size_t, size_t>>& bwMatches, F callback) const
	{
		const size_t maxChainDiagonalDifference = 500;
		size_t result = 0;
		std::sort(bwMatches.begin(), bwMatches.end(), [](std::pair<size_t, size_t> left, std::pair<size_t, size_t> right) { return (int64_t)left.second - (int64_t)left.first < (int64_t)right.second - (int64_t)right.first; });
		if (fwMatches.size() > 0)
		{
			std::sort(fwMatches.begin(), fwMatches.end(), [](std::pair<size_t, size_t> left, std::pair<size_t, size_t> right) { return (int64_t)left.second - (int64_t)left.first < (int64_t)right.second - (int64_t)right.first; });
			size_t leftStart = fwMatches[0].first;
			size_t leftEnd = fwMatches[0].first;
			size_t rightStart = fwMatches[0].second;
			size_t rightEnd = fwMatches[0].second;
			for (size_t i = 1; i < fwMatches.size(); i++)
			{
				if ((int64_t)fwMatches[i].second - (int64_t)fwMatches[i].first > (int64_t)fwMatches[i+1].second - (int64_t)fwMatches[i+1].first + (int64_t)maxChainDiagonalDifference)
				{
					callback(leftread, leftStart, leftEnd, true, rightread, rightStart, rightEnd, true);
					result += 1;
					leftStart = fwMatches[i].first;
					leftEnd = fwMatches[i].first;
					rightStart = fwMatches[i].second;
					rightEnd = fwMatches[i].second;
				}
				else
				{
					leftStart = std::min(leftStart, fwMatches[i].first);
					leftEnd = std::max(leftEnd, fwMatches[i].first);
					rightStart = std::min(rightStart, fwMatches[i].second);
					rightEnd = std::max(rightEnd, fwMatches[i].second);
				}
			}
			callback(leftread, leftStart, leftEnd, true, rightread, rightStart, rightEnd, true);
			result += 1;
		}
		if (bwMatches.size() > 0)
		{
			std::sort(bwMatches.begin(), bwMatches.end(), [](std::pair<size_t, size_t> left, std::pair<size_t, size_t> right) { return (int64_t)left.second - (int64_t)left.first < (int64_t)right.second - (int64_t)right.first; });
			size_t leftStart = bwMatches[0].first;
			size_t leftEnd = bwMatches[0].first;
			size_t rightStart = bwMatches[0].second;
			size_t rightEnd = bwMatches[0].second;
			for (size_t i = 1; i < bwMatches.size(); i++)
			{
				if ((int64_t)bwMatches[i].second - (int64_t)bwMatches[i].first > (int64_t)bwMatches[i+1].second - (int64_t)bwMatches[i+1].first + (int64_t)maxChainDiagonalDifference)
				{
					callback(leftread, leftStart, leftEnd, true, rightread, rightStart, rightEnd, false);
					result += 1;
					leftStart = bwMatches[i].first;
					leftEnd = bwMatches[i].first;
					rightStart = bwMatches[i].second;
					rightEnd = bwMatches[i].second;
				}
				else
				{
					leftStart = std::min(leftStart, bwMatches[i].first);
					leftEnd = std::max(leftEnd, bwMatches[i].first);
					rightStart = std::min(rightStart, bwMatches[i].second);
					rightEnd = std::max(rightEnd, bwMatches[i].second);
				}
			}
			callback(leftread, leftStart, leftEnd, true, rightread, rightStart, rightEnd, false);
			result += 1;
		}
		return result;
	}
	template <typename F>
	IterationInfo iterateMatchNames(const std::vector<std::string>& names, const std::vector<size_t>& rawReadLengths, F callback) const
	{
		return iterateMatchChains(rawReadLengths, [&names, callback](size_t left, size_t leftStart, size_t leftEnd, bool leftFw, size_t right, size_t rightStart, size_t rightEnd, bool rightFw)
		{
			callback(names[left], leftStart, leftEnd, leftFw, names[right], rightStart, rightEnd, rightFw);
		});
	}
private:
	size_t k;
	size_t numWindows;
	size_t windowSize;
	size_t numReads;
	ReadIdContainer idContainer;
};

#endif
