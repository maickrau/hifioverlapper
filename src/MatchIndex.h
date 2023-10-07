#ifndef MatchIndex_h
#define MatchIndex_h

#include <string>
#include <cstdint>
#include <vector>
#include <mutex>
#include <tuple>
#include <unordered_set>
#include <phmap.h>

class ReadIdContainer
{
public:
	ReadIdContainer() = default;
	void addNumber(uint64_t key, __uint128_t value);
	const std::vector<std::vector<__uint128_t>>& getMultiNumbers() const;
	size_t size() const;
private:
	static constexpr __uint128_t FirstBit = (__uint128_t)1 << (__uint128_t)127;
	phmap::flat_hash_map<uint64_t, __uint128_t> firstNumberOrVectorIndex;
	std::vector<std::vector<__uint128_t>> numbers;
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
	void addMatch(uint64_t hash, __uint128_t readKey);
	void addMatchesFromRead(uint32_t readKey, std::mutex& indexMutex, const std::string& readSequence);
	template <typename F>
	IterationInfo iterateMatches(bool alsoSmaller, F callback) const
	{
		IterationInfo result;
		result.numberReads = numReads;
		result.numerWindowChunks = idContainer.size();
		std::vector<std::vector<std::tuple<uint32_t, uint32_t, uint32_t>>> hashesPerRead;
		hashesPerRead.resize(numReads);
		size_t maxPerChunk = 0;
		const auto& numbers = idContainer.getMultiNumbers();
		size_t uniqueChunks = idContainer.size() - numbers.size();
		size_t totalReadChunkMatches = 0;
		for (size_t i = 0; i < numbers.size(); i++)
		{
			totalReadChunkMatches += numbers[i].size();
			maxPerChunk = std::max(maxPerChunk, numbers[i].size());
			for (__uint128_t readkey : numbers[i])
			{
				uint32_t readname = readkey >> (__uint128_t)64;
				uint32_t readStartPos = readkey >> (__uint128_t)32;
				uint32_t readEndPos = readkey;
				hashesPerRead[readname].emplace_back(i, readStartPos, readEndPos);
			}
		}
		result.totalReadChunkMatches = totalReadChunkMatches;
		result.uniqueChunks = uniqueChunks;
		size_t readsWithMatch = 0;
		size_t totalMatches = 0;
		size_t readPairMatches = 0;
		for (size_t i = 0; i < hashesPerRead.size(); i++)
		{
			std::vector<std::tuple<uint32_t, uint32_t, uint32_t, uint32_t, uint32_t>> matches;
			for (std::tuple<uint32_t, uint32_t, uint32_t> t : hashesPerRead[i])
			{
				uint32_t hash = std::get<0>(t);
				uint32_t startpos = std::get<1>(t);
				uint32_t endpos = std::get<2>(t);
				for (__uint128_t readkey : numbers[hash])
				{
					uint32_t read = readkey >> (__uint128_t)64;
					uint32_t otherStartPos = readkey >> (__uint128_t)32;
					uint32_t otherEndPos = readkey;
					if (read == i) continue;
					if (!alsoSmaller && read < i) continue;
					matches.emplace_back(read, startpos, endpos, otherStartPos, otherEndPos);
				}
			}
			if (matches.size() > 0) readsWithMatch += 1;
			std::sort(matches.begin(), matches.end());
			totalMatches += matches.size();
			phmap::flat_hash_set<uint32_t> matchingReads;
			for (auto match : matches)
			{
				size_t otherread = std::get<0>(match);
				uint32_t thisStartPos = std::get<1>(match);
				uint32_t thisEndPos = std::get<2>(match);
				uint32_t otherStartPos = std::get<3>(match);
				uint32_t otherEndPos = std::get<4>(match);
				bool thisbw = thisStartPos & 0x80000000;
				bool otherbw = otherStartPos & 0x80000000;
				thisStartPos &= 0x7FFFFFFF;
				otherStartPos &= 0x7FFFFFFF;
				thisEndPos &= 0x7FFFFFFF;
				otherEndPos &= 0x7FFFFFFF;
				assert(thisEndPos > thisStartPos);
				assert(otherEndPos > otherStartPos);
				callback(i, thisStartPos, thisEndPos, !thisbw, otherread, otherStartPos, otherEndPos, !otherbw);
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
	IterationInfo iterateMatchReadPairs(F callback) const
	{
		size_t currentRead = 0;
		std::unordered_set<size_t> matches;
		auto result = iterateMatches(true, [&currentRead, &matches, callback](size_t left, size_t leftStartPos, size_t leftEndPos, bool leftFw, size_t right, size_t rightStartPos, size_t rightEndPos, bool rightFw)
		{
			if (left != currentRead)
			{
				if (matches.size() > 0) callback(currentRead, matches);
				currentRead = left;
				matches.clear();
			}
			matches.insert(right);
		});
		if (matches.size() > 0) callback(currentRead, matches);
		return result;
	}
	template <typename F>
	IterationInfo iterateMatchChains(const std::vector<size_t>& rawReadLengths, F callback) const
	{
		size_t currentLeftRead = std::numeric_limits<size_t>::max();
		size_t currentRightRead = std::numeric_limits<size_t>::max();
		std::vector<std::tuple<size_t, size_t, size_t, size_t>> currentFwFwMatches;
		std::vector<std::tuple<size_t, size_t, size_t, size_t>> currentFwBwMatches;
		std::vector<std::tuple<size_t, size_t, size_t, size_t>> currentBwFwMatches;
		std::vector<std::tuple<size_t, size_t, size_t, size_t>> currentBwBwMatches;
		size_t totalChains = 0;
		auto result = iterateMatches(false, [this, &totalChains, &currentLeftRead, &currentRightRead, &currentFwFwMatches, &currentFwBwMatches, &currentBwFwMatches, &currentBwBwMatches, &rawReadLengths, callback](size_t leftread, size_t leftStartPos, size_t leftEndPos, bool leftFw, size_t rightread, size_t rightStartPos, size_t rightEndPos, bool rightFw)
		{
			if (leftread != currentLeftRead || rightread != currentRightRead)
			{
				totalChains += iterateChains(currentLeftRead, currentRightRead, currentFwFwMatches, true, true, callback);
				totalChains += iterateChains(currentLeftRead, currentRightRead, currentFwBwMatches, true, false, callback);
				totalChains += iterateChains(currentLeftRead, currentRightRead, currentBwFwMatches, false, true, callback);
				totalChains += iterateChains(currentLeftRead, currentRightRead, currentBwBwMatches, false, false, callback);
				currentFwFwMatches.clear();
				currentBwFwMatches.clear();
				currentFwBwMatches.clear();
				currentBwBwMatches.clear();
				currentLeftRead = leftread;
				currentRightRead = rightread;
			}
			assert(leftEndPos > leftStartPos);
			assert(rightEndPos > rightStartPos);
			if (rightFw)
			{
				if (leftFw)
				{
					currentFwFwMatches.emplace_back(leftStartPos, leftEndPos, rightStartPos, rightEndPos);
				}
				else
				{
					currentBwFwMatches.emplace_back(leftStartPos, leftEndPos, rightStartPos, rightEndPos);
				}
			}
			else
			{
				if (leftFw)
				{
					currentFwBwMatches.emplace_back(leftStartPos, leftEndPos, rightStartPos, rightEndPos);
				}
				else
				{
					currentBwBwMatches.emplace_back(leftStartPos, leftEndPos, rightStartPos, rightEndPos);
				}
			}
		});
		totalChains += iterateChains(currentLeftRead, currentRightRead, currentFwFwMatches, true, true, callback);
		totalChains += iterateChains(currentLeftRead, currentRightRead, currentFwBwMatches, true, false, callback);
		totalChains += iterateChains(currentLeftRead, currentRightRead, currentBwFwMatches, false, true, callback);
		totalChains += iterateChains(currentLeftRead, currentRightRead, currentBwBwMatches, false, false, callback);
		result.readChainMatches = totalChains;
		return result;
	}
	template <typename F>
	size_t iterateChains(size_t leftread, size_t rightread, std::vector<std::tuple<size_t, size_t, size_t, size_t>>& matches, bool leftFw, bool rightFw, F callback) const
	{
		const size_t maxChainDiagonalDifference = 500;
		if (matches.size() == 0) return 0;
		size_t result = 0;
		// sort by diagonal
		std::sort(matches.begin(), matches.end(), [](std::tuple<size_t, size_t, size_t, size_t> left, std::tuple<size_t, size_t, size_t, size_t> right) { return (int64_t)std::get<2>(left) - (int64_t)std::get<0>(left) < (int64_t)std::get<2>(right) - (int64_t)std::get<0>(right); });
		size_t leftStart = std::get<0>(matches[0]);
		size_t leftEnd = std::get<1>(matches[0]);
		size_t rightStart = std::get<2>(matches[0]);
		size_t rightEnd = std::get<3>(matches[0]);
		for (size_t i = 1; i < matches.size(); i++)
		{
			if ((int64_t)std::get<2>(matches[i]) - (int64_t)std::get<0>(matches[i]) > (int64_t)std::get<2>(matches[i-1]) - (int64_t)std::get<0>(matches[i-1]) + (int64_t)maxChainDiagonalDifference)
			{
				callback(leftread, leftStart, leftEnd, leftFw, rightread, rightStart, rightEnd, rightFw);
				result += 1;
				leftStart = std::get<0>(matches[i]);
				leftEnd = std::get<1>(matches[i]);
				rightStart = std::get<2>(matches[i]);
				rightEnd = std::get<3>(matches[i]);
			}
			else
			{
				leftStart = std::min(leftStart, std::get<0>(matches[i]));
				leftEnd = std::max(leftEnd, std::get<1>(matches[i]));
				rightStart = std::min(rightStart, std::get<2>(matches[i]));
				rightEnd = std::max(rightEnd, std::get<3>(matches[i]));
			}
		}
		assert(leftEnd > leftStart);
		assert(rightEnd > rightStart);
		callback(leftread, leftStart, leftEnd, leftFw, rightread, rightStart, rightEnd, rightFw);
		result += 1;
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
	void addMatchesFromReadOneWay(uint32_t readKey, std::mutex& indexMutex, const std::string& readSequence, bool fw);
	size_t k;
	size_t numWindows;
	size_t windowSize;
	size_t numReads;
	ReadIdContainer idContainer;
};

#endif
