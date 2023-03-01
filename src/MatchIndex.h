#ifndef MatchIndex_h
#define MatchIndex_h

#include <string>
#include <cstdint>
#include <vector>
#include <mutex>
#include <phmap.h>

class ReadIdContainer
{
public:
	ReadIdContainer() = default;
	void addNumber(uint64_t key, uint32_t value);
	const std::vector<std::vector<uint32_t>>& getMultiNumbers() const;
	size_t size() const;
private:
	phmap::flat_hash_map<uint64_t, uint32_t> firstNumberOrVectorIndex;
	std::vector<std::vector<uint32_t>> numbers;
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
};

class MatchIndex
{
public:
	MatchIndex(size_t k, size_t numWindows, size_t windowSize);
	void addMatch(uint64_t hash, uint32_t readKey);
	void addMatchesFromRead(uint32_t readKey, std::mutex& indexMutex, const std::string& readSequence);
	template <typename F>
	IterationInfo iterateMatches(F callback) const
	{
		IterationInfo result;
		result.numberReads = numReads;
		result.numerWindowChunks = idContainer.size();
		std::vector<std::vector<uint32_t>> hashesPerRead;
		hashesPerRead.resize(numReads);
		size_t maxPerChunk = 0;
		const auto& numbers = idContainer.getMultiNumbers();
		size_t uniqueChunks = idContainer.size() - numbers.size();
		size_t totalReadChunkMatches = 0;
		for (size_t i = 0; i < numbers.size(); i++)
		{
			totalReadChunkMatches += numbers[i].size();
			maxPerChunk = std::max(maxPerChunk, numbers[i].size());
			for (auto readname : numbers[i])
			{
				hashesPerRead[readname].push_back(i);
			}
		}
		result.totalReadChunkMatches = totalReadChunkMatches;
		result.uniqueChunks = uniqueChunks;
		size_t readsWithMatch = 0;
		size_t totalMatches = 0;
		for (size_t i = 0; i < hashesPerRead.size(); i++)
		{
			phmap::flat_hash_set<size_t> matches;
			for (auto hash : hashesPerRead[i])
			{
				for (auto read : numbers[hash])
				{
					if (read == i) continue;
					if (read < i) continue;
					matches.emplace(read);
				}
			}
			if (matches.size() > 0) readsWithMatch += 1;
			totalMatches += matches.size();
			for (auto match : matches)
			{
				callback(i, match);
			}
		}
		result.readsWithMatch = readsWithMatch;
		result.maxPerChunk = maxPerChunk;
		result.totalMatches = totalMatches;
		return result;
	}
	template <typename F>
	IterationInfo iterateMatchNames(const std::vector<std::string>& names, F callback) const
	{
		return iterateMatches([&names, callback](size_t left, size_t right)
		{
			callback(names[left], names[right]);
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
