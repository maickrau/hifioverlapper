#ifndef MatchIndex_h
#define MatchIndex_h

#include <string>
#include <cstdint>
#include <vector>
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
	void addMatchesFromFile(size_t numThreads, std::string file);
	template <typename F>
	IterationInfo iterateMatches(F callback) const
	{
		IterationInfo result;
		result.numberReads = names.size();
		result.numerWindowChunks = idContainer.size();
		std::vector<std::vector<uint32_t>> hashesPerRead;
		hashesPerRead.resize(names.size());
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
				callback(names[i], names[match]);
			}
		}
		result.readsWithMatch = readsWithMatch;
		result.maxPerChunk = maxPerChunk;
		result.totalMatches = totalMatches;
		return result;
	}
private:
	size_t k;
	size_t numWindows;
	size_t windowSize;
	std::vector<std::string> names;
	ReadIdContainer idContainer;
};

#endif
