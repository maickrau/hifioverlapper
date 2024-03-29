#ifndef MatchIndex_h
#define MatchIndex_h

#include <string>
#include <cstdint>
#include <vector>
#include <mutex>
#include <tuple>
#include <thread>
#include <unordered_set>
#include <phmap.h>
#include "ReadMatchposStorage.h"
#include "ReadHelper.h"
#include "MinimizerIterator.h"

class ReadIdContainer
{
public:
	ReadIdContainer() = default;
	void clearConstructionVariablesAndCompact();
	void addNumber(uint64_t key, __uint128_t value);
	const std::vector<ReadMatchposStorage>& getMultiNumbers() const;
	__uint128_t getFirstNumberOrBucketIndex(uint64_t key) const;
	size_t size() const;
	static constexpr __uint128_t FirstBit = (__uint128_t)1 << (__uint128_t)127;
private:
	phmap::parallel_flat_hash_map<uint64_t, __uint128_t> firstNumberOrVectorIndex;
	std::vector<ReadMatchposStorage> numbers;
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

class Match
{
public:
	Match() = default;
	Match(size_t leftStartPos, size_t leftEndPos, bool leftFw, size_t rightStartPos, size_t rightEndPos, bool rightFw);
	size_t leftStartPos;
	size_t leftEndPos;
	size_t rightStartPos;
	size_t rightEndPos;
	bool leftFw;
	bool rightFw;
};

class MatchIndex
{
public:
	MatchIndex(size_t k, size_t numWindows, size_t windowSize);
	void clearConstructionVariablesAndCompact();
	void addMatch(uint64_t hash, __uint128_t readKey);
	void addMatchesFromRead(uint32_t readKey, std::mutex& indexMutex, const std::string& readSequence);
	size_t numWindowChunks() const;
	size_t numUniqueChunks() const;
	template <typename F>
	void iterateWindowChunksFromRead(const std::string& readSequence, F callback) const
	{
		iterateWindowChunksFromReadOneWay(readSequence, callback);
		auto rev = MBG::revCompRaw(readSequence);
		iterateWindowChunksFromReadOneWay(rev, [callback](uint32_t startPos, uint32_t endPos, uint64_t hash)
		{
			startPos += 0x80000000;
			endPos += 0x80000000;
			callback(startPos, endPos, hash);
		});
	}
	template <typename F>
	void iterateWindowChunksFromReadOneWay(const std::string& readSequence, F callback) const
	{
		MBG::ErrorMasking errorMasking = MBG::ErrorMasking::Microsatellite;
		std::vector<std::string> readFiles { };
		MBG::ReadpartIterator partIterator { 31, 1, errorMasking, 1, readFiles, false, "" };
		partIterator.iteratePartsOfRead("", readSequence, [this, callback](const MBG::ReadInfo& read, const MBG::SequenceCharType& seq, const MBG::SequenceLengthType& poses, const std::string& raw)
		{
			if (seq.size() < numWindows * windowSize + k) return;
			phmap::flat_hash_set<std::tuple<uint64_t, uint32_t, uint32_t>> hashesHere;
			size_t rawReadLength = raw.size();
			size_t compressedReadLength = seq.size();
			iterateWindowchunks(seq, k, numWindows, windowSize, [this, &hashesHere, &poses, rawReadLength, compressedReadLength](const std::vector<uint64_t>& hashes, const size_t startPos, const size_t endPos)
			{
				uint64_t totalhash = 0;
				for (auto hash : hashes)
				{
					totalhash *= 3;
					totalhash += hash;
				}
				assert(endPos > startPos);
				assert(startPos < poses.size());
				assert(endPos < poses.size());
				assert(endPos+1 < poses.size());
				uint32_t realStartPos = poses[startPos];
				uint32_t realEndPos = poses[endPos+1]-1;
				assert(realEndPos > realStartPos);
				assert((realEndPos & 0x7FFFFFFF) < rawReadLength);
				assert(realStartPos < 0x80000000);
				assert(realEndPos < 0x80000000);
				hashesHere.emplace(totalhash, realStartPos, realEndPos);
			});
			for (auto t : hashesHere)
			{
				callback(std::get<1>(t), std::get<2>(t), std::get<0>(t));
			}
		});
	}
	template <typename F>
	void iterateMatchesOneRead(const size_t maxLengthDifference, const std::string& readSequence, F callback) const
	{
		const auto& numbers = idContainer.getMultiNumbers();
		phmap::flat_hash_map<uint32_t, std::vector<Match>> matchesPerRead;
		iterateWindowChunksFromRead(readSequence, [this, &numbers, &matchesPerRead, maxLengthDifference](uint32_t startpos, uint32_t endpos, uint64_t hash)
		{
			__uint128_t info = idContainer.getFirstNumberOrBucketIndex(hash);
			if (info == 0) return;
			bool thisFw = (startpos & 0x80000000) == 0;
			if (info & ReadIdContainer::FirstBit)
			{
				uint32_t read = (uint32_t)(info >> (__uint128_t)64);
				uint32_t otherStartPos = (uint32_t)(info >> (__uint128_t)32);
				uint32_t otherEndPos = (uint32_t)(info >> (__uint128_t)0);
				if ((endpos-startpos) > (otherEndPos-otherStartPos) + maxLengthDifference) return;
				if ((otherEndPos-otherStartPos) > (endpos-startpos) + maxLengthDifference) return;
				bool otherFw = (otherStartPos & 0x80000000) == 0;
				matchesPerRead[read].emplace_back(startpos & 0x7FFFFFFF, endpos & 0x7FFFFFFF, thisFw, otherStartPos & 0x7FFFFFFF, otherEndPos & 0x7FFFFFFF, otherFw);
			}
			else
			{
				for (std::tuple<uint32_t, uint32_t, uint32_t> readkey : numbers[info])
				{
					uint32_t read = std::get<0>(readkey);
					uint32_t otherStartPos = std::get<1>(readkey);
					uint32_t otherEndPos = std::get<2>(readkey);
					if ((endpos-startpos) > (otherEndPos-otherStartPos) + maxLengthDifference) continue;
					if ((otherEndPos-otherStartPos) > (endpos-startpos) + maxLengthDifference) continue;
					bool otherFw = (otherStartPos & 0x80000000) == 0;
					matchesPerRead[read].emplace_back(startpos & 0x7FFFFFFF, endpos & 0x7FFFFFFF, thisFw, otherStartPos & 0x7FFFFFFF, otherEndPos & 0x7FFFFFFF, otherFw);
				}
			}
		});
		for (const auto& pair : matchesPerRead)
		{
			callback(pair.first, pair.second);
		}
	}
	template <typename F>
	void iterateChunks(F callback) const
	{
		const auto& numbers = idContainer.getMultiNumbers();
		for (size_t i = 0; i < numbers.size(); i++)
		{
			callback(i, numbers[i]);
		}
	}
	template <typename F>
	IterationInfo iterateMatches(const size_t numThreads, const size_t minCoverage, const size_t maxCoverage, const size_t maxLengthDifference, bool alsoSmaller, const std::vector<bool>& useThese, F callback) const
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
			if (numbers[i].size() < minCoverage) continue;
			if (numbers[i].size() > maxCoverage) continue;
			if (i < useThese.size() && !useThese[i]) continue;
			totalReadChunkMatches += numbers[i].size();
			maxPerChunk = std::max(maxPerChunk, numbers[i].size());
			phmap::flat_hash_set<uint32_t> readsHere;
			for (std::tuple<uint32_t, uint32_t, uint32_t> readkey : numbers[i])
			{
				uint32_t read = std::get<0>(readkey);
				readsHere.insert(read);
			}
			for (uint32_t read : readsHere)
			{
				hashesPerRead[read].emplace_back(i);
			}
		}
		result.totalReadChunkMatches = totalReadChunkMatches;
		result.uniqueChunks = uniqueChunks;
		size_t readsWithMatch = 0;
		size_t totalMatches = 0;
		size_t readPairMatches = 0;
		std::vector<std::thread> threads;
		size_t nextReadIndex = 0;
		std::mutex readIndexMutex;
		for (size_t threadi = 0; threadi < numThreads; threadi++)
		{
			threads.emplace_back([this, &nextReadIndex, &readIndexMutex, &hashesPerRead, &readPairMatches, &totalMatches, &readsWithMatch, &numbers, alsoSmaller, maxLengthDifference, callback]()
			{
				size_t totalMatchesThisThread = 0;
				size_t readsWithMatchThisThread = 0;
				size_t readPairMatchesThisThread = 0;
				while (true)
				{
					size_t i = 0;
					{
						std::lock_guard<std::mutex> lock { readIndexMutex };
						i = nextReadIndex;
						nextReadIndex += 1;
					}
					if (i >= hashesPerRead.size()) break;
					phmap::flat_hash_map<uint32_t, std::vector<Match>> matchesPerRead;
					bool hasMatch = false;
					for (size_t hash : hashesPerRead[i])
					{
						std::vector<std::tuple<uint32_t, uint32_t>> posesForThisRead;
						for (std::tuple<uint32_t, uint32_t, uint32_t> readkey : numbers[hash])
						{
							uint32_t read = std::get<0>(readkey);
							if (read != i) continue;
							uint32_t startPos = std::get<1>(readkey);
							uint32_t endPos = std::get<2>(readkey);
							posesForThisRead.emplace_back(startPos, endPos);
						}
						for (std::tuple<uint32_t, uint32_t, uint32_t> readkey : numbers[hash])
						{
							uint32_t read = std::get<0>(readkey);
							uint32_t otherStartPos = std::get<1>(readkey);
							uint32_t otherEndPos = std::get<2>(readkey);
							bool otherFw = (otherStartPos & 0x80000000) == 0;
							if (read == i) continue;
							hasMatch = true;
							if (!alsoSmaller && read < i) continue;
							for (auto pos : posesForThisRead)
							{
								uint32_t startpos = std::get<0>(pos);
								uint32_t endpos = std::get<1>(pos);
								if ((endpos-startpos) > (otherEndPos-otherStartPos) + maxLengthDifference) continue;
								if ((otherEndPos-otherStartPos) > (endpos-startpos) + maxLengthDifference) continue;
								bool thisFw = (startpos & 0x80000000) == 0;
								matchesPerRead[read].emplace_back(startpos & 0x7FFFFFFF, endpos & 0x7FFFFFFF, thisFw, otherStartPos & 0x7FFFFFFF, otherEndPos & 0x7FFFFFFF, otherFw);
								totalMatchesThisThread += 1;
							}
						}
					}
					if (hasMatch) readsWithMatchThisThread += 1;
					for (const auto& pair : matchesPerRead)
					{
						callback(i, pair.first, pair.second);
					}
					readPairMatchesThisThread += matchesPerRead.size();
				}
				std::lock_guard<std::mutex> lock { readIndexMutex };
				totalMatches += totalMatchesThisThread;
				readsWithMatch += readsWithMatchThisThread;
				readPairMatches += readPairMatchesThisThread;
			});
		}
		for (size_t i = 0; i < numThreads; i++)
		{
			threads[i].join();
		}
		result.readPairMatches = readPairMatches;
		result.readsWithMatch = readsWithMatch;
		result.maxPerChunk = maxPerChunk;
		result.totalMatches = totalMatches;
		return result;
	}
	template <typename F>
	IterationInfo iterateMatchReadPairs(const std::vector<bool>& useThese, F callback) const
	{
		size_t currentRead = 0;
		std::unordered_set<size_t> currentMatches;
		auto result = iterateMatches(1, 2, std::numeric_limits<size_t>::max(), 10000, true, useThese, [&currentRead, &currentMatches, callback](size_t leftread, size_t rightread, const std::vector<Match>& matches)
		{
			if (leftread != currentRead)
			{
				if (currentMatches.size() > 0) callback(currentRead, currentMatches);
				currentRead = leftread;
				currentMatches.clear();
			}
			currentMatches.insert(rightread);
		});
		if (currentMatches.size() > 0) callback(currentRead, currentMatches);
		return result;
	}
	template <typename F>
	IterationInfo iterateMatchChains(const size_t numThreads, const size_t minCoverage, const size_t maxCoverage, const size_t maxLengthDifference, const std::vector<size_t>& rawReadLengths, const std::vector<bool>& useThese, F callback) const
	{
		std::atomic<size_t> totalChains = 0;
		auto result = iterateMatches(numThreads, minCoverage, maxCoverage, maxLengthDifference, false, useThese, [this, &totalChains, &rawReadLengths, callback](size_t leftRead, size_t rightRead, const std::vector<Match>& matches)
		{
			std::vector<std::tuple<size_t, size_t, size_t, size_t>> currentFwMatches;
			std::vector<std::tuple<size_t, size_t, size_t, size_t>> currentBwMatches;
			for (auto match : matches)
			{
				if (!match.leftFw)
				{
					match.leftFw = true;
					std::swap(match.leftStartPos, match.leftEndPos);
					match.leftStartPos = rawReadLengths[leftRead] - 1 - match.leftStartPos; // end-inclusive positions!
					match.leftEndPos = rawReadLengths[leftRead] - 1 - match.leftEndPos; // end-inclusive positions!
					match.rightFw = !match.rightFw;
					std::swap(match.rightStartPos, match.rightEndPos);
					match.rightStartPos = rawReadLengths[rightRead] - 1 - match.rightStartPos; // end-inclusive positions!
					match.rightEndPos = rawReadLengths[rightRead] - 1 - match.rightEndPos; // end-inclusive positions!
				}
				if (match.rightFw)
				{
					currentFwMatches.emplace_back(match.leftStartPos, match.leftEndPos, match.rightStartPos, match.rightEndPos);
				}
				else
				{
					currentBwMatches.emplace_back(match.leftStartPos, match.leftEndPos, match.rightStartPos, match.rightEndPos);
				}
			}
			size_t totalChainsPerThread = 0;
			totalChainsPerThread += iterateChains(leftRead, rightRead, currentFwMatches, true, true, callback);
			totalChainsPerThread += iterateChains(leftRead, rightRead, currentBwMatches, true, false, callback);
			totalChains += totalChainsPerThread;
		});
		result.readChainMatches = totalChains;
		return result;
	}
	template <typename F>
	void iterateMatchChainsOneRead(const size_t maxLengthDifference, const std::vector<size_t>& rawReadLengths, const std::string& readSequence, F callback) const
	{
		std::atomic<size_t> totalChains = 0;
		iterateMatchesOneRead(maxLengthDifference, readSequence, [this, &totalChains, &rawReadLengths, &readSequence, callback](size_t rightRead, const std::vector<Match>& matches)
		{
			std::vector<std::tuple<size_t, size_t, size_t, size_t>> currentFwMatches;
			std::vector<std::tuple<size_t, size_t, size_t, size_t>> currentBwMatches;
			for (auto match : matches)
			{
				if (!match.leftFw)
				{
					match.leftFw = true;
					std::swap(match.leftStartPos, match.leftEndPos);
					match.leftStartPos = readSequence.size() - 1 - match.leftStartPos; // end-inclusive positions!
					match.leftEndPos = readSequence.size() - 1 - match.leftEndPos; // end-inclusive positions!
					match.rightFw = !match.rightFw;
					std::swap(match.rightStartPos, match.rightEndPos);
					match.rightStartPos = rawReadLengths[rightRead] - 1 - match.rightStartPos; // end-inclusive positions!
					match.rightEndPos = rawReadLengths[rightRead] - 1 - match.rightEndPos; // end-inclusive positions!
				}
				if (match.rightFw)
				{
					currentFwMatches.emplace_back(match.leftStartPos, match.leftEndPos, match.rightStartPos, match.rightEndPos);
				}
				else
				{
					currentBwMatches.emplace_back(match.leftStartPos, match.leftEndPos, match.rightStartPos, match.rightEndPos);
				}
			}
			size_t totalChainsPerThread = 0;
			totalChainsPerThread += iterateChains(0, rightRead, currentFwMatches, true, true, callback);
			totalChainsPerThread += iterateChains(0, rightRead, currentBwMatches, true, false, callback);
			totalChains += totalChainsPerThread;
		});
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
	void iterateMatchNamesOneRead(const size_t maxLengthDifference, const std::vector<std::string>& names, const std::vector<size_t>& rawReadLengths, const std::string& readName, const std::string& readSequence, F callback) const
	{
		iterateMatchChainsOneRead(maxLengthDifference, rawReadLengths, readSequence, [&names, &rawReadLengths, &readName, &readSequence, callback](size_t left, size_t leftStart, size_t leftEnd, bool leftFw, size_t right, size_t rightStart, size_t rightEnd, bool rightFw)
		{
			assert(left == 0);
			callback(readName, readSequence.size(), leftStart, leftEnd, leftFw, names[right], rawReadLengths[right], rightStart, rightEnd, rightFw);
		});
	}
	template <typename F>
	IterationInfo iterateMatchNames(const size_t numThreads, const size_t minCoverage, const size_t maxCoverage, const size_t maxLengthDifference, const std::vector<std::string>& names, const std::vector<size_t>& rawReadLengths, F callback) const
	{
		std::vector<bool> useThese;
		return iterateMatchChains(numThreads, minCoverage, maxCoverage, maxLengthDifference, rawReadLengths, useThese, [&names, &rawReadLengths, callback](size_t left, size_t leftStart, size_t leftEnd, bool leftFw, size_t right, size_t rightStart, size_t rightEnd, bool rightFw)
		{
			callback(names[left], rawReadLengths[left], leftStart, leftEnd, leftFw, names[right], rawReadLengths[right], rightStart, rightEnd, rightFw);
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
