#ifndef KmerCorrector_h
#define KmerCorrector_h

#include <string>
#include <tuple>
#include "ReadStorage.h"
#include "SparseEdgeContainer.h"
#include "HashList.h"
#include "RankBitvector.h"
#include "ReadHelper.h"

class ConcatenatedStringStorage
{
public:
	void resize(size_t numItems, size_t k);
	std::string getSequence(size_t index) const;
	void setSequence(size_t index, std::string seq);
	size_t size() const;
private:
	std::string sequence;
	size_t k;
};

class KmerCorrector
{
public:
	KmerCorrector(size_t kmerSize, size_t minSolidCoverage, size_t minAmbiguousCoverage);
	void buildGraph(const ReadpartIterator& iterator);
	void buildGraph(const ReadStorage& iterator);
	void filterToReadCoverage(const ReadStorage& iterator);
	void initializeThreadCoverageFullGraph();
	const HashList& getHashlist() const;
	std::pair<std::string, bool> getCorrectedSequence(const std::string& rawSeq, const std::vector<size_t>& positions, const std::vector<HashType>& hashes) const;
private:
	template <typename F>
	void loadReadsAsHashesAndKmerSequencesMultithread(HashList& result, const size_t kmerSize, F readGetter)
	{
		hasSequence.resize(reads.size());
		size_t totalWithSequence = 0;
		std::vector<std::pair<size_t, std::string>> sequences;
		std::mutex sequenceMutex;
		std::vector<bool> hasFwCoverage;
		std::vector<bool> hasBwCoverage;
		readGetter([this, &result, kmerSize, &sequenceMutex, &totalWithSequence, &hasFwCoverage, &hasBwCoverage, &sequences](const std::string& rawSeq, const std::vector<size_t>& positions, const std::vector<HashType>& hashes)
		{
			size_t lastMinimizerPosition = std::numeric_limits<size_t>::max();
			std::pair<size_t, bool> last { std::numeric_limits<size_t>::max(), true };
			assert(positions.size() == hashes.size());
			for (size_t i = 0; i < positions.size(); i++)
			{
				const auto pos = positions[i];
				const HashType fwHash = hashes[i];
				assert(last.first == std::numeric_limits<size_t>::max() || pos - lastMinimizerPosition <= kmerSize);
				std::pair<size_t, bool> current = result.addNode(fwHash);
				size_t overlap = lastMinimizerPosition + kmerSize - pos;
				assert(pos+kmerSize <= rawSeq.size());
				assert(lastMinimizerPosition == std::numeric_limits<size_t>::max() || pos - lastMinimizerPosition < kmerSize);
				if (last.first != std::numeric_limits<size_t>::max())
				{
					assert(lastMinimizerPosition + kmerSize >= pos);
					result.addSequenceOverlap(last, current, overlap);
					auto pair = canon(last, current);
					result.addEdgeCoverage(pair.first, pair.second);
				}
				lastMinimizerPosition = pos;
				last = current;
				{
					std::lock_guard lock { sequenceMutex };
					while (current.first >= hasSequence.size())
					{
						hasSequence.push_back(false);
						hasFwCoverage.push_back(false);
						hasBwCoverage.push_back(false);
					}
					assert(current.first < hasSequence.size());
					if (!hasSequence.get(current.first))
					{
						if (current.second)
						{
							hasFwCoverage[current.first] = true;
						}
						else
						{
							hasBwCoverage[current.first] = true;
						}
						if (hasFwCoverage[current.first] && hasBwCoverage[current.first] && result.coverage.get(current.first) >= minAmbiguousCoverage)
						{
							totalWithSequence += 1;
							hasSequence.set(current.first, true);
							std::string seqHere = rawSeq.substr(pos, kmerSize);
							assert(seqHere.size() == kmerSize);
							if (!current.second) seqHere = revCompRaw(seqHere);
							sequences.emplace_back(current.first, seqHere);
						}
					}
				}
			};
		});
		assert(hasSequence.size() == result.size());
		hasSequence.buildRanks();
		kmerSequences.resize(totalWithSequence, kmerSize);
		assert(sequences.size() == totalWithSequence);
		for (size_t i = 0; i < sequences.size(); i++)
		{
			size_t sequenceIndex = hasSequence.getRank(sequences[i].first);
			kmerSequences.setSequence(sequenceIndex, sequences[i].second);
		}
	}
	size_t getCoverage(size_t kmer) const;
	std::pair<size_t, bool> findBubble(std::pair<size_t, bool> start) const;
	void forbidPathNodes(const std::vector<std::pair<size_t, bool>>& path);
	void allowPathNodes(const std::vector<std::pair<size_t, bool>>& path);
	size_t getPathCoverage(const std::vector<std::pair<size_t, bool>>& path) const;
	std::pair<std::string, std::vector<size_t>> getHomopolymerCompressedPathSequence(const std::vector<std::pair<size_t, bool>>& path) const;
	void forbidHomopolymerAlleles(const std::pair<size_t, bool> start, const std::pair<size_t, bool> end);
	void enumeratePathsRecursion(std::vector<std::vector<std::pair<size_t, bool>>>& result, std::vector<std::pair<size_t, bool>>& currentPath, const std::pair<size_t, bool> end, const size_t maxCount) const;
	std::vector<std::vector<std::pair<size_t, bool>>> enumeratePaths(const std::pair<size_t, bool> start, const std::pair<size_t, bool> end, const size_t maxCount) const;
	void forbidHomopolymerErrors();
	size_t kmerSize;
	size_t minSolidCoverage;
	size_t minAmbiguousCoverage;
	SparseEdgeContainer edges;
	ConcatenatedStringStorage kmerSequences;
	RankBitvector hasSequence;
	HashList reads;
	std::vector<bool> removedHomopolymerError;
};

#endif
