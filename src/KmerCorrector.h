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
	void buildGraph(const ReadpartIterator& iterator, size_t numThreads);
	void buildGraph(const ReadStorage& iterator, size_t numThreads);
	std::pair<std::string, bool> getCorrectedSequence(const std::string& rawSeq, const std::vector<size_t>& positions, const std::vector<HashType>& hashes) const;
private:
	void loadReadsAsHashesAndKmerSequencesMultithread(HashList& result, const size_t kmerSize, const ReadStorage& storage, const size_t numThreads);
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
	std::vector<bool> hasFwCoverage;
	std::vector<bool> hasBwCoverage;
	std::vector<bool> removedHomopolymerError;
};

#endif
