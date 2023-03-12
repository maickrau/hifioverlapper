#ifndef KmerCorrector_h
#define KmerCorrector_h

#include <string>
#include <tuple>
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
	std::pair<std::string, bool> getCorrectedSequence(const std::string& rawSeq, const std::vector<size_t>& positions, const std::vector<HashType>& hashes) const;
private:
	size_t kmerSize;
	size_t minSolidCoverage;
	size_t minAmbiguousCoverage;
	SparseEdgeContainer edges;
	ConcatenatedStringStorage kmerSequences;
	RankBitvector hasSequence;
	HashList reads;
};

#endif
