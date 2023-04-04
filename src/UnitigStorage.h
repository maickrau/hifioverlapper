#ifndef UnitigStorage_h
#define UnitigStorage_h

#include <tuple>
#include <vector>
#include <string>
#include <mutex>
#include <phmap.h>
#include "MBGCommon.h"
#include "VectorWithDirection.h"

class UnitigStorage
{
public:
	UnitigStorage(size_t k);
	std::pair<size_t, bool> getNode(HashType hash);
	std::pair<size_t, bool> getNodeOrNull(HashType hash) const;
	void addEdge(std::pair<size_t, bool> from, std::pair<size_t, bool> to, size_t overlap);
	std::tuple<size_t, size_t, std::vector<std::pair<size_t, bool>>> getPath(const std::vector<HashType>& hashes) const;
	std::string getSequence(const std::vector<std::pair<size_t, bool>>& path, size_t leftClip, size_t rightClip) const;
	void buildUnitigGraph();
	void addKmerSequence(std::pair<size_t, bool> kmer, const std::string& seq, size_t start, size_t end);
	void finalizeSequences();
	size_t numHashes() const;
	size_t numUnitigs() const;
private:
	void beginUnitig(std::pair<size_t, bool> start);
	std::pair<size_t, bool> numToPair(size_t num) const;
	size_t pairToNum(std::pair<size_t, bool> p) const;
	size_t kmerSize;
	phmap::flat_hash_map<HashType, size_t> hashToNode;
	phmap::flat_hash_map<std::pair<size_t, size_t>, size_t> unitigEdgeOverlap;
	phmap::flat_hash_map<std::pair<size_t, size_t>, size_t> kmerOverlap;
	VectorWithDirection<size_t> uniqueEdge;
	std::vector<std::pair<size_t, size_t>> kmerPositionInUnitig;
	std::vector<std::vector<size_t>> kmerBpOffsetInsideUnitig;
	std::vector<std::string> unitigSequences;
	std::vector<bool> kmerSequenceLoaded;
	std::vector<size_t> unitigLength;
	std::mutex addMutex;
};

#endif
