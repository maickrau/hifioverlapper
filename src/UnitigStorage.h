#ifndef UnitigStorage_h
#define UnitigStorage_h

#include <tuple>
#include <vector>
#include <string>
#include <mutex>
#include <phmap.h>
#include "MBGCommon.h"
#include "VectorWithDirection.h"
#include "TwobitString.h"

class UnitigStorage
{
public:
	UnitigStorage(size_t k);
	std::pair<size_t, bool> getNode(HashType hash);
	std::pair<size_t, bool> getNodeOrNull(HashType hash) const;
	void addEdge(std::pair<size_t, bool> from, std::pair<size_t, bool> to, size_t overlap);
	std::tuple<size_t, size_t, std::vector<std::pair<size_t, bool>>> getPath(const std::vector<HashType>& hashes) const;
	std::string getSequenceBpClip(const std::vector<std::pair<size_t, bool>>& path, size_t leftClipBp, size_t rightClipBp) const;
	std::string getSequence(const std::vector<std::pair<size_t, bool>>& path, size_t leftClipKmer, size_t rightClipKmer) const;
	size_t unitigSize(size_t unitig) const;
	size_t unitigMinusEdgeLength(const std::pair<size_t, bool> fromUnitig, const std::pair<size_t, bool> toUnitig) const;
	void buildUnitigGraph();
	void addKmerSequence(std::pair<size_t, bool> kmer, const std::string& seq, size_t start, size_t end);
	void finalizeSequences();
	size_t numHashes() const;
	size_t numUnitigs() const;
	size_t totalBps() const;
private:
	void beginUnitig(std::pair<size_t, bool> start);
	std::pair<size_t, bool> numToPair(size_t num) const;
	size_t pairToNum(std::pair<size_t, bool> p) const;
	size_t kmerSize;
	phmap::flat_hash_map<HashType, uint32_t> hashToNode;
	phmap::flat_hash_map<std::pair<uint32_t, uint32_t>, uint16_t> unitigEdgeOverlap;
	phmap::flat_hash_map<std::pair<uint32_t, uint32_t>, uint16_t> kmerOverlap;
	VectorWithDirection<uint32_t> uniqueEdge;
	std::vector<std::pair<uint32_t, uint32_t>> kmerPositionInUnitig;
	std::vector<std::vector<uint32_t>> kmerBpOffsetInsideUnitig;
	std::vector<TwobitString> unitigSequences;
	std::vector<bool> kmerSequenceLoaded;
	std::vector<uint32_t> unitigLength;
	std::mutex addMutex;
};

#endif
