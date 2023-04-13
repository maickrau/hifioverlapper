#ifndef UnitigKmerCorrector_h
#define UnitigKmerCorrector_h

#include <vector>
#include <string>
#include <tuple>
#include <map>
#include <phmap.h>
#include "UnitigStorage.h"
#include "ReadHelper.h"
#include "VectorWithDirection.h"

class UnitigKmerCorrector
{
	class Read
	{
	public:
		std::string name;
		std::vector<std::pair<size_t, bool>> unitigPath;
		size_t leftClip;
		size_t rightClip;
		std::string leftHanger;
		std::string rightHanger;
	};
public:
	class LocalGraph
	{
	public:
		std::vector<size_t> globalToLocal;
		std::vector<size_t> localToGlobal;
		std::vector<bool> safeNode;
		std::vector<bool> ambiguousNode;
		VectorWithDirection<std::vector<std::pair<size_t, bool>>> safeEdges;
		VectorWithDirection<std::vector<std::pair<size_t, bool>>> ambiguousEdges;
		size_t size() const;
	};
	UnitigKmerCorrector(size_t k);
	void build(const ReadpartIterator& iterator);
	template <typename F>
	void iterateRawSequences(F callback) const
	{
		for (size_t i = 0; i < reads.size(); i++)
		{
			std::string corrected = getRawSequence(i);
			callback(i, reads[i].name, corrected);
		}
	}
	std::pair<std::string, bool> getCorrectedSequence(size_t readIndex, const std::vector<size_t>& context, size_t minAmbiguousCoverage, size_t minSafeCoverage) const;
	std::vector<size_t> filterDifferentHaplotypesOut(size_t readIndex, const std::vector<size_t>& context, size_t minAmbiguousCoverage, size_t minSafeCoverage) const;
	std::string getRawSequence(size_t index) const;
	size_t numReads() const;
	const std::string& getName(size_t index) const;
private:
	std::vector<std::vector<std::pair<size_t, bool>>> getPossiblePaths(const std::pair<size_t, bool> start, const std::pair<size_t, bool> end, const phmap::flat_hash_set<std::pair<size_t, bool>>& bwVisited, const std::vector<bool>& allowedNodes, const VectorWithDirection<std::vector<std::pair<size_t, bool>>>& allowedEdges, size_t maxPaths) const;
	void addPathsRecursive(std::vector<std::vector<std::pair<size_t, bool>>>& result, std::vector<std::pair<size_t, bool>>& currentPath, const std::pair<size_t, bool> pos, const std::pair<size_t, bool> end, const phmap::flat_hash_set<std::pair<size_t, bool>>& bwVisited, const std::vector<bool>& allowedNodes, const VectorWithDirection<std::vector<std::pair<size_t, bool>>>& allowedEdges, size_t maxPaths) const;
	void clearLocalGraph(LocalGraph& graph) const;
	void assignLocalGraph(LocalGraph& graph, const std::vector<size_t>& context, size_t minAmbiguousCoverage, size_t minSafeCoverage) const;
	std::vector<std::pair<size_t, bool>> getUniqueReplacementPath(std::pair<size_t, bool> start, std::pair<size_t, bool> end, const std::vector<bool>& allowedNodes, const VectorWithDirection<std::vector<std::pair<size_t, bool>>>& allowedEdges, const std::vector<size_t>& localToGlobal, size_t maxLength) const;
	void assignReadsToAlleles(const std::vector<size_t>& context, const std::vector<size_t>& localToGlobal, std::vector<std::vector<std::vector<size_t>>>& result, const std::vector<std::vector<size_t>>& alleles) const;
	void forbidOtherHaplotypes(phmap::flat_hash_set<size_t>& forbiddenReads, size_t readIndex, const std::vector<std::vector<size_t>>& leftAlleles, const std::vector<std::vector<size_t>>& rightAlleles) const;
	size_t kmerSize;
	UnitigStorage unitigs;
	std::vector<Read> reads;
};

#endif
