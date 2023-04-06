#ifndef UnitigKmerCorrector_h
#define UnitigKmerCorrector_h

#include <vector>
#include <string>
#include <tuple>
#include <map>
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
	std::string getRawSequence(size_t index) const;
	size_t numReads() const;
	const std::string& getName(size_t index) const;
private:
	std::vector<std::pair<size_t, bool>> getUniqueReplacementPath(std::pair<size_t, bool> start, std::pair<size_t, bool> end, const std::vector<bool>& allowedNodes, const VectorWithDirection<std::vector<std::pair<size_t, bool>>>& allowedEdges, const std::vector<size_t>& localToGlobal, size_t maxLength) const;
	size_t kmerSize;
	UnitigStorage unitigs;
	std::vector<Read> reads;
};

#endif
