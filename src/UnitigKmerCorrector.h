#ifndef UnitigKmerCorrector_h
#define UnitigKmerCorrector_h

#include <vector>
#include <string>
#include <tuple>
#include "UnitigStorage.h"
#include "ReadHelper.h"

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
	void iterateCorrectedSequences(F callback) const
	{
		for (size_t i = 0; i < reads.size(); i++)
		{
			std::string corrected = getCorrected(i);
			callback(i, reads[i].name, corrected);
		}
	}
private:
	std::string getCorrected(size_t index) const;
	size_t kmerSize;
	UnitigStorage unitigs;
	std::vector<Read> reads;
};

#endif
