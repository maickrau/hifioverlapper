#include "ReadStorage.h"

void ReadStorage::storeReadsFromFile(const std::string& filename, bool includeSequences)
{
	FastQ::streamFastqFromFile(filename, false, [this](FastQ& read)
	{
		names.push_back(read.seq_id);
		sequences.push_back(read.sequence);
	});
}

const std::vector<std::string>& ReadStorage::getNames() const
{
	return names;
}

const std::vector<size_t>& ReadStorage::getRawReadLengths() const
{
	return rawReadLengths;
}