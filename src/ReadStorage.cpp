#include "ReadStorage.h"

thread_local std::vector<size_t> ReadStorageMemoryIterables;

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

std::pair<std::string, std::string> ReadStorage::getRead(size_t i) const
{
	return std::make_pair(names[i], sequences[i]);
}

size_t ReadStorage::size() const
{
	return names.size();
}

const std::vector<size_t>& ReadStorage::getPositions(size_t i) const
{
	return positions[i];
}

const std::vector<HashType>& ReadStorage::getHashes(size_t i) const
{
	return hashes[i];
}

void ReadStorage::buildHashes(const ReadpartIterator& partIterator)
{
	for (size_t i = 0; i < sequences.size(); i++)
	{
		partIterator.iterateHashesOfRead(names[i], sequences[i], [this](const ReadInfo& read, const SequenceCharType& seq, const SequenceLengthType& poses, const std::string& rawSeq, const std::vector<size_t>& positions, const std::vector<HashType>& hashes)
		{
			this->positions.push_back(positions);
			this->hashes.push_back(hashes);
		});
	}
	assert(positions.size() == sequences.size());
	assert(hashes.size() == sequences.size());
}

void ReadStorage::setMemoryIterables(const std::vector<size_t>& iterables)
{
	ReadStorageMemoryIterables = iterables;
}

const std::string& ReadStorage::getSequence(size_t i) const
{
	return sequences[i];
}
