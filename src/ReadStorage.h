#ifndef ReadStorage_h
#define ReadStorage_h

#include <string>
#include <vector>
#include <mutex>
#include "HashList.h"
#include "TwobitString.h"
#include "ReadHelper.h"
#include "fastqloader.h"

extern thread_local std::vector<size_t> ReadStorageMemoryIterables;

class ReadStorage
{
public:
	void storeReadsFromFile(const std::string& filename, bool includeSequences);
	template <typename F>
	void iterateReadsFromFile(const std::string& filename, size_t numThreads, bool store, F callback)
	{
		std::mutex nameMutex;
		std::vector<std::string> files { filename };
		iterateReadsMultithreaded(files, numThreads, [this, store, &nameMutex, callback](const ReadInfo& info, const std::string& sequence)
		{
			size_t name = 0;
			{
				std::lock_guard<std::mutex> guard { nameMutex };
				name = names.size();
				names.push_back(info.readName.first);
				rawReadLengths.push_back(sequence.size());
				if (store) sequences.emplace_back(sequence);
			}
			callback(name, sequence);
		});
	}
	template <typename F>
	void iterateReadsFromStorage(F callback) const
	{
		assert(names.size() == sequences.size());
		for (size_t i = 0; i < names.size(); i++)
		{
			callback(i, sequences[i].toString());
		}
	}
	template <typename F>
	void iterateReadsAndHashesFromStorage(F callback) const
	{
		assert(names.size() == sequences.size());
		assert(positions.size() == sequences.size());
		assert(hashes.size() == positions.size());
		if (ReadStorageMemoryIterables.size() > 0)
		{
			for (auto i : ReadStorageMemoryIterables)
			{
				callback(i, sequences[i].toString(), positions[i], hashes[i]);
			}
		}
		else
		{
			for (size_t i = 0; i < names.size(); i++)
			{
				callback(i, sequences[i].toString(), positions[i], hashes[i]);
			}
		}
	}
	template <typename F>
	void iterateKmersFromStorage(F callback) const
	{
		assert(names.size() == sequences.size());
		assert(positions.size() == sequences.size());
		assert(kmers.size() == positions.size());
		if (ReadStorageMemoryIterables.size() > 0)
		{
			for (auto i : ReadStorageMemoryIterables)
			{
				callback(i, positions[i], kmers[i]);
			}
		}
		else
		{
			for (size_t i = 0; i < names.size(); i++)
			{
				callback(i, positions[i], kmers[i]);
			}
		}
	}
	std::pair<std::string, std::string> getRead(size_t i) const;
	const std::vector<std::string>& getNames() const;
	std::string getSequence(size_t i) const;
	const std::vector<size_t>& getRawReadLengths() const;
	const std::vector<size_t>& getPositions(size_t i) const;
	const std::vector<HashType>& getHashes(size_t i) const;
	void buildHashes(const ReadpartIterator& partIterator);
	void buildKmers(const HashList& hashlist);
	size_t size() const;
	void setMemoryIterables(const std::vector<size_t>& iterables);
private:
	std::vector<std::string> names;
	std::vector<size_t> rawReadLengths;
	std::vector<TwobitString> sequences;
	std::vector<std::vector<size_t>> positions;
	std::vector<std::vector<HashType>> hashes;
	std::vector<std::vector<std::pair<size_t, bool>>> kmers;
};

#endif
