#include <fstream>
#include <mutex>
#include <iostream>
#include <vector>
#include <string>
#include "ReadStorage.h"
#include "MatchIndex.h"

int main(int argc, char** argv)
{
	size_t numThreads = std::stoi(argv[1]);
	size_t k = std::stoull(argv[2]);
	size_t numWindows = std::stoull(argv[3]);
	size_t windowSize = std::stoull(argv[4]);
	std::vector<std::string> readFiles;
	for (size_t i = 5; i < argc; i++)
	{
		readFiles.emplace_back(argv[i]);
	}
	size_t numReads = 0;
	size_t numHashes = 0;
	size_t numIndexedHashes = 0;
	phmap::flat_hash_set<uint64_t> seenOnce;
	phmap::flat_hash_map<uint64_t, uint32_t> hashToIndex;
	std::ofstream tmpCollectionFile { "index.tmp", std::ofstream::binary };
	MatchIndex matchIndex { k, numWindows, windowSize };
	ReadStorage storage;
	for (auto file : readFiles)
	{
		std::mutex indexMutex;
		storage.iterateReadsFromFile(file, numThreads, false, [&matchIndex, &indexMutex, &tmpCollectionFile, &seenOnce, &hashToIndex, &numReads, &numHashes, &numIndexedHashes](size_t readName, const std::string& sequence)
		{
			std::vector<std::tuple<uint64_t, uint32_t, uint32_t>> hashes;
			phmap::flat_hash_set<uint64_t> hashesHere;
			matchIndex.iterateWindowChunksFromRead(sequence, [&tmpCollectionFile, &seenOnce, &hashToIndex, &numReads, &numHashes, &numIndexedHashes, &hashesHere, &hashes](uint32_t startPos, uint32_t endPos, uint64_t hash)
			{
				hashesHere.emplace(hash);
				hashes.emplace_back(hash, startPos, endPos);
			});
			std::lock_guard<std::mutex> lock { indexMutex };
			numReads += 1;
			numHashes += hashes.size();
			for (auto hash : hashesHere)
			{
				if (seenOnce.count(hash) == 0)
				{
					seenOnce.insert(hash);
				}
				else
				{
					if (hashToIndex.count(hash) == 0)
					{
						size_t index = hashToIndex.size();
						hashToIndex[hash] = index;
					}
				}
			}
			if (hashes.size() == 0) return;
			uint32_t read = readName;
			uint32_t count = hashes.size();
			tmpCollectionFile.write((char*)&read, 4);
			tmpCollectionFile.write((char*)&count, 4);
			for (auto t : hashes)
			{
				tmpCollectionFile.write((char*)&std::get<0>(t), 8);
				tmpCollectionFile.write((char*)&std::get<1>(t), 4);
				tmpCollectionFile.write((char*)&std::get<2>(t), 4);
			}
		});
	}
	const std::vector<size_t>& readLengths = storage.getRawReadLengths();
	const std::vector<std::string>& readNames = storage.getNames();
	tmpCollectionFile.close();
	std::cerr << numReads << " reads" << std::endl;
	std::cerr << seenOnce.size() << " distinct hashes" << std::endl;
	std::cerr << hashToIndex.size() << " indexed hashes" << std::endl;
	std::cerr << numHashes << " total positions" << std::endl;
	{
		decltype(seenOnce) tmp;
		std::swap(tmp, seenOnce);
	}
	{
		std::ofstream indexMetadata { "index.metadata.index", std::ofstream::binary };
		size_t countHashes = hashToIndex.size();
		size_t countReads = readLengths.size();
		indexMetadata.write((char*)&countHashes, sizeof(size_t));
		indexMetadata.write((char*)&countReads, sizeof(size_t));
		for (size_t read = 0; read < countReads; read++)
		{
			size_t nameLength = readNames[read].size();
			size_t readLength = readLengths[read];
			indexMetadata.write((char*)&readLength, sizeof(size_t));
			indexMetadata.write((char*)&nameLength, sizeof(size_t));
			indexMetadata.write((char*)readNames[read].data(), nameLength);
		}
	}
	std::ofstream indexFile { "index.positions.index", std::ofstream::binary };
	std::ifstream redoHashes { "index.tmp", std::ifstream::binary };
	while (redoHashes.good())
	{
		uint32_t read;
		uint32_t count;
		redoHashes.read((char*)&read, 4);
		redoHashes.read((char*)&count, 4);
		if (!redoHashes.good()) break;
		std::vector<std::tuple<uint32_t, uint32_t, uint32_t>> hashes;
		for (size_t i = 0; i < count; i++)
		{
			uint64_t hash;
			uint32_t startPos, endPos;
			redoHashes.read((char*)&hash, 8);
			redoHashes.read((char*)&startPos, 4);
			redoHashes.read((char*)&endPos, 4);
			if (hashToIndex.count(hash) == 0) continue;
			hashes.emplace_back(hashToIndex.at(hash), startPos, endPos);
		}
		if (hashes.size() == 0) continue;
		count = hashes.size();
		indexFile.write((char*)&read, 4);
		indexFile.write((char*)&count, 4);
		numIndexedHashes += count;
		for (auto t : hashes)
		{
			indexFile.write((char*)&std::get<0>(t), 4);
			indexFile.write((char*)&std::get<1>(t), 4);
			indexFile.write((char*)&std::get<2>(t), 4);
		}
	}
	redoHashes.close();
	std::remove("index.tmp");
	std::cerr << numIndexedHashes << " indexed positions" << std::endl;
}
