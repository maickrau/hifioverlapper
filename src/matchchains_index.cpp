#include <fstream>
#include <filesystem>
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
	size_t numHashPasses = std::stoull(argv[5]);
	std::string indexPrefix = argv[6];
	std::vector<std::string> readFiles;
	for (size_t i = 7; i < argc; i++)
	{
		readFiles.emplace_back(argv[i]);
	}
	size_t numReads = 0;
	size_t numHashes = 0;
	size_t numIndexedHashes = 0;
	std::ofstream tmpPositionsFile { indexPrefix + ".tmp1", std::ofstream::binary };
	std::ofstream tmpHashesFile { indexPrefix + ".tmp2", std::ofstream::binary };
	MatchIndex matchIndex { k, numWindows, windowSize };
	ReadStorage storage;
	for (auto file : readFiles)
	{
		std::mutex indexMutex;
		storage.iterateReadsFromFile(file, numThreads, false, [&matchIndex, &indexMutex, &tmpPositionsFile, &tmpHashesFile, &numReads, &numHashes, &numIndexedHashes](size_t readName, const std::string& sequence)
		{
			std::vector<std::tuple<uint64_t, uint32_t, uint32_t>> hashes;
			phmap::flat_hash_set<uint64_t> hashesHere;
			matchIndex.iterateWindowChunksFromRead(sequence, [&tmpPositionsFile, &tmpHashesFile, &numReads, &numHashes, &numIndexedHashes, &hashesHere, &hashes](uint32_t startPos, uint32_t endPos, uint64_t hash)
			{
				hashesHere.emplace(hash);
				hashes.emplace_back(hash, startPos, endPos);
			});
			if (hashes.size() == 0) return;
			std::lock_guard<std::mutex> lock { indexMutex };
			numReads += 1;
			numHashes += hashes.size();
			uint32_t read = readName;
			uint32_t count = hashes.size();
			tmpPositionsFile.write((char*)&read, 4);
			tmpPositionsFile.write((char*)&count, 4);
			for (auto t : hashes)
			{
				tmpHashesFile.write((char*)&std::get<0>(t), 8);
				tmpPositionsFile.write((char*)&std::get<1>(t), 4);
				tmpPositionsFile.write((char*)&std::get<2>(t), 4);
			}
		});
	}
	const std::vector<size_t>& readLengths = storage.getRawReadLengths();
	const std::vector<std::string>& readNames = storage.getNames();
	tmpPositionsFile.close();
	tmpHashesFile.close();
	std::cerr << numReads << " reads" << std::endl;
	std::cerr << numHashes << " total positions" << std::endl;
	phmap::flat_hash_map<uint64_t, uint32_t> hashToIndex;
	size_t totalDistinctHashes = 0;
	for (size_t i = 0; i < numHashPasses; i++)
	{
		std::ifstream hashPass { indexPrefix + ".tmp2", std::ios::binary };
		phmap::flat_hash_set<uint64_t> seenOnce;
		while (hashPass.good())
		{
			uint64_t hash;
			hashPass.read((char*)&hash, 8);
			if (!hashPass.good()) break;
			if (hash % numHashPasses != i) continue;
			if (seenOnce.count(hash) == 1)
			{
				if (hashToIndex.count(hash) == 0)
				{
					size_t index = hashToIndex.size();
					hashToIndex[hash] = index;
				}
			}
			else
			{
				seenOnce.insert(hash);
			}
		}
		totalDistinctHashes += seenOnce.size();
	}
	std::cerr << totalDistinctHashes << " distinct hashes" << std::endl;
	std::cerr << hashToIndex.size() << " indexed hashes" << std::endl;
	{
		std::ofstream indexMetadata { indexPrefix + ".metadata", std::ofstream::binary };
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
	std::ofstream indexFile { indexPrefix + ".positions", std::ofstream::binary };
	std::ifstream redoPositions { indexPrefix + ".tmp1", std::ifstream::binary };
	std::ifstream redoHashes { indexPrefix + ".tmp2", std::ifstream::binary };
	while (redoPositions.good())
	{
		uint32_t read;
		uint32_t count;
		redoPositions.read((char*)&read, 4);
		redoPositions.read((char*)&count, 4);
		if (!redoPositions.good()) break;
		std::vector<std::tuple<uint32_t, uint32_t, uint32_t>> hashes;
		for (size_t i = 0; i < count; i++)
		{
			uint64_t hash;
			uint32_t startPos, endPos;
			assert(redoHashes.good());
			redoHashes.read((char*)&hash, 8);
			redoPositions.read((char*)&startPos, 4);
			redoPositions.read((char*)&endPos, 4);
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
	redoPositions.close();
	redoHashes.close();
	std::filesystem::remove(indexPrefix + ".tmp1");
	std::filesystem::remove(indexPrefix + ".tmp2");
	std::cerr << numIndexedHashes << " indexed positions" << std::endl;
}
