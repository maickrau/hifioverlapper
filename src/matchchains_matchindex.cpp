#include <fstream>
#include <mutex>
#include <iostream>
#include <vector>
#include <string>
#include "ReadStorage.h"
#include "MatchIndex.h"

int main(int argc, char** argv)
{
	std::string indexPrefix = argv[1];
	size_t indexSplitThis = std::stoull(argv[2]);
	size_t indexSplitTotalCount = std::stoull(argv[3]);
	std::vector<size_t> readLengths;
	std::vector<std::string> readNames;
	size_t countHashes = 0;
	size_t firstIndexedRead;
	size_t firstNonindexedRead;
	{
		std::ifstream file { indexPrefix + ".metadata", std::ios::binary };
		size_t readCount = 0;
		file.read((char*)&countHashes, sizeof(size_t));
		file.read((char*)&readCount, sizeof(size_t));
		firstIndexedRead = (double)indexSplitThis/(double)indexSplitTotalCount * readCount;
		firstNonindexedRead = (double)(indexSplitThis+1)/(double)indexSplitTotalCount * readCount;
		if (indexSplitThis == indexSplitTotalCount-1) firstNonindexedRead = readCount;
		readLengths.resize(readCount);
		readNames.resize(readCount);
		for (size_t i = 0; i < readLengths.size(); i++)
		{
			file.read((char*)&readLengths[i], sizeof(size_t));
			size_t size = 0;
			file.read((char*)&size, sizeof(size_t));
			readNames[i].resize(size);
			file.read((char*)readNames[i].data(), size);
		}
	}
	std::cerr << readLengths.size() << " reads" << std::endl;
	std::cerr << countHashes << " indexed hashes" << std::endl;
	MatchIndex matchIndex { 0, 0, 0 };
	{
		std::vector<size_t> countHashesPerBucket;
		countHashesPerBucket.resize(countHashes, 0);
		std::ifstream file { indexPrefix + ".positions", std::ios::binary };
		while (file.good())
		{
			uint32_t read = 0;
			uint32_t count = 0;
			file.read((char*)&read, 4);
			file.read((char*)&count, 4);
			if (!file.good()) break;
			assert(count >= 1);
			assert(read < readLengths.size());
			std::vector<std::tuple<uint32_t, uint32_t, uint32_t>> hashes;
			for (size_t i = 0; i < count; i++)
			{
				uint32_t hash;
				uint32_t startPos;
				uint32_t endPos;
				file.read((char*)&hash, 4);
				file.read((char*)&startPos, 4);
				file.read((char*)&endPos, 4);
				countHashesPerBucket[hash] += 1;
			}
		}
		matchIndex.initBuckets(countHashes, countHashesPerBucket);
	}
	{
		std::ifstream file { indexPrefix + ".positions", std::ios::binary };
		while (file.good())
		{
			uint32_t read = 0;
			uint32_t count = 0;
			file.read((char*)&read, 4);
			file.read((char*)&count, 4);
			if (!file.good()) break;
			assert(count >= 1);
			assert(read < readLengths.size());
			std::vector<std::tuple<uint32_t, uint32_t, uint32_t>> hashes;
			for (size_t i = 0; i < count; i++)
			{
				uint32_t hash;
				uint32_t startPos;
				uint32_t endPos;
				file.read((char*)&hash, 4);
				file.read((char*)&startPos, 4);
				file.read((char*)&endPos, 4);
				hashes.emplace_back(hash, startPos, endPos);
			}
			if (read >= firstIndexedRead && read < firstNonindexedRead) matchIndex.addHashes(read, hashes);
		}
	}
	size_t countMatches = 0;
	{
		std::mutex printMutex;
		std::ifstream file { indexPrefix + ".positions", std::ios::binary };
		while (file.good())
		{
			uint32_t read = 0;
			uint32_t count = 0;
			file.read((char*)&read, 4);
			file.read((char*)&count, 4);
			if (!file.good()) break;
			assert(count >= 1);
			std::vector<std::tuple<uint32_t, uint32_t, uint32_t>> hashes;
			for (size_t i = 0; i < count; i++)
			{
				uint32_t hash;
				uint32_t startPos;
				uint32_t endPos;
				file.read((char*)&hash, 4);
				file.read((char*)&startPos, 4);
				file.read((char*)&endPos, 4);
				hashes.emplace_back(hash, startPos, endPos);
			}
			assert(file.good());
			matchIndex.iterateMatchNamesOneRead(read, 10000, readNames, readLengths, hashes, [&printMutex, &countMatches](const std::string& left, const size_t leftlen, const size_t leftstart, const size_t leftend, const bool leftFw, const std::string& right, const size_t rightlen, const size_t rightstart, const size_t rightend, const bool rightFw)
			{
				std::lock_guard<std::mutex> lock { printMutex };
				std::cout << left << "\t" << leftlen << "\t" << leftstart << "\t" << leftend << "\t" << (leftFw ? "fw" : "bw") << "\t" << right << "\t" << rightlen << "\t" << rightstart << "\t" << rightend << "\t" << (rightFw ? "fw" : "bw") << std::endl;
				countMatches += 1;
			});
		}
	}
	std::cerr << countMatches << " chain matches" << std::endl;
}
