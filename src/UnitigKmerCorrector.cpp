#include <cassert>
#include <mutex>
#include "UnitigKmerCorrector.h"

UnitigKmerCorrector::UnitigKmerCorrector(size_t k) :
	kmerSize(k),
	unitigs(k)
{
}

void UnitigKmerCorrector::build(const ReadpartIterator& partIterator)
{
	std::mutex mutex;
	std::cerr << "collect k-mers" << std::endl;
	partIterator.iterateHashes([this, &mutex](const ReadInfo& read, const SequenceCharType& seq, const SequenceLengthType& poses, const std::string& rawSeq, const std::vector<size_t>& positions, const std::vector<HashType>& hashes)
	{
		std::lock_guard lock { mutex };
		std::pair<size_t, bool> last { std::numeric_limits<size_t>::max(), true };
		for (size_t i = 0; i < hashes.size(); i++)
		{
			std::pair<size_t, bool> thisKmer = unitigs.getNode(hashes[i]);
			if (last.first != std::numeric_limits<size_t>::max())
			{
				unitigs.addEdge(last, thisKmer, kmerSize - (positions[i] - positions[i-1]));
			}
			last = thisKmer;
		}
	});
	std::cerr << unitigs.numHashes() << " hashes" << std::endl;
	std::cerr << "build unitigs" << std::endl;
	unitigs.buildUnitigGraph();
	std::cerr << unitigs.numUnitigs() << " unitigs" << std::endl;
	std::cerr << "build unitig sequences and collect read paths" << std::endl;
	partIterator.iterateHashes([this, &mutex](const ReadInfo& read, const SequenceCharType& seq, const SequenceLengthType& poses, const std::string& rawSeq, const std::vector<size_t>& positions, const std::vector<HashType>& hashes)
	{
		std::lock_guard lock { mutex };
		for (size_t i = 0; i < hashes.size(); i++)
		{
			std::pair<size_t, bool> thisKmer = unitigs.getNodeOrNull(hashes[i]);
			unitigs.addKmerSequence(thisKmer, rawSeq, positions[i], positions[i]+kmerSize);
		}
		reads.emplace_back();
		reads.back().name = read.readName.first;
		std::tie(reads.back().leftClip, reads.back().rightClip, reads.back().unitigPath) = unitigs.getPath(hashes);
		if (reads.back().unitigPath.size() == 0)
		{
			reads.back().leftHanger = rawSeq;
		}
		else
		{
			reads.back().leftHanger = rawSeq.substr(0, positions[0]);
			reads.back().rightHanger = rawSeq.substr(positions.back() + kmerSize);
		}
	});
	std::cerr << reads.size() << " reads" << std::endl;
	std::cerr << "finalize" << std::endl;
	unitigs.finalizeSequences();
	std::cerr << unitigs.totalBps() << " bp" << std::endl;
}

void UnitigKmerCorrector::getAnchorSpanners(phmap::flat_hash_map<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>, std::map<std::vector<std::pair<size_t, bool>>, size_t>>& anchorSpannerCounts, const phmap::flat_hash_set<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>& adjacentAnchors, const phmap::flat_hash_set<std::pair<size_t, bool>>& isAnchor, size_t read) const
{
	size_t lastAnchor = std::numeric_limits<size_t>::max();
	for (size_t i = 0; i < reads[read].unitigPath.size(); i++)
	{
		if (isAnchor.count(reads[read].unitigPath[i]) == 0) continue;
		if (lastAnchor == std::numeric_limits<size_t>::max())
		{
			lastAnchor = i;
			continue;
		}
		std::pair<size_t, bool> oldAnchor = reads[read].unitigPath[lastAnchor];
		std::pair<size_t, bool> thisAnchor = reads[read].unitigPath[i];
		if (adjacentAnchors.count(std::make_pair(oldAnchor, thisAnchor)) == 1)
		{
			std::vector<std::pair<size_t, bool>> subpath { reads[read].unitigPath.begin()+lastAnchor+1, reads[read].unitigPath.begin()+i };
			anchorSpannerCounts[std::make_pair(oldAnchor, thisAnchor)][subpath] += 1;
		}
		if (adjacentAnchors.count(std::make_pair(reverse(thisAnchor), reverse(oldAnchor))) == 1)
		{
			std::vector<std::pair<size_t, bool>> subpath { reads[read].unitigPath.begin()+lastAnchor+1, reads[read].unitigPath.begin()+i };
			std::reverse(subpath.begin(), subpath.end());
			for (size_t j = 0; j < subpath.size(); j++)
			{
				subpath[j] = reverse(subpath[j]);
			}
			anchorSpannerCounts[std::make_pair(reverse(thisAnchor), reverse(oldAnchor))][subpath] += 1;
		}
	}
}

std::string UnitigKmerCorrector::getCorrectedSequence(size_t readIndex, const std::vector<size_t>& context, size_t minAmbiguousCoverage, size_t minSafeCoverage) const
{
	if (reads[readIndex].unitigPath.size() == 1) return getRawSequence(readIndex);
	phmap::flat_hash_map<size_t, size_t> keyReadUnitigCoverage;
	phmap::flat_hash_map<size_t, size_t> unitigCoverage;
	for (size_t read : context)
	{
		for (const auto& node : reads[read].unitigPath)
		{
			unitigCoverage[node.first] += 1;
		}
	}
	for (const auto& node : reads[readIndex].unitigPath)
	{
		keyReadUnitigCoverage[node.first] += 1;
	}
	std::vector<size_t> anchorIndices;
	phmap::flat_hash_set<std::pair<size_t, bool>> isAnchor;
	for (size_t i = 0; i < reads[readIndex].unitigPath.size(); i++)
	{
		if (unitigCoverage[reads[readIndex].unitigPath[i].first] < minSafeCoverage) continue;
		anchorIndices.push_back(i);
		isAnchor.insert(reads[readIndex].unitigPath[i]);
	}
	if (anchorIndices.size() < 2) return getRawSequence(readIndex);
	phmap::flat_hash_set<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>> adjacentAnchors;
	for (size_t i = 1; i < anchorIndices.size(); i++)
	{
		adjacentAnchors.emplace(reads[readIndex].unitigPath[anchorIndices[i-1]], reads[readIndex].unitigPath[anchorIndices[i]]);
	}
	phmap::flat_hash_map<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>, std::map<std::vector<std::pair<size_t, bool>>, size_t>> anchorSpannerCounts;
	for (size_t read : context)
	{
		getAnchorSpanners(anchorSpannerCounts, adjacentAnchors, isAnchor, read);
	}
	std::vector<std::pair<size_t, bool>> correctedLocalPath;
	correctedLocalPath.insert(correctedLocalPath.end(), reads[readIndex].unitigPath.begin(), reads[readIndex].unitigPath.begin() + anchorIndices[0]+1);
	for (size_t i = 1; i < anchorIndices.size(); i++)
	{
		size_t safeCount = 0;
		size_t ambiguousCount = 0;
		for (const auto& pair : anchorSpannerCounts[std::make_pair(reads[readIndex].unitigPath[anchorIndices[i-1]], reads[readIndex].unitigPath[anchorIndices[i]])])
		{
			if (pair.second >= minAmbiguousCoverage) ambiguousCount += 1;
			if (pair.second >= minSafeCoverage) safeCount += 1;
		}
		if (ambiguousCount != 1 || safeCount != 1)
		{
			correctedLocalPath.insert(correctedLocalPath.end(), reads[readIndex].unitigPath.begin()+anchorIndices[i-1]+1, reads[readIndex].unitigPath.begin()+anchorIndices[i]+1);
			continue;
		}
		assert(safeCount == 1);
		assert(ambiguousCount == 1);
		for (const auto& pair : anchorSpannerCounts[std::make_pair(reads[readIndex].unitigPath[anchorIndices[i-1]], reads[readIndex].unitigPath[anchorIndices[i]])])
		{
			if (pair.second < minSafeCoverage) continue;
			correctedLocalPath.insert(correctedLocalPath.end(), pair.first.begin(), pair.first.end());
		}
		correctedLocalPath.emplace_back(reads[readIndex].unitigPath[anchorIndices[i]]);
	}
	correctedLocalPath.insert(correctedLocalPath.end(), reads[readIndex].unitigPath.begin()+anchorIndices.back()+1, reads[readIndex].unitigPath.end());
	assert(correctedLocalPath.size() >= 2);
	assert(correctedLocalPath[0] == reads[readIndex].unitigPath[0]);
	assert(correctedLocalPath.back() == reads[readIndex].unitigPath.back());
	std::string result = reads[readIndex].leftHanger + unitigs.getSequence(correctedLocalPath, reads[readIndex].leftClip, reads[readIndex].rightClip) + reads[readIndex].rightHanger;
	return result;
}

std::string UnitigKmerCorrector::getRawSequence(size_t index) const
{
	std::string result;
	result = reads[index].leftHanger + unitigs.getSequence(reads[index].unitigPath, reads[index].leftClip, reads[index].rightClip) + reads[index].rightHanger;
	return result;
}

size_t UnitigKmerCorrector::numReads() const
{
	return reads.size();
}

const std::string& UnitigKmerCorrector::getName(size_t index) const
{
	return reads[index].name;
}
