#include <mutex>
#include <iostream>
#include <vector>
#include <string>
#include <tuple>
#include <set>
#include "ReadStorage.h"
#include "MatchIndex.h"
#include "UnitigResolver.h"
#include "UnitigGraph.h"
#include "UnitigHelper.h"
#include "HPCConsensus.h"

class ConcatenatedStringStorage
{
public:
	void resize(size_t numItems, size_t k)
	{
		this->k = k;
		sequence.resize(numItems * k);
	}
	std::string getSequence(size_t index) const
	{
		return sequence.substr(index*k, k);
	}
	void setSequence(size_t index, std::string seq)
	{
		assert(index*k+k <= sequence.size());
		for (size_t i = 0; i < k; i++)
		{
			sequence[index*k+i] = seq[i];
		}
	}
	size_t size() const
	{
		return sequence.size() / k;
	}
private:
	std::string sequence;
	size_t k;
};

// defined in MBG.cpp but not exposed, so declare them here too
void loadReadsAsHashesMultithread(HashList& result, const size_t kmerSize, const ReadpartIterator& partIterator, const size_t numThreads, std::ostream& log);
UnitigGraph getUnitigGraph(HashList& hashlist, const size_t minCoverage, const double minUnitigCoverage, const bool keepGaps, const bool oneCovHeuristic);
std::vector<ReadPath> getReadPaths(const UnitigGraph& graph, const HashList& hashlist, const size_t numThreads, const ReadpartIterator& partIterator, const size_t kmerSize);
SparseEdgeContainer getCoveredEdges(const HashList& hashlist, size_t minCoverage);

void loadKmerSequences(const ReadpartIterator& partIterator, const size_t kmerSize, const HashList& reads, const RankBitvector& hasSequence, ConcatenatedStringStorage& kmerSequences)
{
	std::vector<bool> loaded;
	loaded.resize(kmerSequences.size(), false);
	partIterator.iterateHashes([&hasSequence, &kmerSequences, &reads, &loaded, kmerSize](const ReadInfo& read, const SequenceCharType& seq, const SequenceLengthType& poses, const std::string& rawSeq, const std::vector<size_t>& positions, const std::vector<HashType>& hashes)
	{
		for (size_t i = 0; i < positions.size(); i++)
		{
			const HashType fwHash = hashes[i];
			std::pair<size_t, bool> current = reads.getNodeOrNull(fwHash);
			if (!hasSequence.get(current.first)) continue;
			size_t sequenceIndex = hasSequence.getRank(current.first);
			if (loaded[sequenceIndex]) continue;
			std::string sequence = rawSeq.substr(positions[i], kmerSize);
			assert(sequence.size() == kmerSize);
			if (!current.second) sequence = revCompRaw(sequence);
			kmerSequences.setSequence(sequenceIndex, sequence);
			loaded[sequenceIndex] = true;
		}
	});
}

std::vector<std::pair<size_t, bool>> getUniqueAltPath(const HashList& reads, const SparseEdgeContainer& edges, const std::unordered_set<size_t>& pathHashes, const std::pair<size_t, bool> start, const std::pair<size_t, bool> end)
{
	if (end == start) return std::vector<std::pair<size_t, bool>>{};
	std::set<std::pair<size_t, bool>> bwVisited;
	std::vector<std::pair<size_t, bool>> checkStack;
	checkStack.emplace_back(reverse(end));
	while (checkStack.size() > 0)
	{
		auto top = checkStack.back();
		checkStack.pop_back();
		if (bwVisited.count(top) == 1) continue;
		bwVisited.insert(top);
		if (bwVisited.size() > 1000) return std::vector<std::pair<size_t, bool>>{};
		for (auto edge : edges[top])
		{
			if (edge == start) return std::vector<std::pair<size_t, bool>>{};
			if (pathHashes.count(edge.first) == 1) continue;
			checkStack.push_back(edge);
		}
	}
	std::vector<std::pair<size_t, bool>> result;
	result.push_back(start);
	std::set<std::pair<size_t, bool>> visited;
	while (result.back() != end)
	{
		std::pair<size_t, bool> uniqueOption { std::numeric_limits<size_t>::max(), false };
		for (auto edge : edges[result.back()])
		{
			if (bwVisited.count(reverse(edge)) == 0) continue;
			if (uniqueOption.first == std::numeric_limits<size_t>::max())
			{
				uniqueOption = edge;
			}
			else
			{
				return std::vector<std::pair<size_t, bool>>{};
			}
		}
		if (uniqueOption.first == std::numeric_limits<size_t>::max()) return std::vector<std::pair<size_t, bool>>{};
		assert(visited.count(uniqueOption) == 0);
		result.push_back(uniqueOption);
		visited.insert(uniqueOption);
	}
	assert(result.size() > 1);
	assert(result[0] == start);
	assert(result.back() == end);
	return result;
}

std::string getPathSequence(const HashList& reads, const size_t kmerSize, const RankBitvector& hasSequence, const ConcatenatedStringStorage& kmerSequences, const std::vector<std::pair<size_t, bool>>& path)
{
	std::string result;
	for (size_t i = 0; i < path.size(); i++)
	{
		assert(hasSequence.get(path[i].first));
		std::string add = kmerSequences.getSequence(hasSequence.getRank(path[i].first));
		if (!path[i].second) add = revCompRaw(add);
		if (i == 0)
		{
			assert(result.size() == 0);
			result = add;
		}
		else
		{
			size_t overlap = reads.getOverlap(path[i-1], path[i]);
			result.insert(result.end(), add.begin()+overlap, add.end());
		}
	}
	return result;
}

std::string getFixedSequence(const std::string& rawSeq, const HashList& reads, const RankBitvector& hasSequence, const ConcatenatedStringStorage& kmerSequences, const size_t kmerSize, const std::vector<std::tuple<size_t, size_t, std::vector<std::pair<size_t, bool>>>>& fixlist)
{
	assert(fixlist.size() > 0);
	std::string result = rawSeq;
	for (size_t i = fixlist.size()-1; i < fixlist.size(); i--)
	{
		std::string partSequence = getPathSequence(reads, kmerSize, hasSequence, kmerSequences, std::get<2>(fixlist[i]));
		size_t replaceStart = std::get<0>(fixlist[i]);
		size_t replaceEnd = std::get<1>(fixlist[i]);
		assert(replaceEnd > replaceStart);
		assert(partSequence != result.substr(replaceStart, replaceEnd-replaceStart));
		result = result.substr(0, replaceStart) + partSequence + result.substr(replaceEnd);
	}
	assert(result != rawSeq);
	return result;
}

template <typename T>
bool vecmatch(const std::vector<T>& left, const std::vector<T>& right, size_t rightstart, size_t rightend)
{
	if (left.size() != rightend - rightstart) return false;
	for (size_t i = 0; i < left.size(); i++)
	{
		if (left[i] != right[rightstart+i]) return false;
	}
	return true;
}

template <typename F>
void iterateCorrectedSequences(const ReadpartIterator& partIterator, const size_t kmerSize, const HashList& reads, const RankBitvector& hasSequence, const ConcatenatedStringStorage& kmerSequences, const size_t minCoverage, const size_t ambiguousCoverage, F callback)
{
	auto edges = getCoveredEdges(reads, ambiguousCoverage);
	partIterator.iterateHashes([&hasSequence, &kmerSequences, &reads, &edges, kmerSize, minCoverage, ambiguousCoverage, callback](const ReadInfo& read, const SequenceCharType& seq, const SequenceLengthType& poses, const std::string& rawSeq, const std::vector<size_t>& positions, const std::vector<HashType>& hashes)
	{
		if (positions.size() == 0)
		{
			callback(read.readName.first, rawSeq, false);
			return;
		}
		std::vector<std::pair<size_t, bool>> rawPath;
		std::unordered_set<size_t> pathHashes;
		for (size_t i = 0; i < positions.size(); i++)
		{
			const HashType fwHash = hashes[i];
			rawPath.push_back(reads.getNodeOrNull(fwHash));
			assert(rawPath.back().first < reads.size());
			pathHashes.insert(rawPath.back().first);
		}
		assert(rawPath.size() == positions.size());
		std::vector<std::tuple<size_t, size_t, std::vector<std::pair<size_t, bool>>>> fixlist;
		size_t lastSolid = std::numeric_limits<size_t>::max();
		for (size_t i = 0; i < positions.size(); i++)
		{
			if (reads.coverage.get(rawPath[i].first) < ambiguousCoverage) continue;
			if (lastSolid == std::numeric_limits<size_t>::max())
			{
				assert(fixlist.size() == 0);
				lastSolid = i;
				continue;
			}
			if (i == lastSolid+1)
			{
				lastSolid = i;
				continue;
			}
			std::vector<std::pair<size_t, bool>> uniqueAltPath = getUniqueAltPath(reads, edges, pathHashes, rawPath[lastSolid], rawPath[i]);
			if (uniqueAltPath.size() == 0)
			{
				lastSolid = i;
				continue;
			}
			assert(uniqueAltPath[0] == rawPath[lastSolid]);
			assert(uniqueAltPath.back() == rawPath[i]);
			assert(!vecmatch(uniqueAltPath, rawPath, lastSolid, i+1));
			bool validAltPath = true;
			for (size_t j = 0; j < uniqueAltPath.size(); j++)
			{
				if (reads.coverage.get(uniqueAltPath[j].first) < minCoverage) validAltPath = false;
			}
			if (validAltPath) fixlist.emplace_back(positions[lastSolid], positions[i]+kmerSize, uniqueAltPath);
			lastSolid = i;
		}
		if (fixlist.size() == 0)
		{
			callback(read.readName.first, rawSeq, false);
			return;
		}
		std::string fixedSequence = getFixedSequence(rawSeq, reads, hasSequence, kmerSequences, kmerSize, fixlist);
		callback(read.readName.first, fixedSequence, true);
	});
}

int main(int argc, char** argv)
{
	size_t numThreads = std::stoi(argv[1]);
	size_t kmerSize = std::stoull(argv[2]);
	std::vector<std::string> readFiles;
	for (size_t i = 3; i < argc; i++)
	{
		readFiles.emplace_back(argv[i]);
	}
	size_t windowSize = kmerSize-30;
	size_t solidCoverage = 5;
	size_t ambiguousCoverage = 2;
	ReadpartIterator partIterator { kmerSize, windowSize, ErrorMasking::No, numThreads, readFiles, false, "" };
	HashList reads { kmerSize };
	loadReadsAsHashesMultithread(reads, kmerSize, partIterator, numThreads, std::cerr);
	RankBitvector hasSequence { reads.size() };
	size_t totalWithSequence = 0;
	for (size_t i = 0; i < reads.size(); i++)
	{
		if (reads.coverage.get(i) < solidCoverage) continue;
		totalWithSequence += 1;
		hasSequence.set(i, true);
	}
	std::cerr << totalWithSequence << " solid k-mers" << std::endl;
	hasSequence.buildRanks();
	ConcatenatedStringStorage kmerSequences;
	kmerSequences.resize(totalWithSequence, kmerSize);
	std::cerr << "loading k-mer sequences" << std::endl;
	loadKmerSequences(partIterator, kmerSize, reads, hasSequence, kmerSequences);
	std::cerr << "correcting reads" << std::endl;
	std::mutex writeMutex;
	size_t countCorrected = 0;
	size_t countNotCorrected = 0;
	iterateCorrectedSequences(partIterator, kmerSize, reads, hasSequence, kmerSequences, solidCoverage, ambiguousCoverage, [&writeMutex, &countCorrected, &countNotCorrected](const std::string& readname, const std::string& readseq, const bool corrected)
	{
		std::lock_guard<std::mutex> lock { writeMutex };
		if (corrected)
		{
			countCorrected += 1;
		}
		else
		{
			countNotCorrected += 1;
		}
		std::cout << ">" << readname << std::endl;
		std::cout << readseq << std::endl;
	});
	std::cerr << countCorrected << " reads corrected" << std::endl;
	std::cerr << countNotCorrected << " reads not corrected" << std::endl;
}
