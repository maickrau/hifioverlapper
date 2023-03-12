#include <unordered_set>
#include <vector>
#include <set>
#include "KmerCorrector.h"

// defined in MBG.cpp but not exposed, so declare them here too
void loadReadsAsHashesMultithread(HashList& result, const size_t kmerSize, const ReadpartIterator& partIterator, const size_t numThreads, std::ostream& log);
SparseEdgeContainer getCoveredEdges(const HashList& hashlist, size_t minCoverage);

void ConcatenatedStringStorage::resize(size_t numItems, size_t k)
{
	this->k = k;
	sequence.resize(numItems * k);
}

std::string ConcatenatedStringStorage::getSequence(size_t index) const
{
	return sequence.substr(index*k, k);
}

void ConcatenatedStringStorage::setSequence(size_t index, std::string seq)
{
	assert(index*k+k <= sequence.size());
	for (size_t i = 0; i < k; i++)
	{
		sequence[index*k+i] = seq[i];
	}
}

size_t ConcatenatedStringStorage::size() const
{
	return sequence.size() / k;
}

KmerCorrector::KmerCorrector(size_t kmerSize, size_t minSolidCoverage, size_t minAmbiguousCoverage) :
kmerSize(kmerSize),
minSolidCoverage(minSolidCoverage),
minAmbiguousCoverage(minAmbiguousCoverage),
reads(kmerSize)
{

}

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

void KmerCorrector::buildGraph(const ReadpartIterator& partIterator, size_t numThreads)
{
	loadReadsAsHashesMultithread(reads, kmerSize, partIterator, numThreads, std::cerr);
	hasSequence.resize(reads.size());
	size_t totalWithSequence = 0;
	for (size_t i = 0; i < reads.size(); i++)
	{
		if (reads.coverage.get(i) < minSolidCoverage) continue;
		totalWithSequence += 1;
		hasSequence.set(i, true);
	}
	std::cerr << totalWithSequence << " solid k-mers" << std::endl;
	hasSequence.buildRanks();
	kmerSequences.resize(totalWithSequence, kmerSize);
	std::cerr << "loading k-mer sequences" << std::endl;
	loadKmerSequences(partIterator, kmerSize, reads, hasSequence, kmerSequences);
	edges = getCoveredEdges(reads, minAmbiguousCoverage);
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

std::pair<std::string, bool> KmerCorrector::getCorrectedSequence(const std::string& rawSeq, const std::vector<size_t>& positions, const std::vector<HashType>& hashes) const
{
	if (positions.size() == 0) return std::make_pair(rawSeq, false);
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
		if (reads.coverage.get(rawPath[i].first) < minAmbiguousCoverage) continue;
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
			if (reads.coverage.get(uniqueAltPath[j].first) < minSolidCoverage) validAltPath = false;
		}
		if (validAltPath) fixlist.emplace_back(positions[lastSolid], positions[i]+kmerSize, uniqueAltPath);
		lastSolid = i;
	}
	if (fixlist.size() == 0) return std::make_pair(rawSeq, false);
	std::string fixedSequence = getFixedSequence(rawSeq, reads, hasSequence, kmerSequences, kmerSize, fixlist);
	return std::make_pair(fixedSequence, true);
}
