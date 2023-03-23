#include <queue>
#include <unordered_set>
#include <vector>
#include <set>
#include "KmerCorrector.h"

thread_local std::vector<uint8_t> localCoverage;
thread_local std::vector<bool> hasLocalFw;
thread_local std::vector<bool> hasLocalBw;

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

class DistantKmerComparer
{
public:
	bool operator()(const std::pair<size_t, std::pair<size_t, bool>> left, const std::pair<size_t, std::pair<size_t, bool>> right) const
	{
		return left.first > right.first;
	}
};

std::vector<std::pair<size_t, bool>> getUniqueAltPath(const HashList& reads, const size_t kmerSize, const SparseEdgeContainer& edges, const RankBitvector& hasSequence, const std::pair<size_t, bool> start, const std::pair<size_t, bool> end, const size_t maxLength)
{
	if (end == start) return std::vector<std::pair<size_t, bool>>{};
	std::set<std::pair<size_t, bool>> bwVisited;
	std::priority_queue<std::pair<size_t, std::pair<size_t, bool>>, std::vector<std::pair<size_t, std::pair<size_t, bool>>>, DistantKmerComparer> checkStack;
	checkStack.emplace(0, reverse(end));
	while (checkStack.size() > 0)
	{
		auto topandlen = checkStack.top();
		checkStack.pop();
		size_t distance = topandlen.first;
		std::pair<size_t, bool> top = topandlen.second;
		if (distance > maxLength) continue;
		if (bwVisited.count(top) == 1) continue;
		bwVisited.insert(top);
		if (bwVisited.size() > 1000) return std::vector<std::pair<size_t, bool>>{};
		for (auto edge : edges[top])
		{
			if (edge == start) return std::vector<std::pair<size_t, bool>>{};
			if (!hasSequence.get(edge.first)) continue;
			size_t extra = kmerSize - reads.getOverlap(top, edge);
			checkStack.emplace(distance + extra, edge);
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

void KmerCorrector::initializeThreadCoverageFullGraph()
{
	localCoverage.resize(0, 0);
}

void KmerCorrector::filterToReadCoverage(const ReadStorage& iterator)
{
	localCoverage.resize(removedHomopolymerError.size(), 0);
	hasLocalBw.resize(removedHomopolymerError.size(), false);
	hasLocalFw.resize(removedHomopolymerError.size(), false);
	for (size_t i = 0; i < localCoverage.size(); i++)
	{
		localCoverage[i] = 0;
		hasLocalFw[i] = false;
		hasLocalBw[i] = false;
	}
	iterator.iterateReadsAndHashesFromStorage([this](const size_t readId, const std::string& rawSeq, const std::vector<size_t>& positions, const std::vector<HashType>& hashes)
	{
		for (auto hash : hashes)
		{
			std::pair<size_t, bool> current = reads.getNodeOrNull(hash);
			if (current.first >= hasSequence.size()) continue;
			if (!hasSequence.get(current.first)) continue;
			size_t rank = hasSequence.getRank(current.first);
			if (localCoverage[rank] < std::numeric_limits<uint8_t>::max()) localCoverage[rank] += 1;
			if (current.second) hasLocalFw[rank] = true;
			if (!current.second) hasLocalBw[rank] = true;
		}
	});
}

std::pair<size_t, bool> KmerCorrector::findBubble(std::pair<size_t, bool> start) const
{
	std::vector<std::pair<size_t, bool>> S;
	S.push_back(start);
	std::unordered_set<std::pair<size_t, bool>> visited;
	std::unordered_set<std::pair<size_t, bool>> seen;
	seen.insert(start);
	std::pair<size_t, bool> bubbleEnd { std::numeric_limits<size_t>::max(), false };
	while (S.size() > 0)
	{
		const std::pair<size_t, bool> v = S.back();
		S.pop_back();
		assert(seen.count(v) == 1);
		seen.erase(v);
		if (visited.count(reverse(v)) == 1) return bubbleEnd;
		assert(visited.count(v) == 0);
		visited.insert(v);
		if (edges[v].size() == 0) return bubbleEnd;
		bool hasEdge = false;
		for (auto u : edges[v])
		{
			if (!hasSequence.get(u.first)) continue;
			if (u == v) return bubbleEnd;
			if (u == start) return bubbleEnd;
			assert(visited.count(u) == 0);
			seen.insert(u);
			hasEdge = true;
			bool hasUnvisitedInneighbor = false;
			for (auto w : edges[reverse(u)])
			{
				if (w == u) return bubbleEnd;
				if (!hasSequence.get(w.first)) continue;
				if (visited.count(reverse(w)) == 0)
				{
					hasUnvisitedInneighbor = true;
				}
			}
			if (!hasUnvisitedInneighbor) S.push_back(u);
		}
		if (!hasEdge) return bubbleEnd;
		if (S.size() == 1 && seen.size() == 1 && seen.count(S[0]) == 1)
		{
			bubbleEnd = S[0];
			break;
		}
	}
	return bubbleEnd;
}

void KmerCorrector::forbidPathNodes(const std::vector<std::pair<size_t, bool>>& path)
{
	for (auto node : path)
	{
		assert(hasSequence.get(node.first));
		size_t rank = hasSequence.getRank(node.first);
		removedHomopolymerError[rank] = true;
	}
}

void KmerCorrector::allowPathNodes(const std::vector<std::pair<size_t, bool>>& path)
{
	for (auto node : path)
	{
		assert(hasSequence.get(node.first));
		size_t rank = hasSequence.getRank(node.first);
		removedHomopolymerError[rank] = false;
	}
}

bool homopolymerAdjacent(const std::vector<size_t>& left, const std::vector<size_t>& right)
{
	assert(left.size() == right.size());
	bool foundDifference = false;
	for (size_t i = 0; i < left.size(); i++)
	{
		if (left[i] == right[i]) continue;
		if (left[i] == right[i]+1)
		{
			if (foundDifference) return false;
			foundDifference = true;
		}
		if (left[i]+1 == right[i])
		{
			if (foundDifference) return false;
			foundDifference = true;
		}
		return false;
	}
	return true;
}

size_t KmerCorrector::getCoverage(size_t kmer) const
{
	if (localCoverage.size() == 0) return reads.coverage.get(kmer);
	if (!hasSequence.get(kmer)) return 0;
	size_t rank = hasSequence.getRank(kmer);
	if (!hasLocalBw[rank] || !hasLocalFw[rank]) return 0;
	return localCoverage[rank];
}

size_t KmerCorrector::getPathCoverage(const std::vector<std::pair<size_t, bool>>& path) const
{
	assert(path.size() >= 1);
	size_t result = getCoverage(path[0].first);
	for (size_t i = 1; i < path.size(); i++)
	{
		result = std::min(result, getCoverage(path[i].first));
	}
	return result;
}

std::pair<std::string, std::vector<size_t>> KmerCorrector::getHomopolymerCompressedPathSequence(const std::vector<std::pair<size_t, bool>>& path) const
{
	std::string resultExpanded = getPathSequence(reads, kmerSize, hasSequence, kmerSequences, path);
	std::string resultCompressed;
	std::vector<size_t> resultLengths;
	assert(resultExpanded.size() >= 1);
	resultCompressed.push_back(resultExpanded[0]);
	resultLengths.push_back(1);
	for (size_t i = 1; i < resultExpanded.size(); i++)
	{
		if (resultExpanded[i] == resultExpanded[i-1])
		{
			resultLengths.back() += 1;
		}
		else
		{
			resultCompressed.push_back(resultExpanded[i]);
			resultLengths.push_back(1);
		}
	}
	assert(resultCompressed.size() == resultLengths.size());
	return std::make_pair(resultCompressed, resultLengths);
}

void KmerCorrector::enumeratePathsRecursion(std::vector<std::vector<std::pair<size_t, bool>>>& result, std::vector<std::pair<size_t, bool>>& currentPath, const std::pair<size_t, bool> end, const size_t maxCount) const
{
	if (currentPath.size() > 2000)
	{
		// fill to make sure nothing is returned
		for (size_t i = 0; i < maxCount; i++)
		{
			result.emplace_back();
		}
		return;
	}
	assert(currentPath.size() < edges.size()+5);
	assert(currentPath.size() >= 1);
	assert(currentPath.back().first < edges.size());
	for (auto node : edges[currentPath.back()])
	{
		if (!hasSequence.get(node.first)) continue;
		currentPath.push_back(node);
		if (node == end)
		{
			result.push_back(currentPath);
			currentPath.pop_back();
			if (result.size() >= maxCount) return;
		}
		else
		{
			enumeratePathsRecursion(result, currentPath, end, maxCount);
			currentPath.pop_back();
			if (result.size() >= maxCount) return;
		}
	}
}

std::vector<std::vector<std::pair<size_t, bool>>> KmerCorrector::enumeratePaths(const std::pair<size_t, bool> start, const std::pair<size_t, bool> end, const size_t maxCount) const
{
	std::vector<std::vector<std::pair<size_t, bool>>> result;
	std::vector<std::pair<size_t, bool>> currentPath;
	currentPath.push_back(start);
	enumeratePathsRecursion(result, currentPath, end, maxCount);
	if (result.size() >= maxCount)
	{
		result.clear();
	}
	return result;
}

void KmerCorrector::forbidHomopolymerAlleles(const std::pair<size_t, bool> start, const std::pair<size_t, bool> end)
{
	assert(start.first < edges.size());
	assert(end.first < edges.size());
	assert(end != start);
	std::vector<std::vector<std::pair<size_t, bool>>> paths = enumeratePaths(start, end, 20);
	if (paths.size() < 2) return;
	std::string hpcSequence;
	std::vector<std::vector<size_t>> pathLengths;
	for (size_t i = 0; i < paths.size(); i++)
	{
		auto seqHere = getHomopolymerCompressedPathSequence(paths[i]);
		if (i == 0)
		{
			hpcSequence = seqHere.first;
		}
		else
		{
			if (seqHere.first != hpcSequence) return;
		}
		pathLengths.emplace_back(seqHere.second);
	}
	std::vector<size_t> parent;
	parent.resize(paths.size());
	std::vector<size_t> pathCoverages;
	pathCoverages.resize(paths.size());
	for (size_t i = 0; i < paths.size(); i++)
	{
		parent[i] = i;
		pathCoverages[i] = getPathCoverage(paths[i]);
	}
	for (size_t i = 0; i < paths.size(); i++)
	{
		for (size_t j = i+1; j < paths.size(); j++)
		{
			if (homopolymerAdjacent(pathLengths[i], pathLengths[j]))
			{
				while (parent[i] != parent[parent[i]]) parent[i] = parent[parent[i]];
				while (parent[j] != parent[parent[j]]) parent[j] = parent[parent[j]];
				size_t leftParent = parent[i];
				size_t rightParent = parent[j];
				if (pathCoverages[leftParent] > pathCoverages[rightParent])
				{
					parent[rightParent] = leftParent;
				}
				else
				{
					parent[leftParent] = rightParent;
				}
			}
		}
	}
	std::vector<std::vector<size_t>> clusters;
	clusters.resize(paths.size());
	for (size_t i = 0; i < paths.size(); i++)
	{
		while (parent[i] != parent[parent[i]]) parent[i] = parent[parent[i]];
		assert(parent[i] == i || pathCoverages[parent[i]] >= pathCoverages[i]);
		clusters[parent[i]].push_back(i);
	}
	std::vector<size_t> keepMaxes;
	for (size_t i = 0; i < clusters.size(); i++)
	{
		size_t clusterMax = 0;
		for (auto path : clusters[i])
		{
			if (pathCoverages[path] > pathCoverages[clusterMax]) clusterMax = path;
		}
		bool coverageValid = true;
		for (auto path : clusters[i])
		{
			if (path != clusterMax && pathCoverages[path]*2 > pathCoverages[clusterMax]) coverageValid = false;
		}
		if (!coverageValid)
		{
			for (auto path : clusters[i]) keepMaxes.push_back(path);
			continue;
		}
		for (auto path : clusters[i])
		{
			forbidPathNodes(paths[path]);
		}
		keepMaxes.push_back(clusterMax);
	}
	for (auto path : keepMaxes)
	{
		allowPathNodes(paths[path]);
	}
}

void KmerCorrector::forbidHomopolymerErrors()
{
	for (size_t i = 0; i < hasSequence.size(); i++)
	{
		if (!hasSequence.get(i)) continue;
		auto bubbleEndFw = findBubble(std::make_pair(i, true));
		if (bubbleEndFw.first != std::numeric_limits<size_t>::max())
		{
			forbidHomopolymerAlleles(std::make_pair(i, true), bubbleEndFw);
		}
		auto bubbleEndBw = findBubble(std::make_pair(i, false));
		if (bubbleEndBw.first != std::numeric_limits<size_t>::max())
		{
			forbidHomopolymerAlleles(std::make_pair(i, false), bubbleEndBw);
		}
	}
}

void KmerCorrector::buildGraph(const ReadpartIterator& partIterator)
{
	loadReadsAsHashesAndKmerSequencesMultithread(reads, kmerSize, [&partIterator](auto callback)
	{
		partIterator.iterateHashes([callback](const ReadInfo& read, const SequenceCharType& seq, const SequenceLengthType& poses, const std::string& rawSeq, const std::vector<size_t>& positions, const std::vector<HashType>& hashes)
		{
			callback(rawSeq, positions, hashes);
		});
	});
	edges = getCoveredEdges(reads, minAmbiguousCoverage);
	removedHomopolymerError.resize(kmerSequences.size(), false);
	forbidHomopolymerErrors();
}

void KmerCorrector::buildGraph(const ReadStorage& storage)
{
	loadReadsAsHashesAndKmerSequencesMultithread(reads, kmerSize, [&storage](auto callback)
	{
		storage.iterateReadsAndHashesFromStorage([callback](const size_t readId, const std::string& rawSeq, const std::vector<size_t>& positions, const std::vector<HashType>& hashes)
		{
			callback(rawSeq, positions, hashes);
		});
	});
	edges = getCoveredEdges(reads, minAmbiguousCoverage);
	removedHomopolymerError.resize(kmerSequences.size(), false);
	forbidHomopolymerErrors();
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
	for (size_t i = 0; i < positions.size(); i++)
	{
		const HashType fwHash = hashes[i];
		rawPath.push_back(reads.getNodeOrNull(fwHash));
		assert(rawPath.back().first < reads.size());
	}
	assert(rawPath.size() == positions.size());
	std::vector<std::tuple<size_t, size_t, std::vector<std::pair<size_t, bool>>>> fixlist;
	size_t lastSolid = std::numeric_limits<size_t>::max();
	for (size_t i = 0; i < positions.size(); i++)
	{
		if (getCoverage(rawPath[i].first) < minSolidCoverage) continue;
		if (!hasSequence.get(rawPath[i].first)) continue;
		size_t rank = hasSequence.getRank(rawPath[i].first);
		if (removedHomopolymerError[rank]) continue;
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
		size_t maxLength = positions[i] - positions[lastSolid] + 100;
		std::vector<std::pair<size_t, bool>> uniqueAltPath = getUniqueAltPath(reads, kmerSize, edges, hasSequence, rawPath[lastSolid], rawPath[i], maxLength);
		if (uniqueAltPath.size() == 0 || vecmatch(uniqueAltPath, rawPath, lastSolid, i+1))
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
			if (getCoverage(uniqueAltPath[j].first) < minSolidCoverage) validAltPath = false;
			if (!hasSequence.get(uniqueAltPath[j].first)) validAltPath = false;
		}
		if (validAltPath) fixlist.emplace_back(positions[lastSolid], positions[i]+kmerSize, uniqueAltPath);
		lastSolid = i;
	}
	if (fixlist.size() == 0) return std::make_pair(rawSeq, false);
	std::string fixedSequence = getFixedSequence(rawSeq, reads, hasSequence, kmerSequences, kmerSize, fixlist);
	return std::make_pair(fixedSequence, true);
}
