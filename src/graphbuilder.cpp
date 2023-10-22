#include <mutex>
#include <iostream>
#include <vector>
#include <string>
#include "ReadStorage.h"
#include "MatchIndex.h"
#include "SparseEdgeContainer.h"
#include "MostlySparse2DHashmap.h"
#include "CommonUtils.h"
#include "ReadHelper.h"

const size_t minCoverage = 3;

void getUnitig(std::vector<std::pair<size_t, bool>>& belongsToUnitig, const std::vector<bool>& keptNodes, const SparseEdgeContainer& edges, std::pair<size_t, bool> start, std::vector<std::pair<size_t, bool>>& unitigEnd, std::vector<std::pair<size_t, bool>>& unitigStart)
{
	assert(belongsToUnitig[start.first].first == std::numeric_limits<size_t>::max());
	start = reverse(start);
	auto pos = start;
	while (edges.getEdges(pos).size() == 1)
	{
		std::pair<size_t, bool> next = edges.getEdges(pos)[0];
		if (next.first == pos.first) break;
		if (next == start) break;
		if (edges.getEdges(reverse(next)).size() != 1) break;
		pos = next;
		assert(belongsToUnitig[pos.first].first == std::numeric_limits<size_t>::max());
	}
	unitigStart.push_back(pos);
	pos = reverse(pos);
	start = pos;
	assert(belongsToUnitig[pos.first].first == std::numeric_limits<size_t>::max());
	belongsToUnitig[pos.first] = std::make_pair(unitigEnd.size(), pos.second);
	while (edges.getEdges(pos).size() == 1)
	{
		std::pair<size_t, bool> next = edges.getEdges(pos)[0];
		if (edges.getEdges(reverse(next)).size() != 1) break;
		if (next == start) break;
		if (next.first == pos.first) break;
		pos = next;
		assert(belongsToUnitig[pos.first].first == std::numeric_limits<size_t>::max());
		belongsToUnitig[pos.first] = std::make_pair(unitigEnd.size(), pos.second);
	}
	unitigEnd.push_back(pos);
}

void removeArtefactEdges(std::vector<std::vector<size_t>>& outEdges, std::vector<std::vector<size_t>>& inEdges)
{
	std::vector<std::pair<size_t, size_t>> removeEdges;
	for (size_t i = 0; i < outEdges.size(); i++)
	{
		if (outEdges[i].size() != 2) continue;
		if (outEdges[outEdges[i][0]].size() == 1 && outEdges[outEdges[i][0]][0] == outEdges[i][1])
		{
			removeEdges.emplace_back(i, outEdges[i][0]);
		}
		if (outEdges[outEdges[i][1]].size() == 1 && outEdges[outEdges[i][1]][0] == outEdges[i][0])
		{
			removeEdges.emplace_back(i, outEdges[i][1]);
		}
	}
	if (removeEdges.size() == 0) return;
	for (auto edge : removeEdges)
	{
		for (size_t i = outEdges[edge.first].size()-1; i < outEdges[edge.first].size(); i--)
		{
			if (outEdges[edge.first][i] == edge.second)
			{
				std::swap(outEdges[edge.first][i], outEdges[edge.first].back());
				outEdges[edge.first].pop_back();
			}
		}
	}
	inEdges.clear();
	inEdges.resize(outEdges.size());
	for (size_t i = 0; i < outEdges.size(); i++)
	{
		for (size_t j = 0; j < outEdges[i].size(); j++)
		{
			inEdges[outEdges[i][j]].push_back(i);
		}
	}
}

template <typename F>
void iterateKmerChunkhashes(HashList& hashList, const size_t kmerSize, const phmap::flat_hash_map<uint64_t, size_t>& hashCounts, const ReadInfo& read, const SequenceCharType& seq, const SequenceLengthType& poses, const std::string& rawSeq, const std::vector<size_t>& positions, const std::vector<HashType>& hashes, const MatchIndex& matchIndex, F callback)
{
	std::vector<uint64_t> fwHashChunks;
	std::vector<uint32_t> fwStarts;
	std::vector<uint32_t> fwEnds;
	std::vector<uint64_t> bwHashChunks;
	std::vector<uint32_t> bwStarts;
	std::vector<uint32_t> bwEnds;
	matchIndex.iterateHashesFromRead(rawSeq, [&fwHashChunks, &fwStarts, &fwEnds, &hashCounts](uint64_t hash, size_t start, size_t end)
	{
		if (hashCounts.at(hash) < 2) return;
		fwHashChunks.push_back(hash);
		fwStarts.push_back(start);
		fwEnds.push_back(end);
	});
	matchIndex.iterateHashesFromRead(CommonUtils::ReverseComplement(rawSeq), [&bwHashChunks, &bwStarts, &bwEnds, &hashCounts, &rawSeq](uint64_t hash, size_t start, size_t end)
	{
		if (hashCounts.at(hash) < 2) return;
		bwHashChunks.push_back(hash);
		bwStarts.push_back(rawSeq.size()-1-end);
		bwEnds.push_back(rawSeq.size()-1-start);
	});
	std::reverse(bwHashChunks.begin(), bwHashChunks.end());
	std::reverse(bwStarts.begin(), bwStarts.end());
	std::reverse(bwEnds.begin(), bwEnds.end());
	for (size_t i = 1; i < fwHashChunks.size(); i++)
	{
		assert(fwStarts[i] >= fwStarts[i-1]);
	}
	for (size_t i = 1; i < bwHashChunks.size(); i++)
	{
		assert(bwStarts[i] >= bwStarts[i-1]);
	}
	assert(positions.size() == hashes.size());
	size_t fwStart = 0;
	size_t fwEnd = 0;
	size_t bwStart = 0;
	size_t bwEnd = 0;
	for (size_t i = 0; i < positions.size(); i++)
	{
		const auto pos = positions[i];
		const HashType fwHash = hashes[i];
		std::pair<size_t, bool> current = hashList.addNode(fwHash);
		while (fwStart < fwEnds.size() && fwEnds[fwStart] < poses[positions[i]+kmerSize]) fwStart += 1;
		fwEnd = std::max(fwEnd, fwStart);
		while (fwEnd < fwStarts.size() && fwStarts[fwEnd] < poses[positions[i]]) fwEnd += 1;
		while (bwStart < bwEnds.size() && bwEnds[bwStart] < poses[positions[i]+kmerSize]) bwStart += 1;
		bwEnd = std::max(bwEnd, bwStart);
		while (bwEnd < bwStarts.size() && bwStarts[bwEnd] < poses[positions[i]]) bwEnd += 1;
		callback(current, read.readName, poses[pos], fwHashChunks, fwStart, fwEnd, bwHashChunks, bwStart, bwEnd);
	}
}

template <typename F1, typename F2>
void iterateKmerChunkhashesFromFiles(const std::vector<std::string>& readFiles, const ReadpartIterator& partIterator, const MatchIndex& matchIndex, HashList& hashlist, const size_t kmerSize, const phmap::flat_hash_map<uint64_t, size_t>& hashCounts, F1 callback, F2 nextReadCallback)
{
	for (auto file : readFiles)
	{
		std::mutex indexMutex;
		partIterator.iterateHashes([callback, nextReadCallback, &hashlist, kmerSize, &hashCounts, &matchIndex](const ReadInfo& read, const SequenceCharType& seq, const SequenceLengthType& poses, const std::string& rawSeq, const std::vector<size_t>& positions, const std::vector<HashType>& hashes)
		{
			nextReadCallback();
			iterateKmerChunkhashes(hashlist, kmerSize, hashCounts, read, seq, poses, rawSeq, positions, hashes, matchIndex, callback);
		});
	}
}

uint64_t find(phmap::flat_hash_map<uint64_t, uint64_t>& parent, uint64_t key)
{
	if (parent.count(key) == 0)
	{
		parent[key] = key;
		return key;
	}
	while (parent.count(key) == 1 && parent.at(key) != key)
	{
		key = parent.at(key);
	}
	parent[key] = key;
	return key;
}

void merge(phmap::flat_hash_map<uint64_t, uint64_t>& parent, uint64_t left, uint64_t right)
{
	left = find(parent, left);
	right = find(parent, right);
	assert(parent.at(left) == left);
	assert(parent.at(right) == right);
	parent[right] = left;
}

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
	MatchIndex matchIndex { k, numWindows, windowSize };
	ReadStorage storage;
	ReadpartIterator partIterator { k, k-5, ErrorMasking::No, numThreads, readFiles, false, "" };
	phmap::flat_hash_map<uint64_t, size_t> hashCounts;
	std::mutex indexMutex;
	iterateReadsMultithreaded(readFiles, numThreads, [&matchIndex, &indexMutex, &hashCounts](ReadInfo readInfo, const std::string& sequence)
	{
		std::vector<uint64_t> hashes;
		matchIndex.iterateHashesFromRead(sequence, [&hashes](uint64_t hash, size_t start, size_t end) {
			hashes.push_back(hash);
		});
		matchIndex.iterateHashesFromRead(CommonUtils::ReverseComplement(sequence), [&hashes](uint64_t hash, size_t start, size_t end) {
			hashes.push_back(hash);
		});
		std::lock_guard<std::mutex> lock { indexMutex };
		for (auto hash : hashes)
		{
			hashCounts[hash] += 1;
		}
	});
	HashList hashList { k };
	std::vector<phmap::flat_hash_map<uint64_t, uint64_t>> parent;
	iterateKmerChunkhashesFromFiles(readFiles, partIterator, matchIndex, hashList, k, hashCounts, [&parent](std::pair<size_t, bool> kmer, const ReadName& readName, const size_t kmerPos, const std::vector<uint64_t>& fwHashes, size_t fwStart, size_t fwEnd, const std::vector<uint64_t>& bwHashes, size_t bwStart, size_t bwEnd)
	{
		if (kmer.first == parent.size()) parent.emplace_back();
		uint64_t anyHash = 0;
		if (fwEnd > fwStart)
		{
			anyHash = fwHashes[fwStart];
		}
		else if (bwEnd > bwStart)
		{
			anyHash = bwHashes[bwStart];
		}
		else
		{
			return;
		}
		for (size_t i = fwStart; i < fwEnd; i++)
		{
			merge(parent[kmer.first], fwHashes[i], anyHash);
		}
		for (size_t i = bwStart; i < bwEnd; i++)
		{
			merge(parent[kmer.first], bwHashes[i], anyHash);
		}
	}, [](){});
	std::vector<size_t> chunkKmerCoverage;
	std::vector<phmap::flat_hash_map<uint64_t, size_t>> pairToChunkKmer;
	pairToChunkKmer.resize(parent.size());
	for (size_t i = 0; i < parent.size(); i++)
	{
		phmap::flat_hash_set<uint64_t> clusters;
		phmap::flat_hash_set<uint64_t> keys;
		for (auto pair : parent[i])
		{
			keys.insert(pair.first);
		}
		for (auto key : keys)
		{
			clusters.insert(find(parent[i], key));
		}
		for (auto key : clusters)
		{
			pairToChunkKmer[i][key] = chunkKmerCoverage.size();
			chunkKmerCoverage.emplace_back(0);
		}
	}
	std::cerr << chunkKmerCoverage.size() << " chunk-kmers unfiltered" << std::endl;
	MostlySparse2DHashmap<uint8_t, size_t> edgeCoverages;
	edgeCoverages.resize(chunkKmerCoverage.size());
	std::pair<size_t, bool> last { std::numeric_limits<size_t>::max(), true };
	iterateKmerChunkhashesFromFiles(readFiles, partIterator, matchIndex, hashList, k, hashCounts, [&parent, &chunkKmerCoverage, &edgeCoverages, &last, &pairToChunkKmer](std::pair<size_t, bool> kmer, const ReadName& readName, const size_t kmerPos, const std::vector<uint64_t>& fwHashes, size_t fwStart, size_t fwEnd, const std::vector<uint64_t>& bwHashes, size_t bwStart, size_t bwEnd)
	{
		uint64_t anyHash = 0;
		if (fwEnd > fwStart)
		{
			anyHash = fwHashes[fwStart];
		}
		else if (bwEnd > bwStart)
		{
			anyHash = bwHashes[bwStart];
		}
		else
		{
			last.first = std::numeric_limits<size_t>::max();
			return;
		}
		uint64_t hashchunk = find(parent[kmer.first], anyHash);
		std::pair<size_t, bool> current { pairToChunkKmer[kmer.first].at(hashchunk), kmer.second };
		chunkKmerCoverage[current.first] += 1;
		if (last.first != std::numeric_limits<size_t>::max())
		{
			auto key = canon(last, current);
			size_t old = 0;
			if (edgeCoverages.hasValue(key.first, key.second)) old = edgeCoverages.get(key.first, key.second);
			edgeCoverages.set(key.first, key.second, old+1);
		}
		last = current;
	}, [&last]() { last.first = std::numeric_limits<size_t>::max(); });
	std::vector<bool> keptNodes;
	SparseEdgeContainer keptEdges;
	keptNodes.resize(chunkKmerCoverage.size(), false);
	keptEdges.resize(chunkKmerCoverage.size());
	size_t numKeptNodes = 0;
	size_t numKeptEdges = 0;
	for (size_t i = 0; i < chunkKmerCoverage.size(); i++)
	{
		if (chunkKmerCoverage[i] < minCoverage) continue;
		keptNodes[i] = true;
		numKeptNodes += 1;
	}
	std::cerr << numKeptNodes << " chunk-kmers filtered" << std::endl;
	for (size_t i = 0; i < chunkKmerCoverage.size(); i++)
	{
		if (!keptNodes[i]) continue;
		for (auto value : edgeCoverages.getValues(std::make_pair(i, true)))
		{
			if (value.second < minCoverage) continue;
			if (!keptNodes[value.first.first]) continue;
			keptEdges.addEdge(std::make_pair(i, true), value.first);
			keptEdges.addEdge(reverse(value.first), std::make_pair(i, false));
			numKeptEdges += 1;
		}
		for (auto value : edgeCoverages.getValues(std::make_pair(i, false)))
		{
			if (value.second < minCoverage) continue;
			if (!keptNodes[value.first.first]) continue;
			keptEdges.addEdge(std::make_pair(i, false), value.first);
			keptEdges.addEdge(reverse(value.first), std::make_pair(i, true));
			numKeptEdges += 1;
		}
	}
	std::cerr << numKeptEdges << " chunk-kmers edges filtered" << std::endl;
	std::vector<std::pair<size_t, bool>> belongsToUnitig;
	belongsToUnitig.resize(chunkKmerCoverage.size(), std::make_pair(std::numeric_limits<size_t>::max(), false));
	std::vector<std::pair<size_t, bool>> unitigEnd;
	std::vector<std::pair<size_t, bool>> unitigStart;
	for (size_t i = 0; i < belongsToUnitig.size(); i++)
	{
		if (!keptNodes[i]) continue;
		if (belongsToUnitig[i].first != std::numeric_limits<size_t>::max()) continue;
		getUnitig(belongsToUnitig, keptNodes, keptEdges, std::make_pair(i, true), unitigEnd, unitigStart);
		assert(belongsToUnitig[i].first != std::numeric_limits<size_t>::max());
	}
	std::vector<size_t> unitigCoverageSum;
	std::vector<size_t> unitigCoverageDiv;
	unitigCoverageDiv.resize(unitigEnd.size(), 0);
	unitigCoverageSum.resize(unitigEnd.size(), 0);
	for (size_t i = 0; i < belongsToUnitig.size(); i++)
	{
		if (!keptNodes[i]) continue;
		assert(belongsToUnitig[i].first != std::numeric_limits<size_t>::max());
		assert(belongsToUnitig[i].first < unitigEnd.size());
		unitigCoverageSum[belongsToUnitig[i].first] += chunkKmerCoverage[i];
		unitigCoverageDiv[belongsToUnitig[i].first] += 1;
	}
	size_t unitigEdges = 0;
	size_t unitigCount = 0;
	std::ofstream graphfile { "graph.gfa" };
	for (size_t i = 0; i < unitigEnd.size(); i++)
	{
		// if (unitigCoverageDiv[i] == 0) continue;
		size_t size = k + (unitigCoverageDiv[i]-1)*(k-5)/2;
		size_t coverage = unitigCoverageSum[i] / unitigCoverageDiv[i];
		// if (coverage < minCoverage) continue;
		unitigCount += 1;
		graphfile << "S\t" << i << "\t*\tLN:i:" << size << "\tll:f:" << coverage << "\tFC:i:" << (coverage*size) << std::endl;
		for (auto edge : keptEdges.getEdges(unitigEnd[i]))
		{
			auto targetUnitig = belongsToUnitig[edge.first];
			if (!edge.second) targetUnitig = reverse(targetUnitig);
			graphfile << "L\t" << i << "\t+\t" << targetUnitig.first << "\t" << (targetUnitig.second ? "+" : "-") << "\t" << int((k-5)/2) << "M" << std::endl;
			unitigEdges += 1;
		}
		for (auto edge : keptEdges.getEdges(unitigStart[i]))
		{
			auto targetUnitig = belongsToUnitig[edge.first];
			if (!edge.second) targetUnitig = reverse(targetUnitig);
			graphfile << "L\t" << i << "\t-\t" << targetUnitig.first << "\t" << (targetUnitig.second ? "+" : "-") << "\t" << int((k-5)/2) << "M" << std::endl;
			unitigEdges += 1;
		}
	}
	std::cerr << unitigCount << " unitigs" << std::endl;
	std::cerr << unitigEdges << " unitig-edges" << std::endl;
}
