#include <iostream>
#include <limits>
#include <cassert>
#include "UnitigStorage.h"

HashType reverseHash(HashType fwHash)
{
	return (fwHash << 64) + (fwHash >> 64);
}

UnitigStorage::UnitigStorage(size_t k) :
	kmerSize(k)
{
}

std::pair<size_t, bool> UnitigStorage::getNode(HashType fwHash)
{
	auto revhash = reverseHash(fwHash);
	std::pair<size_t, bool> kmer;
	if (revhash < fwHash)
	{
		std::lock_guard lock { addMutex };
		auto found = hashToNode.find(revhash);
		if (found != hashToNode.end())
		{
			return std::make_pair(found->second, false);
		}
		kmer = std::make_pair(hashToNode.size(), false);
		hashToNode[revhash] = kmer.first;
		uniqueEdge.emplace_back(std::numeric_limits<uint32_t>::max());
		return kmer;
	}
	else
	{
		std::lock_guard lock { addMutex };
		auto found = hashToNode.find(fwHash);
		if (found != hashToNode.end())
		{
			return std::make_pair(found->second, true);
		}
		kmer = std::make_pair(hashToNode.size(), true);
		hashToNode[fwHash] = kmer.first;
		uniqueEdge.emplace_back(std::numeric_limits<uint32_t>::max());
		return kmer;
	}
}

std::pair<size_t, bool> UnitigStorage::getNodeOrNull(HashType fwHash) const
{
	auto revhash = reverseHash(fwHash);
	if (revhash < fwHash)
	{
		auto found = hashToNode.find(revhash);
		if (found == hashToNode.end()) return std::make_pair(std::numeric_limits<size_t>::max(), true);
		return std::make_pair(found->second, false);
	}
	auto found = hashToNode.find(fwHash);
	if (found == hashToNode.end()) return std::make_pair(std::numeric_limits<size_t>::max(), true);
	return std::make_pair(found->second, true);
}

std::tuple<size_t, size_t, std::vector<std::pair<size_t, bool>>> UnitigStorage::getPath(const std::vector<HashType>& hashes) const
{
	std::vector<std::pair<size_t, bool>> unitigPath;
	size_t leftClip = 0;
	size_t rightClip = 0;
	assert(hashes.size() > 0);
	std::pair<size_t, bool> currentUnitig { std::numeric_limits<size_t>::max(), true };
	size_t currentUnitigOffset = std::numeric_limits<size_t>::max()-1;
	for (size_t i = 0; i < hashes.size(); i++)
	{
		std::pair<size_t, bool> previousUnitig = currentUnitig;
		size_t previousUnitigOffset = currentUnitigOffset;
		std::pair<size_t, bool> kmer = getNodeOrNull(hashes[i]);
		assert(kmer.first != std::numeric_limits<size_t>::max());
		currentUnitig = numToPair(kmerPositionInUnitig[kmer.first].first);
		currentUnitigOffset = kmerPositionInUnitig[kmer.first].second;
		if (!kmer.second)
		{
			currentUnitig = reverse(currentUnitig);
			currentUnitigOffset = unitigLength[currentUnitig.first] - 1 - currentUnitigOffset;
		}
		if (currentUnitig != previousUnitig || currentUnitigOffset != previousUnitigOffset+1)
		{
			assert(unitigPath.size() == 0 || previousUnitigOffset == unitigLength[previousUnitig.first]-1);
			assert(unitigPath.size() == 0 || currentUnitigOffset == 0);
			unitigPath.push_back(currentUnitig);
		}
		if (i == 0) leftClip = currentUnitigOffset;
		if (i == hashes.size()-1) rightClip = unitigLength[currentUnitig.first] - 1 - currentUnitigOffset;
	}
	assert(unitigPath.size() >= 1);
	return std::make_tuple(leftClip, rightClip, unitigPath);
}

void UnitigStorage::addEdge(std::pair<size_t, bool> from, std::pair<size_t, bool> to, size_t overlap)
{
	if (uniqueEdge[from] == std::numeric_limits<uint32_t>::max())
	{
		uniqueEdge[from] = pairToNum(to);
	}
	else if (uniqueEdge[from] != pairToNum(to))
	{
		uniqueEdge[from] = std::numeric_limits<uint32_t>::max()-1;
	}
	if (uniqueEdge[reverse(to)] == std::numeric_limits<uint32_t>::max())
	{
		uniqueEdge[reverse(to)] = pairToNum(reverse(from));
	}
	else if (uniqueEdge[reverse(to)] != pairToNum(reverse(from)))
	{
		uniqueEdge[reverse(to)] = std::numeric_limits<uint32_t>::max()-1;
	}
	auto c = canon(from, to);
	auto key = std::make_pair(pairToNum(c.first), pairToNum(c.second));
	assert(kmerOverlap.count(key) == 0 || kmerOverlap.at(key) == overlap);
	kmerOverlap[key] = overlap;
}

std::string UnitigStorage::getSequence(const std::vector<std::pair<size_t, bool>>& path, size_t leftClip, size_t rightClip) const
{
	std::string result;
	for (size_t i = 0; i < path.size(); i++)
	{
		std::string add;
		size_t leftClipBp = 0;
		size_t rightClipBp = 0;
		if (i == 0)
		{
			if (path[i].second)
			{
				leftClipBp = kmerBpOffsetInsideUnitig[path[i].first][leftClip];
			}
			else
			{
				rightClipBp = kmerBpOffsetInsideUnitig[path[i].first].back() - kmerBpOffsetInsideUnitig[path[i].first][kmerBpOffsetInsideUnitig[path[i].first].size()-1-leftClip];
			}
		}
		else
		{
			auto c = canon(path[i-1], path[i]);
			size_t overlap = unitigEdgeOverlap.at(std::make_pair(pairToNum(c.first), pairToNum(c.second)));
			if (path[i].second)
			{
				leftClipBp = overlap;
			}
			else
			{
				rightClipBp = overlap;
			}
		}
		if (i == path.size()-1)
		{
			if (path[i].second)
			{
				rightClipBp = kmerBpOffsetInsideUnitig[path[i].first].back() - kmerBpOffsetInsideUnitig[path[i].first][kmerBpOffsetInsideUnitig[path[i].first].size()-1-rightClip];
			}
			else
			{
				leftClipBp = kmerBpOffsetInsideUnitig[path[i].first][rightClip];
			}
		}
		add = unitigSequences[path[i].first].substr(leftClipBp, unitigSequences[path[i].first].size() - leftClipBp - rightClipBp);
		if (!path[i].second) add = revCompRaw(add);
		result += add;
	}
	return result;
}

void UnitigStorage::buildUnitigGraph()
{
	assert(hashToNode.size() == uniqueEdge.size());
	kmerPositionInUnitig.resize(hashToNode.size(), std::make_pair(std::numeric_limits<uint32_t>::max(), std::numeric_limits<uint32_t>::max()));
	for (size_t i = 0; i < kmerPositionInUnitig.size(); i++)
	{
		if (kmerPositionInUnitig[i].first != std::numeric_limits<uint32_t>::max()) continue;
		std::pair<size_t, bool> fw { i, true };
		std::pair<size_t, bool> bw { i, false };
		if (uniqueEdge[fw] >= std::numeric_limits<uint32_t>::max()-1)
		{
			beginUnitig(bw);
			continue;
		}
		if (uniqueEdge[bw] >= std::numeric_limits<uint32_t>::max()-1)
		{
			beginUnitig(fw);
			continue;
		}
		if (uniqueEdge[fw] < std::numeric_limits<uint32_t>::max()-1 && uniqueEdge[reverse(numToPair(uniqueEdge[fw]))] != pairToNum(bw))
		{
			beginUnitig(bw);
			continue;
		}
		if (uniqueEdge[bw] < std::numeric_limits<uint32_t>::max()-1 && uniqueEdge[reverse(numToPair(uniqueEdge[bw]))] != pairToNum(fw))
		{
			beginUnitig(fw);
			continue;
		}
	}
	for (size_t i = 0; i < kmerPositionInUnitig.size(); i++)
	{
		if (kmerPositionInUnitig[i].first != std::numeric_limits<uint32_t>::max()) continue;
		beginUnitig(std::make_pair(i, true));
	}
	{
		VectorWithDirection<uint32_t> tmp;
		std::swap(tmp, uniqueEdge);
	}
	kmerSequenceLoaded.resize(hashToNode.size(), false);
	assert(kmerSequenceLoaded.size() == hashToNode.size());
	assert(unitigSequences.size() == unitigLength.size());
	std::vector<bool> unitigTip;
	unitigTip.resize(kmerPositionInUnitig.size(), false);
	size_t foundTips = 0;
	for (size_t i = 0; i < kmerPositionInUnitig.size(); i++)
	{
		assert(kmerPositionInUnitig[i].first != std::numeric_limits<uint32_t>::max());
		if (kmerPositionInUnitig[i].second == 0 || kmerPositionInUnitig[i].second == unitigLength[numToPair(kmerPositionInUnitig[i].first).first]-1)
		{
			unitigTip[i] = true;
			foundTips += 1;
			if (unitigLength[numToPair(kmerPositionInUnitig[i].first).first] == 1) foundTips += 1;
		}
	}
	assert(foundTips == unitigLength.size()*2);
	for (auto pair : kmerOverlap)
	{
		std::pair<size_t, bool> kmerFrom = numToPair(pair.first.first);
		if (!unitigTip[kmerFrom.first]) continue;
		std::pair<size_t, bool> kmerTo = numToPair(pair.first.second);
		if (!unitigTip[kmerTo.first]) continue;
		std::pair<size_t, bool> unitigFrom, unitigTo;
		unitigFrom = numToPair(kmerPositionInUnitig[kmerFrom.first].first);
		size_t fromPos = kmerPositionInUnitig[kmerFrom.first].second;
		if (!kmerFrom.second)
		{
			unitigFrom = reverse(unitigFrom);
			fromPos = unitigLength[unitigFrom.first] - 1 - fromPos;
		}
		assert(fromPos == 0 || fromPos == unitigLength[unitigFrom.first]-1);
		if (fromPos != unitigLength[unitigFrom.first]-1) continue;
		unitigTo = numToPair(kmerPositionInUnitig[kmerTo.first].first);
		size_t toPos = kmerPositionInUnitig[kmerTo.first].second;
		if (!kmerTo.second)
		{
			unitigTo = reverse(unitigTo);
			toPos = unitigLength[unitigTo.first] - 1 - toPos;
		}
		assert(toPos == 0 || toPos == unitigLength[unitigTo.first]-1);
		if (toPos != 0) continue;
		auto c = canon(unitigFrom, unitigTo);
		auto key = std::make_pair(pairToNum(c.first), pairToNum(c.second));
		assert(unitigEdgeOverlap.count(key) == 0);
		unitigEdgeOverlap[key] = pair.second;
	}
	{
		phmap::flat_hash_map<std::pair<uint32_t, uint32_t>, uint16_t> tmp;
		std::swap(tmp, kmerOverlap);
	}
	assert(kmerBpOffsetInsideUnitig.size() == unitigLength.size());
	assert(unitigSequences.size() == unitigLength.size());
}

void UnitigStorage::beginUnitig(std::pair<size_t, bool> start)
{
	std::pair<size_t, bool> pos = start;
	size_t thisUnitig = unitigLength.size();
	unitigLength.push_back(0);
	unitigSequences.emplace_back();
	kmerBpOffsetInsideUnitig.emplace_back();
	size_t bpOffset = 0;
	std::vector<size_t> swapThese;
	while (true)
	{
		assert(kmerPositionInUnitig[pos.first].first == std::numeric_limits<uint32_t>::max());
		kmerPositionInUnitig[pos.first] = std::make_pair(pairToNum(std::make_pair(thisUnitig, pos.second)), unitigLength[thisUnitig]);
		kmerBpOffsetInsideUnitig[thisUnitig].push_back(bpOffset);
		if (!pos.second) swapThese.push_back(pos.first);
		unitigLength[thisUnitig] += 1;
		if (uniqueEdge[pos] >= std::numeric_limits<uint32_t>::max()-1) break;
		std::pair<size_t, bool> next = numToPair(uniqueEdge[pos]);
		if (uniqueEdge[reverse(next)] != pairToNum(reverse(pos))) break;
		if (next.first == pos.first) break;
		if (next.first == start.first) break;
		auto c = canon(pos, next);
		auto key = std::make_pair(pairToNum(c.first), pairToNum(c.second));
		bpOffset += kmerSize - kmerOverlap.at(key);
		pos = next;
	}
	for (auto index : swapThese)
	{
		assert(numToPair(kmerPositionInUnitig[index].first).first == thisUnitig);
		kmerPositionInUnitig[index].second = unitigLength[thisUnitig] - 1 - kmerPositionInUnitig[index].second;
	}
	assert(kmerBpOffsetInsideUnitig[thisUnitig].size() == unitigLength[thisUnitig]);
	assert(kmerPositionInUnitig[start.first].second == 0 || kmerPositionInUnitig[start.first].second == unitigLength[thisUnitig] - 1);
	assert(kmerPositionInUnitig[pos.first].second == 0 || kmerPositionInUnitig[pos.first].second == unitigLength[thisUnitig] - 1);
	assert(kmerPositionInUnitig[start.first].second == 0 || kmerPositionInUnitig[start.first].second == unitigLength[numToPair(kmerPositionInUnitig[start.first].first).first] - 1);
	assert(kmerPositionInUnitig[pos.first].second == 0 || kmerPositionInUnitig[pos.first].second == unitigLength[numToPair(kmerPositionInUnitig[pos.first].first).first] - 1);
	unitigSequences[thisUnitig].resize(kmerBpOffsetInsideUnitig[thisUnitig].back()+kmerSize);
}

void UnitigStorage::addKmerSequence(std::pair<size_t, bool> kmer, const std::string& seq, size_t start, size_t end)
{
	assert(end-start == kmerSize);
	assert(kmer.first < kmerSequenceLoaded.size());
	if (kmerSequenceLoaded[kmer.first]) return;
	std::pair<size_t, bool> unitig = numToPair(kmerPositionInUnitig[kmer.first].first);
	size_t offset = kmerPositionInUnitig[kmer.first].second;
	bool revcomp = !kmer.second;
	if (!unitig.second)
	{
		revcomp = !revcomp;
		offset = unitigLength[unitig.first]-1-offset;
	}
	for (size_t i = 0; i < kmerSize; i++)
	{
		char c = seq[start+i];
		assert(c == 'A' || c == 'C' || c == 'G' || c == 'T');
		switch(c)
		{
			case 'A':
				c = 0;
				break;
			case 'C':
				c = 1;
				break;
			case 'G':
				c = 2;
				break;
			case 'T':
				c = 3;
				break;
		}
		size_t bpIndex = kmerBpOffsetInsideUnitig[unitig.first][offset]+i;
		if (revcomp)
		{
			bpIndex = kmerBpOffsetInsideUnitig[unitig.first][offset]+kmerSize-1-i;
			c = 3-c;
			switch(c)
			{
				case 'A':
					c = 'T';
					break;
				case 'C':
					c = 'G';
					break;
				case 'G':
					c = 'C';
					break;
				case 'T':
					c = 'A';
					break;
			}
		}
		unitigSequences[unitig.first].set(bpIndex, c);
	}
	kmerSequenceLoaded[kmer.first] = true;
}

void UnitigStorage::finalizeSequences()
{
	{
		phmap::flat_hash_map<HashType, uint32_t> tmp;
		std::swap(tmp, hashToNode);
	}
	for (size_t i = 0; i < kmerSequenceLoaded.size(); i++)
	{
		assert(kmerSequenceLoaded[i]);
	}
	for (size_t i = 0; i < unitigLength.size(); i++)
	{
		assert(kmerBpOffsetInsideUnitig[i].size() == unitigLength[i]);
		assert(unitigSequences[i].size() == kmerBpOffsetInsideUnitig[i].back() + kmerSize);
	}
	{
		std::vector<bool> tmp;
		std::swap(tmp, kmerSequenceLoaded);
	}
}

std::pair<size_t, bool> UnitigStorage::numToPair(size_t num) const
{
	return std::make_pair(num / 2, num % 2);
}

size_t UnitigStorage::pairToNum(std::pair<size_t, bool> p) const
{
	return p.first * 2 + (p.second ? 1 : 0);
}

size_t UnitigStorage::numHashes() const
{
	if (hashToNode.size() > 0) return hashToNode.size();
	return kmerPositionInUnitig.size();
}

size_t UnitigStorage::numUnitigs() const
{
	return unitigSequences.size();
}

size_t UnitigStorage::totalBps() const
{
	size_t result = 0;
	for (size_t i = 0; i < unitigSequences.size(); i++)
	{
		result += unitigSequences[i].size();
	}
	return result;
}

size_t UnitigStorage::unitigMinusEdgeLength(const std::pair<size_t, bool> fromUnitig, const std::pair<size_t, bool> toUnitig) const
{
	size_t result = unitigSequences[toUnitig.first].size();
	auto c = canon(fromUnitig, toUnitig);
	size_t overlap = unitigEdgeOverlap.at(std::make_pair(pairToNum(c.first), pairToNum(c.second)));
	assert(overlap < result);
	return result - overlap;
}
