#include <mutex>
#include <iostream>
#include <vector>
#include <string>
#include "ReadStorage.h"
#include "MatchIndex.h"

class Path
{
public:
	std::vector<std::pair<size_t, bool>> nodes;
	std::vector<size_t> lengths;
	double startFractionalClip;
	double endFractionalClip;
	size_t lengthIncludingClippedOutParts() const
	{
		size_t result = 0;
		for (auto l : lengths) result += l;
		return result;
	};
};

Path getPath(const std::vector<size_t>& readLengths, const std::vector<std::vector<std::pair<std::pair<size_t, bool>, size_t>>>& readPaths, const size_t index, const bool fw, const size_t start, const size_t end)
{
	assert(end > start);
	assert(end <= readLengths[index]);
	if (!fw)
	{
		Path partialResult = getPath(readLengths, readPaths, index, true, readLengths[index]-end, readLengths[index]-start);
		std::reverse(partialResult.nodes.begin(), partialResult.nodes.end());
		for (size_t i = 0; i < partialResult.nodes.size(); i++)
		{
			partialResult.nodes[i].second = !partialResult.nodes[i].second;
		}
		std::reverse(partialResult.lengths.begin(), partialResult.lengths.end());
		std::swap(partialResult.startFractionalClip, partialResult.endFractionalClip);
		return partialResult;
	}
	std::vector<std::pair<size_t, bool>> pathnodes;
	std::vector<size_t> lengths;
	double startFractionalClip;
	double endFractionalClip;
	bool startClipFound = false;
	bool endClipFound = false;
	size_t nodeStartPos = 0;
	for (size_t i = 0; i < readPaths[index].size(); i++)
	{
		size_t nodeEndPos = nodeStartPos + readPaths[index][i].second;
		if (nodeEndPos > start && nodeStartPos <= start)
		{
			assert(!startClipFound);
			startFractionalClip = (double)(start - nodeStartPos) / (double)(nodeEndPos - nodeStartPos);
			startClipFound = true;
		}
		assert(nodeStartPos < end);
		if (nodeEndPos > start)
		{
			pathnodes.push_back(readPaths[index][i].first);
			lengths.push_back(readPaths[index][i].second);
		}
		if (nodeEndPos >= end)
		{
			assert(!endClipFound);
			endFractionalClip = (double)(nodeEndPos - end) / (double)(nodeEndPos - nodeStartPos);
			endClipFound = true;
			break;
		}
		nodeStartPos = nodeEndPos;
	}
	assert(nodeStartPos <= readLengths[index]);
	assert(startClipFound);
	assert(endClipFound);
	assert(pathnodes.size() >= 1);
	assert(lengths.size() == pathnodes.size());
	assert(lengths[0] > 0);
	assert(lengths.back() > 0);
	Path result;
	result.nodes = pathnodes;
	result.lengths = lengths;
	result.startFractionalClip = startFractionalClip;
	result.endFractionalClip = endFractionalClip;
	return result;
}

std::vector<std::pair<Path, Path>> getPartialPathsAtNonequalSites(const Path& leftPath, const Path& rightPath)
{
	std::unordered_set<size_t> nodesInLeft;
	std::unordered_set<size_t> nodesInRight;
	for (auto node : leftPath.nodes)
	{
		nodesInLeft.emplace(node.first);
	}
	for (auto node : rightPath.nodes)
	{
		nodesInRight.emplace(node.first);
	}
	std::vector<size_t> leftEqualIndices;
	std::vector<size_t> rightEqualIndices;
	for (size_t i = 0; i < leftPath.nodes.size(); i++)
	{
		if (nodesInRight.count(leftPath.nodes[i].first) == 0) continue;
		leftEqualIndices.push_back(i);
	}
	for (size_t i = 0; i < rightPath.nodes.size(); i++)
	{
		if (nodesInLeft.count(rightPath.nodes[i].first) == 0) continue;
		rightEqualIndices.push_back(i);
	}
	if (leftEqualIndices.size() == 0 && rightEqualIndices.size() == 0)
	{
		return std::vector<std::pair<Path, Path>> { std::make_pair(leftPath, rightPath) };
	}
	if (leftEqualIndices.size() != rightEqualIndices.size())
	{
		std::cerr << "could not merge paths, different equal index sizes" << std::endl;
		return std::vector<std::pair<Path, Path>> {};
	}
	for (size_t i = 0; i < leftEqualIndices.size(); i++)
	{
		if (leftPath.nodes[leftEqualIndices[i]].first != rightPath.nodes[rightEqualIndices[i]].first)
		{
			std::cerr << "could not merge paths, equal indices have different nodes" << std::endl;
			return std::vector<std::pair<Path, Path>> {};
		}
	}
	std::cerr << "partial matches at:";
	for (size_t i = 0; i < leftEqualIndices.size(); i++)
	{
		std::cerr << " " << leftEqualIndices[i] << "," << rightEqualIndices[i];
	}
	std::cerr << std::endl;
	std::vector<std::pair<Path, Path>> result;
	if (leftEqualIndices[0] > 0 && rightEqualIndices[0] > 0)
	{
		Path leftSubpath;
		leftSubpath.startFractionalClip = leftPath.startFractionalClip;
		leftSubpath.endFractionalClip = 0;
		leftSubpath.nodes.insert(leftSubpath.nodes.end(), leftPath.nodes.begin(), leftPath.nodes.begin() + leftEqualIndices[0]);
		leftSubpath.lengths.insert(leftSubpath.lengths.end(), leftPath.lengths.begin(), leftPath.lengths.begin() + leftEqualIndices[0]);
		Path rightSubpath;
		rightSubpath.startFractionalClip = rightPath.startFractionalClip;
		rightSubpath.endFractionalClip = 0;
		rightSubpath.nodes.insert(rightSubpath.nodes.end(), rightPath.nodes.begin(), rightPath.nodes.begin() + rightEqualIndices[0]);
		rightSubpath.lengths.insert(rightSubpath.lengths.end(), rightPath.lengths.begin(), rightPath.lengths.begin() + rightEqualIndices[0]);
		result.emplace_back(leftSubpath, rightSubpath);
	}
	for (size_t i = 1; i < leftEqualIndices.size(); i++)
	{
		if (leftEqualIndices[i-1]+1 == leftEqualIndices[i]) continue;
		if (rightEqualIndices[i-1]+1 == rightEqualIndices[i]) continue;
		Path leftSubpath;
		leftSubpath.startFractionalClip = 0;
		leftSubpath.endFractionalClip = 0;
		leftSubpath.nodes.insert(leftSubpath.nodes.end(), leftPath.nodes.begin() + leftEqualIndices[i-1]+1, leftPath.nodes.begin() + leftEqualIndices[i]);
		leftSubpath.lengths.insert(leftSubpath.lengths.end(), leftPath.lengths.begin() + leftEqualIndices[i-1]+1, leftPath.lengths.begin() + leftEqualIndices[i]);
		Path rightSubpath;
		rightSubpath.startFractionalClip = 0;
		rightSubpath.endFractionalClip = 0;
		rightSubpath.nodes.insert(rightSubpath.nodes.end(), rightPath.nodes.begin() + rightEqualIndices[i-1]+1, rightPath.nodes.begin() + rightEqualIndices[i]);
		rightSubpath.lengths.insert(rightSubpath.lengths.end(), rightPath.lengths.begin() + rightEqualIndices[i-1]+1, rightPath.lengths.begin() + rightEqualIndices[i]);
		result.emplace_back(leftSubpath, rightSubpath);
	}
	if (leftEqualIndices.back()+1 < leftPath.nodes.size() && rightEqualIndices.back()+1 < rightPath.nodes.size())
	{
		Path leftSubpath;
		leftSubpath.startFractionalClip = 0;
		leftSubpath.endFractionalClip = leftPath.endFractionalClip;
		leftSubpath.nodes.insert(leftSubpath.nodes.end(), leftPath.nodes.begin() + leftEqualIndices.back()+1, leftPath.nodes.end());
		leftSubpath.lengths.insert(leftSubpath.lengths.end(), leftPath.lengths.begin() + leftEqualIndices.back()+1, leftPath.lengths.end());
		Path rightSubpath;
		rightSubpath.startFractionalClip = 0;
		rightSubpath.endFractionalClip = rightPath.endFractionalClip;
		rightSubpath.nodes.insert(rightSubpath.nodes.end(), rightPath.nodes.begin() + rightEqualIndices.back()+1, rightPath.nodes.end());
		rightSubpath.lengths.insert(rightSubpath.lengths.end(), rightPath.lengths.begin() + rightEqualIndices.back()+1, rightPath.lengths.end());
		result.emplace_back(leftSubpath, rightSubpath);
	}
	return result;
}

void replaceReadpathNode(const std::vector<size_t>& readLengths, std::vector<std::vector<std::pair<std::pair<size_t, bool>, size_t>>>& readPaths, std::vector<std::vector<std::pair<size_t, size_t>>>& readsPerNode, const size_t node, const std::vector<std::pair<size_t, bool>>& newNodes, const size_t readIndex, const size_t index, const std::vector<double>& nodeFractionalLengths)
{
	assert(newNodes.size() >= 2);
	assert(nodeFractionalLengths.size() == newNodes.size());
	assert(readPaths[readIndex][index].first.first == node);
	std::vector<std::pair<std::pair<size_t, bool>, size_t>> newPath;
	newPath.insert(newPath.end(), readPaths[readIndex].begin(), readPaths[readIndex].begin()+index);
	size_t oldNodeLength = readPaths[readIndex][index].second;
	size_t lengthUsed = 0;
	if (readPaths[readIndex][index].first.second)
	{
		for (size_t i = 0; i < newNodes.size(); i++)
		{
			size_t lengthHere = nodeFractionalLengths[i] * oldNodeLength;
			newPath.emplace_back(newNodes[i], lengthHere);
			lengthUsed += lengthHere;
		}
	}
	else
	{
		for (size_t i = newNodes.size()-1; i < newNodes.size(); i--)
		{
			size_t lengthHere = nodeFractionalLengths[i] * oldNodeLength;
			newPath.emplace_back(newNodes[i], lengthHere);
			newPath.back().first.second = !newPath.back().first.second;
			lengthUsed += lengthHere;
		}
	}
	assert(lengthUsed <= oldNodeLength);
	while (lengthUsed < oldNodeLength)
	{
		for (size_t i = 0; i < newNodes.size() && lengthUsed < oldNodeLength; i++)
		{
			newPath[newPath.size()-newNodes.size()+i].second += 1;
			lengthUsed += 1;
		}
	}
	assert(lengthUsed == oldNodeLength);
	newPath.insert(newPath.end(), readPaths[readIndex].begin()+index+1, readPaths[readIndex].end());
	for (size_t i = index; i < readPaths[readIndex].size(); i++)
	{
		bool found = false;
		for (size_t j = 0; j < readsPerNode[readPaths[readIndex][i].first.first].size(); j++)
		{
			if (readsPerNode[readPaths[readIndex][i].first.first][j].first != readIndex) continue;
			if (readsPerNode[readPaths[readIndex][i].first.first][j].second != i) continue;
			assert(!found);
			std::swap(readsPerNode[readPaths[readIndex][i].first.first][j], readsPerNode[readPaths[readIndex][i].first.first].back());
			readsPerNode[readPaths[readIndex][i].first.first].pop_back();
			found = true;
		}
		assert(found);
	}
	for (size_t i = index; i < newPath.size(); i++)
	{
		readsPerNode[newPath[i].first.first].emplace_back(readIndex, i);
	}
	readPaths[readIndex] = newPath;
}

std::vector<std::pair<size_t, bool>> cutNode(const std::vector<size_t>& readLengths, std::vector<std::vector<std::pair<std::pair<size_t, bool>, size_t>>>& readPaths, std::vector<std::vector<std::pair<size_t, size_t>>>& readsPerNode, const size_t node, const std::vector<double>& nodeFractionalLengths)
{
	assert(nodeFractionalLengths.size() >= 2);
	std::vector<std::pair<size_t, bool>> newNodes;
	for (size_t i = 0; i < nodeFractionalLengths.size(); i++)
	{
		newNodes.emplace_back(readsPerNode.size(), true);
		readsPerNode.emplace_back();
	}
	std::cerr << "cut node " << node << " into " << nodeFractionalLengths.size() << " pieces (";
	for (auto node : newNodes)
	{
		std::cerr << (node.second ? ">" : "<") << node.first;
	}
	std::cerr << ")" << std::endl;
	std::vector<std::pair<size_t, size_t>> relevants = readsPerNode[node];
	//needs to be sorted last pos to first pos so replacement works, otherwise later poses are invalidated by earlier replacements
	std::sort(relevants.begin(), relevants.end(), [](auto left, auto right){ return left.second > right.second; });
	for (auto pos : relevants)
	{
		replaceReadpathNode(readLengths, readPaths, readsPerNode, node, newNodes, pos.first, pos.second, nodeFractionalLengths);
	}
	readsPerNode[node].clear();
	return newNodes;
}

Path cutPath(const std::vector<size_t>& readLengths, std::vector<std::vector<std::pair<std::pair<size_t, bool>, size_t>>>& readPaths, std::vector<std::vector<std::pair<size_t, size_t>>>& readsPerNode, const Path& path, const std::vector<double>& missingFractionalBreakpoints)
{
	if (missingFractionalBreakpoints.size() == 0) return path;
	phmap::flat_hash_map<size_t, std::vector<double>> fractionalCutsPerNode;
	size_t nodeStartPos = 0;
	size_t pathLength = path.lengthIncludingClippedOutParts();
	size_t fractionalIndex = 0;
	for (size_t i = 0; i < path.nodes.size(); i++)
	{
		size_t nodeEndPos = nodeStartPos + path.lengths[i];
		double nodeFractionalStart = (double)nodeStartPos / (double)pathLength;
		double nodeFractionalEnd = (double)nodeEndPos / (double)pathLength;
		while (fractionalIndex < missingFractionalBreakpoints.size() && nodeFractionalStart <= missingFractionalBreakpoints[fractionalIndex] && nodeFractionalEnd > missingFractionalBreakpoints[fractionalIndex])
		{
			double pos = missingFractionalBreakpoints[fractionalIndex] - nodeFractionalStart;
			pos /= nodeFractionalEnd - nodeFractionalStart;
			assert(pos >= 0);
			assert(pos <= 1);
			if (!path.nodes[i].second) pos = 1.0-pos;
			fractionalCutsPerNode[path.nodes[i].first].push_back(pos);
			fractionalIndex += 1;
		}
		if (fractionalIndex == missingFractionalBreakpoints.size()) break;
		nodeStartPos = nodeEndPos;
	}
	assert(fractionalIndex == missingFractionalBreakpoints.size());
	phmap::flat_hash_map<size_t, std::vector<std::pair<size_t, bool>>> newNodesPerNode;
	phmap::flat_hash_map<size_t, std::vector<double>> fractionalSizesPerNode;
	for (auto pair : fractionalCutsPerNode)
	{
		std::sort(pair.second.begin(), pair.second.end());
		for (size_t j = 0; j <= pair.second.size(); j++)
		{
			double cutStart = 0;
			double cutEnd = 1;
			if (j != 0) cutStart = pair.second[j-1];
			if (j < pair.second.size()) cutEnd = pair.second[j];
			fractionalSizesPerNode[pair.first].push_back(cutEnd-cutStart);
		}
		std::vector<std::pair<size_t, bool>> newNodes = cutNode(readLengths, readPaths, readsPerNode, pair.first, fractionalSizesPerNode[pair.first]);
		newNodesPerNode[pair.first] = newNodes;
	}
	Path result;
	result.startFractionalClip = path.startFractionalClip;
	result.endFractionalClip = path.endFractionalClip;
	for (size_t i = 0; i < path.nodes.size(); i++)
	{
		if (fractionalCutsPerNode.count(path.nodes[i].first) == 0)
		{
			result.nodes.push_back(path.nodes[i]);
			result.lengths.push_back(path.lengths[i]);
			continue;
		}
		std::vector<std::pair<size_t, bool>> newNodes = newNodesPerNode.at(path.nodes[i].first);
		std::vector<double> fractionalSizes = fractionalSizesPerNode.at(path.nodes[i].first);
		assert(newNodes.size() >= 2);
		assert(newNodes.size() == fractionalSizes.size());
		if (!path.nodes[i].second)
		{
			std::reverse(newNodes.begin(), newNodes.end());
			for (size_t j = 0; j < newNodes.size(); j++)
			{
				newNodes[j].second = !newNodes[j].second;
			}
			std::reverse(fractionalSizes.begin(), fractionalSizes.end());
		}
		size_t lengthUsed = 0;
		for (size_t j = 0; j < newNodes.size(); j++)
		{
			result.nodes.push_back(newNodes[j]);
			size_t lengthHere = path.lengths[i] * fractionalSizes[j];
			result.lengths.push_back(lengthHere);
			lengthUsed += lengthHere;
		}
		assert(lengthUsed <= path.lengths[i]);
		while (lengthUsed < path.lengths[i])
		{
			for (size_t j = 0; j < newNodes.size() && lengthUsed < path.lengths[i]; j++)
			{
				result.lengths[result.lengths.size()-newNodes.size()+j] += 1;
				lengthUsed += 1;
			}
		}
		assert(lengthUsed == path.lengths[i]);
	}
	assert(result.lengthIncludingClippedOutParts() == path.lengthIncludingClippedOutParts());
	return result;
}

std::vector<double> getMissingFractionalBreakpoints(const std::vector<double>& missingFrom, const std::vector<double>& includeFromThis)
{
	size_t missingIndex = 0;
	size_t includeIndex = 0;
	std::vector<double> result;
	while (missingIndex < missingFrom.size() || includeIndex < includeFromThis.size())
	{
		if (missingIndex == missingFrom.size())
		{
			result.push_back(includeFromThis[includeIndex]);
			includeIndex += 1;
			continue;
		}
		if (includeIndex == includeFromThis.size())
		{
			missingIndex += 1;
			continue;
		}
		if (missingFrom[missingIndex] > includeFromThis[includeIndex] + 0.01)
		{
			result.push_back(includeFromThis[includeIndex]);
			includeIndex += 1;
			continue;
		}
		if (includeFromThis[includeIndex] > missingFrom[missingIndex] + 0.01)
		{
			missingIndex += 1;
			continue;
		}
		assert(!(missingFrom[missingIndex] > includeFromThis[includeIndex] + 0.01));
		assert(!(includeFromThis[includeIndex] > missingFrom[missingIndex] + 0.01));
		missingIndex += 1;
		includeIndex += 1;
	}
	return result;
}

bool cutsDuplicateNodes(const Path& path, const std::vector<double>& fractionalCutPositions)
{
	phmap::flat_hash_map<size_t, size_t> nodeCount;
	for (size_t i = 0; i < path.nodes.size(); i++)
	{
		nodeCount[path.nodes[i].first] += 1;
	}
	size_t startPos = 0;
	size_t length = path.lengthIncludingClippedOutParts();
	size_t cutIndex = 0;
	for (size_t i = 1; i < fractionalCutPositions.size(); i++)
	{
		assert(fractionalCutPositions[i-1] <= fractionalCutPositions[i]);
	}
	for (size_t i = 0; i < path.nodes.size(); i++)
	{
		size_t endPos = startPos + path.lengths[i];
		double endFractionHere = (double)endPos / (double)length;
		while (cutIndex < fractionalCutPositions.size() && fractionalCutPositions[cutIndex] < endFractionHere)
		{
			if (nodeCount.at(path.nodes[i].first) != 1) return true;
			cutIndex += 1;
		}
		startPos = endPos;
	}
	assert(startPos == length);
	assert(cutIndex == fractionalCutPositions.size());
	return false;
}

std::pair<Path, Path> cutPaths(const std::vector<size_t>& readLengths, std::vector<std::vector<std::pair<std::pair<size_t, bool>, size_t>>>& readPaths, std::vector<std::vector<std::pair<size_t, size_t>>>& readsPerNode, const Path& leftPath, const Path& rightPath)
{
	std::vector<double> leftFractionalBreakpoints;
	std::vector<double> rightFractionalBreakpoints;
	size_t leftLength = leftPath.lengthIncludingClippedOutParts();
	size_t rightLength = rightPath.lengthIncludingClippedOutParts();
	std::vector<double> alignmentFractionalBreakpoints;
	double alignmentFractionalBreakpointLeftStart = leftPath.startFractionalClip * (double)leftPath.lengths[0] / (double)(leftLength);
	double alignmentFractionalBreakpointRightStart = rightPath.startFractionalClip * (double)rightPath.lengths[0] / (double)(rightLength);
	double alignmentFractionalBreakpointLeftEnd = 1.0 - leftPath.endFractionalClip * (double)leftPath.lengths.back() / (double)(leftLength);
	double alignmentFractionalBreakpointRightEnd = 1.0 - rightPath.endFractionalClip * (double)rightPath.lengths.back() / (double)(rightLength);
	size_t nodeStartPos = 0;
	assert(leftPath.lengths.back() > 0);
	assert(rightPath.lengths.back() > 0);
	for (size_t i = 0; i < leftPath.nodes.size(); i++)
	{
		if (i > 0) leftFractionalBreakpoints.emplace_back((double)nodeStartPos / (double)leftLength);
		nodeStartPos += leftPath.lengths[i];
	}
	assert(nodeStartPos == leftLength);
	nodeStartPos = 0;
	for (size_t i = 0; i < rightPath.nodes.size(); i++)
	{
		if (i > 0) rightFractionalBreakpoints.emplace_back((double)nodeStartPos / (double)rightLength);
		nodeStartPos += rightPath.lengths[i];
	}
	assert(nodeStartPos == rightLength);
	for (size_t i = 1; i < leftFractionalBreakpoints.size(); i++)
	{
		assert(leftFractionalBreakpoints[i-1] <= leftFractionalBreakpoints[i]);
		assert(leftFractionalBreakpoints[i] < 1);
	}
	for (size_t i = 1; i < rightFractionalBreakpoints.size(); i++)
	{
		assert(rightFractionalBreakpoints[i-1] <= rightFractionalBreakpoints[i]);
		assert(rightFractionalBreakpoints[i] < 1);
	}
	std::vector<double> fractionalBreakpointsMissingInLeft = getMissingFractionalBreakpoints(leftFractionalBreakpoints, rightFractionalBreakpoints);
	std::vector<double> fractionalBreakpointsMissingInRight = getMissingFractionalBreakpoints(rightFractionalBreakpoints, leftFractionalBreakpoints);
	bool addedToStart = false;
	bool addedToEnd = false;
	if (alignmentFractionalBreakpointLeftStart != 0 || alignmentFractionalBreakpointRightStart != 0)
	{
		if (alignmentFractionalBreakpointLeftStart + 0.01 < alignmentFractionalBreakpointRightStart && alignmentFractionalBreakpointRightStart + 0.01 < alignmentFractionalBreakpointLeftStart)
		{
			fractionalBreakpointsMissingInLeft.emplace_back((alignmentFractionalBreakpointRightStart + alignmentFractionalBreakpointLeftStart)/2);
			fractionalBreakpointsMissingInRight.emplace_back((alignmentFractionalBreakpointRightStart + alignmentFractionalBreakpointLeftStart)/2);
		}
		else
		{
			if (alignmentFractionalBreakpointLeftStart > 0)
			{
				fractionalBreakpointsMissingInLeft.emplace_back(alignmentFractionalBreakpointLeftStart);
				fractionalBreakpointsMissingInRight.emplace_back(alignmentFractionalBreakpointLeftStart);
			}
			if (alignmentFractionalBreakpointRightStart > 0)
			{
				fractionalBreakpointsMissingInLeft.emplace_back(alignmentFractionalBreakpointRightStart);
				fractionalBreakpointsMissingInRight.emplace_back(alignmentFractionalBreakpointRightStart);
			}
		}
		addedToStart = true;
	}
	if (alignmentFractionalBreakpointLeftEnd != 1 || alignmentFractionalBreakpointRightEnd != 1)
	{
		if (alignmentFractionalBreakpointLeftEnd + 0.01 < alignmentFractionalBreakpointRightEnd && alignmentFractionalBreakpointRightEnd + 0.01 < alignmentFractionalBreakpointLeftEnd)
		{
			fractionalBreakpointsMissingInLeft.emplace_back((alignmentFractionalBreakpointRightEnd + alignmentFractionalBreakpointLeftEnd)/2);
			fractionalBreakpointsMissingInRight.emplace_back((alignmentFractionalBreakpointRightEnd + alignmentFractionalBreakpointLeftEnd)/2);
		}
		else
		{
			if (alignmentFractionalBreakpointLeftEnd != 1)
			{
				fractionalBreakpointsMissingInLeft.emplace_back(alignmentFractionalBreakpointLeftEnd);
				fractionalBreakpointsMissingInRight.emplace_back(alignmentFractionalBreakpointLeftEnd);
			}
			if (alignmentFractionalBreakpointRightEnd != 1)
			{
				fractionalBreakpointsMissingInLeft.emplace_back(alignmentFractionalBreakpointRightEnd);
				fractionalBreakpointsMissingInRight.emplace_back(alignmentFractionalBreakpointRightEnd);
			}
		}
		addedToEnd = true;
	}
	std::sort(fractionalBreakpointsMissingInLeft.begin(), fractionalBreakpointsMissingInLeft.end());
	std::sort(fractionalBreakpointsMissingInRight.begin(), fractionalBreakpointsMissingInRight.end());
	std::cerr << "left breakpoints";
	for (auto pos : leftFractionalBreakpoints) std::cerr << " " << pos;
	std::cerr << std::endl;
	std::cerr << "right breakpoints";
	for (auto pos : rightFractionalBreakpoints) std::cerr << " " << pos;
	std::cerr << std::endl;
	std::cerr << "missing left breakpoints";
	for (auto pos : fractionalBreakpointsMissingInLeft) std::cerr << " " << pos;
	std::cerr << std::endl;
	std::cerr << "missing right breakpoints";
	for (auto pos : fractionalBreakpointsMissingInRight) std::cerr << " " << pos;
	std::cerr << std::endl;
	assert(leftFractionalBreakpoints.size() + fractionalBreakpointsMissingInLeft.size() == rightFractionalBreakpoints.size() + fractionalBreakpointsMissingInRight.size());
	if (cutsDuplicateNodes(leftPath, fractionalBreakpointsMissingInLeft))
	{
		std::cerr << "cuts duplicate nodes" << std::endl;
		return std::make_pair(Path{}, Path{});
	}
	if (cutsDuplicateNodes(rightPath, fractionalBreakpointsMissingInRight))
	{
		std::cerr << "cuts duplicate nodes" << std::endl;
		return std::make_pair(Path{}, Path{});
	}
	auto newLeftPath = cutPath(readLengths, readPaths, readsPerNode, leftPath, fractionalBreakpointsMissingInLeft);
	auto newRightPath = cutPath(readLengths, readPaths, readsPerNode, rightPath, fractionalBreakpointsMissingInRight);
	std::cerr << "cut path ";
	for (auto node : leftPath.nodes)
	{
		std::cerr << (node.second ? ">" : "<") << node.first;
	}
	std::cerr << " into ";
	for (auto node : newLeftPath.nodes)
	{
		std::cerr << (node.second ? ">" : "<") << node.first;
	}
	std::cerr << std::endl;
	std::cerr << "cut path ";
	for (auto node : rightPath.nodes)
	{
		std::cerr << (node.second ? ">" : "<") << node.first;
	}
	std::cerr << " into ";
	for (auto node : newRightPath.nodes)
	{
		std::cerr << (node.second ? ">" : "<") << node.first;
	}
	std::cerr << std::endl;
	assert(newLeftPath.nodes.size() == newRightPath.nodes.size());
	assert(newLeftPath.nodes.size() > addedToStart + addedToEnd);
	if (addedToStart)
	{
		newLeftPath.nodes.erase(newLeftPath.nodes.begin());
		newLeftPath.lengths.erase(newLeftPath.lengths.begin());
		newLeftPath.startFractionalClip = 0;
		newRightPath.nodes.erase(newRightPath.nodes.begin());
		newRightPath.lengths.erase(newRightPath.lengths.begin());
		newRightPath.startFractionalClip = 0;
	}
	if (addedToEnd)
	{
		newLeftPath.nodes.pop_back();
		newLeftPath.lengths.pop_back();
		newLeftPath.endFractionalClip = 0;
		newRightPath.nodes.pop_back();
		newRightPath.lengths.pop_back();
		newRightPath.endFractionalClip = 0;
	}
	assert(newLeftPath.nodes.size() == newRightPath.nodes.size());
	return std::make_pair(newLeftPath, newRightPath);
}

void mergePaths(const std::vector<size_t>& readLengths, std::vector<std::vector<std::pair<std::pair<size_t, bool>, size_t>>>& readPaths, std::vector<std::vector<std::pair<size_t, size_t>>>& readsPerNode, const Path& leftPath, const Path& rightPath)
{
	std::cerr << "merge nodes ";
	for (auto node : leftPath.nodes)
	{
		std::cerr << (node.second ? ">": "<") << node.first;
	}
	std::cerr << " and ";
	for (auto node : rightPath.nodes)
	{
		std::cerr << (node.second ? ">": "<") << node.first;
	}
	std::cerr << std::endl;
	assert(leftPath.nodes.size() == rightPath.nodes.size());
	for (size_t i = 0; i < leftPath.nodes.size(); i++)
	{
		bool sameOrientation = (leftPath.nodes[i].second == rightPath.nodes[i].second);
		for (auto pos : readsPerNode[rightPath.nodes[i].first])
		{
			assert(readPaths[pos.first][pos.second].first.first == rightPath.nodes[i].first);
			assert(readPaths[pos.first][pos.second].first.first != leftPath.nodes[i].first);
			readPaths[pos.first][pos.second].first.first = leftPath.nodes[i].first;
			if (!sameOrientation) readPaths[pos.first][pos.second].first.second = !readPaths[pos.first][pos.second].first.second;
			readsPerNode[leftPath.nodes[i].first].emplace_back(pos);
		}
		readsPerNode[rightPath.nodes[i].first].clear();
	}
}

void cutAndMergePathPairs(const std::vector<size_t>& readLengths, std::vector<std::vector<std::pair<std::pair<size_t, bool>, size_t>>>& readPaths, std::vector<std::vector<std::pair<size_t, size_t>>>& readsPerNode, const Path& leftPath, const Path& rightPath)
{
	Path newLeftPath, newRightPath;
	std::tie(newLeftPath, newRightPath) = cutPaths(readLengths, readPaths, readsPerNode, leftPath, rightPath);
	if (newLeftPath.nodes.size() == 0 && newRightPath.nodes.size() == 0) return; //something bad happened in cut, don't do merge
	mergePaths(readLengths, readPaths, readsPerNode, newLeftPath, newRightPath);
}

void unitigify(const std::vector<size_t>& readLengths, std::vector<std::vector<std::pair<std::pair<size_t, bool>, size_t>>>& readPaths, std::vector<std::vector<std::pair<size_t, size_t>>>& readsPerNode, const std::vector<std::pair<size_t, bool>>& unitigifyPath)
{
	size_t newNode = readsPerNode.size();
	std::cerr << "unitigify ";
	for (auto node : unitigifyPath)
	{
		std::cerr << (node.second ? ">" : "<") << node.first;
	}
	std::cerr << " into " << newNode << std::endl;
	readsPerNode.emplace_back();
	phmap::flat_hash_set<size_t> relevantReads;
	for (auto node : unitigifyPath)
	{
		for (auto pos : readsPerNode[node.first])
		{
			relevantReads.insert(pos.first);
		}
	}
	for (const size_t read : relevantReads)
	{
		for (size_t i = 0; i < readPaths[read].size(); i++)
		{
			bool found = false;
			for (size_t j = 0; j < readsPerNode[readPaths[read][i].first.first].size(); j++)
			{
				if (readsPerNode[readPaths[read][i].first.first][j].first != read) continue;
				if (readsPerNode[readPaths[read][i].first.first][j].second != i) continue;
				std::swap(readsPerNode[readPaths[read][i].first.first][j], readsPerNode[readPaths[read][i].first.first].back());
				readsPerNode[readPaths[read][i].first.first].pop_back();
				assert(!found);
				found = true;
			}
			assert(found);
		}
		std::vector<std::pair<std::pair<size_t, bool>, size_t>> newPath;
		for (size_t i = 0; i < readPaths[read].size(); i++)
		{
			if (readPaths[read][i].first == unitigifyPath[0])
			{
				assert(i+unitigifyPath.size() <= readPaths[read].size());
				size_t totalLength = 0;
				for (size_t j = 0; j < unitigifyPath.size(); j++)
				{
					assert(readPaths[read][i+j].first == unitigifyPath[j]);
					totalLength += readPaths[read][i+j].second;
				}
				newPath.emplace_back(std::make_pair(std::make_pair(newNode, true), totalLength));
				i += unitigifyPath.size()-1;
			}
			else if (readPaths[read][i].first == reverse(unitigifyPath.back()))
			{
				assert(i+unitigifyPath.size() <= readPaths[read].size());
				size_t totalLength = 0;
				for (size_t j = 0; j < unitigifyPath.size(); j++)
				{
					assert(readPaths[read][i+j].first == reverse(unitigifyPath[unitigifyPath.size()-1-j]));
					totalLength += readPaths[read][i+j].second;
				}
				newPath.emplace_back(std::make_pair(std::make_pair(newNode, false), totalLength));
				i += unitigifyPath.size()-1;
			}
			else
			{
				newPath.emplace_back(readPaths[read][i]);
			}
		}
		readPaths[read] = newPath;
		for (size_t i = 0; i < readPaths[read].size(); i++)
		{
			readsPerNode[readPaths[read][i].first.first].emplace_back(read, i);
		}
	}
}

void tryUnitigifyPath(const std::vector<size_t>& readLengths, std::vector<std::vector<std::pair<std::pair<size_t, bool>, size_t>>>& readPaths, std::vector<std::vector<std::pair<size_t, size_t>>>& readsPerNode, const size_t pathIndex, const bool pathFw, const size_t pathStart, const size_t pathEnd)
{
	auto relevantPath = getPath(readLengths, readPaths, pathIndex, pathFw, pathStart, pathEnd);
	std::cerr << "try unitigify ";
	for (auto node : relevantPath.nodes) std::cerr << (node.second ? ">" : "<") << node.first;
	std::cerr << std::endl;
	size_t unitigStart = 0;
	for (size_t i = 1; i < relevantPath.nodes.size(); i++)
	{
		phmap::flat_hash_set<std::pair<size_t, bool>> prevEdges;
		for (auto pos : readsPerNode[relevantPath.nodes[i-1].first])
		{
			if (readPaths[pos.first][pos.second].first.second == relevantPath.nodes[i-1].second && pos.second+1 < readPaths[pos.first].size())
			{
				prevEdges.insert(readPaths[pos.first][pos.second+1].first);
			}
			else if (readPaths[pos.first][pos.second].first.second != relevantPath.nodes[i-1].second && pos.second != 0)
			{
				prevEdges.insert(reverse(readPaths[pos.first][pos.second-1].first));
			}
		}
		phmap::flat_hash_set<std::pair<size_t, bool>> nextEdges;
		for (auto pos : readsPerNode[relevantPath.nodes[i].first])
		{
			if (readPaths[pos.first][pos.second].first.second == relevantPath.nodes[i].second && pos.second != 0)
			{
				nextEdges.insert(reverse(readPaths[pos.first][pos.second-1].first));
			}
			else if (readPaths[pos.first][pos.second].first.second != relevantPath.nodes[i].second && pos.second+1 < readPaths[pos.first].size())
			{
				nextEdges.insert(readPaths[pos.first][pos.second+1].first);
			}
		}
		if (!(prevEdges.count(relevantPath.nodes[i]) == 1))
		{
			std::cerr << i << std::endl;
		}
		if (!(nextEdges.count(reverse(relevantPath.nodes[i-1])) == 1))
		{
			std::cerr << i << std::endl;
		}
		assert(prevEdges.count(relevantPath.nodes[i]) == 1);
		assert(nextEdges.count(reverse(relevantPath.nodes[i-1])) == 1);
		if (prevEdges.size() == 1 && nextEdges.size() == 1) continue;
		if (i > unitigStart+1)
		{
			std::vector<std::pair<size_t, bool>> unitigifyPath;
			for (size_t j = unitigStart; j < i; j++)
			{
				unitigifyPath.emplace_back(relevantPath.nodes[j]);
			}
			assert(unitigifyPath.size() >= 2);
			unitigify(readLengths, readPaths, readsPerNode, unitigifyPath);
			// recursive just in case the same unitigification appears twice in the same path
			tryUnitigifyPath(readLengths, readPaths, readsPerNode, pathIndex, pathFw, pathStart, pathEnd);
			return;
		}
		unitigStart = i;
	}
	if (unitigStart != relevantPath.nodes.size()-1)
	{
		std::vector<std::pair<size_t, bool>> unitigifyPath;
		for (size_t j = unitigStart; j < relevantPath.nodes.size(); j++)
		{
			unitigifyPath.emplace_back(relevantPath.nodes[j]);
		}
		assert(unitigifyPath.size() >= 2);
		unitigify(readLengths, readPaths, readsPerNode, unitigifyPath);
	}
}

void trimPartialPath(Path& path)
{
	assert(path.lengths.size() >= 1);
	while (path.lengths.back() == 0)
	{
		assert(path.endFractionalClip == 0);
		path.lengths.pop_back();
		path.nodes.pop_back();
		assert(path.lengths.size() >= 1);
	}
	while (path.lengths[0] == 0)
	{
		assert(path.startFractionalClip == 0);
		path.lengths.erase(path.lengths.begin());
		path.nodes.erase(path.nodes.begin());
		assert(path.lengths.size() >= 1);
	}
}

void joinPaths(const std::vector<size_t>& readLengths, std::vector<std::vector<std::pair<std::pair<size_t, bool>, size_t>>>& readPaths, std::vector<std::vector<std::pair<size_t, size_t>>>& readsPerNode, const size_t leftPathIndex, const bool leftPathFw, const size_t leftPathStart, const size_t leftPathEnd, const size_t rightPathIndex, const bool rightPathFw, const size_t rightPathStart, const size_t rightPathEnd)
{
	std::cerr << "join " << leftPathIndex << " " << leftPathFw << " " << leftPathStart << " " << leftPathEnd << " to " << rightPathIndex << " " << rightPathFw << " " << rightPathStart << " " << rightPathEnd << std::endl;
	Path leftPath = getPath(readLengths, readPaths, leftPathIndex, leftPathFw, leftPathStart, leftPathEnd);
	Path rightPath = getPath(readLengths, readPaths, rightPathIndex, rightPathFw, rightPathStart, rightPathEnd);
	std::cerr << "joinfull ";
	for (auto node : leftPath.nodes)
	{
		std::cerr << (node.second ? ">" : "<") << node.first;
	}
	std::cerr << " to ";
	for (auto node : rightPath.nodes)
	{
		std::cerr << (node.second ? ">" : "<") << node.first;
	}
	std::cerr << std::endl;
	auto partialPaths = getPartialPathsAtNonequalSites(leftPath, rightPath);
	for (size_t i = 0; i < partialPaths.size(); i++)
	{
		if (partialPaths[i].first.lengthIncludingClippedOutParts() == 0 || partialPaths[i].second.lengthIncludingClippedOutParts() == 0) continue;
		trimPartialPath(partialPaths[i].first);
		trimPartialPath(partialPaths[i].second);
		std::cerr << "joinpart ";
		for (auto node : partialPaths[i].first.nodes)
		{
			std::cerr << (node.second ? ">" : "<") << node.first;
		}
		std::cerr << " to ";
		for (auto node : partialPaths[i].second.nodes)
		{
			std::cerr << (node.second ? ">" : "<") << node.first;
		}
		std::cerr << std::endl;
		cutAndMergePathPairs(readLengths, readPaths, readsPerNode, partialPaths[i].first, partialPaths[i].second);
	}
	tryUnitigifyPath(readLengths, readPaths, readsPerNode, leftPathIndex, leftPathFw, 0, readLengths[leftPathIndex]);
}

void makeGraph(const std::vector<size_t>& readLengths, const std::vector<std::tuple<size_t, size_t, size_t, size_t, size_t, size_t, bool>>& matches, size_t minCoverage, const std::string& outfile)
{
	std::vector<std::vector<std::pair<std::pair<size_t, bool>, size_t>>> readPaths;
	std::vector<std::vector<std::pair<size_t, size_t>>> readsPerNode;
	readPaths.resize(readLengths.size());
	readsPerNode.resize(readLengths.size()*3);
	for (size_t i = 0; i < readLengths.size(); i++)
	{
		readPaths[i].emplace_back(std::make_pair(i*3, true), 0);
		readPaths[i].emplace_back(std::make_pair(i*3+1, true), readLengths[i]);
		readPaths[i].emplace_back(std::make_pair(i*3+2, true), 0);
		readsPerNode[i*3].emplace_back(i, 0);
		readsPerNode[i*3+1].emplace_back(i, 1);
		readsPerNode[i*3+2].emplace_back(i, 2);
	}
	for (size_t i = 0; i < matches.size(); i++)
	{
		joinPaths(readLengths, readPaths, readsPerNode, std::get<0>(matches[i]), true, std::get<1>(matches[i]), std::get<2>(matches[i]), std::get<3>(matches[i]), std::get<6>(matches[i]), std::get<4>(matches[i]), std::get<5>(matches[i]));
		for (size_t j = 0; j < readsPerNode.size(); j++)
		{
			for (auto pos : readsPerNode[j])
			{
				if (!(readPaths[pos.first][pos.second].first.first == j)) std::cerr << j << std::endl;
				assert(readPaths[pos.first][pos.second].first.first == j);
			}
		}
	}
	std::ofstream graphOut { outfile };
	for (size_t i = 0; i < readsPerNode.size(); i++)
	{
		if (readsPerNode[i].size() < minCoverage) continue;
		size_t totalLength = 0;
		for (auto pos : readsPerNode[i])
		{
			assert(readPaths[pos.first][pos.second].first.first == i);
			totalLength += readPaths[pos.first][pos.second].second;
		}
		size_t nodeLength = totalLength / readsPerNode[i].size();
		size_t coverage = readsPerNode[i].size();
		graphOut << "S\t" << i << "\t*\tLN:i:" << nodeLength << "\tll:f:" << coverage << "\tFC:i:" << coverage*nodeLength << std::endl;
		phmap::flat_hash_map<std::pair<size_t, bool>, size_t> fwEdgeCoverages;
		phmap::flat_hash_map<std::pair<size_t, bool>, size_t> bwEdgeCoverages;
		for (auto pos : readsPerNode[i])
		{
			assert(readPaths[pos.first][pos.second].first.first == i);
			if (pos.second > 0)
			{
				auto prev = reverse(readPaths[pos.first][pos.second-1].first);
				if (readPaths[pos.first][pos.second].first.second)
				{
					bwEdgeCoverages[prev] += 1;
				}
				else
				{
					fwEdgeCoverages[prev] += 1;
				}
			}
			if (pos.second < readPaths[pos.first].size()-1)
			{
				auto next = (readPaths[pos.first][pos.second+1].first);
				if (readPaths[pos.first][pos.second].first.second)
				{
					fwEdgeCoverages[next] += 1;
				}
				else
				{
					bwEdgeCoverages[next] += 1;
				}
			}
		}
		for (auto pair : fwEdgeCoverages)
		{
			if (pair.second < minCoverage) continue;
			graphOut << "L\t" << i << "\t+\t" << pair.first.first << "\t" << (pair.first.second ? "+" : "-") << "\t0M\tec:i:" << pair.second << std::endl;
		}
		for (auto pair : bwEdgeCoverages)
		{
			if (pair.second < minCoverage) continue;
			graphOut << "L\t" << i << "\t-\t" << pair.first.first << "\t" << (pair.first.second ? "+" : "-") << "\t0M\tec:i:" << pair.second << std::endl;
		}
	}
}

int main(int argc, char** argv)
{
	size_t numThreads = std::stoi(argv[1]);
	size_t k = std::stoull(argv[2]);
	size_t numWindows = std::stoull(argv[3]);
	size_t windowSize = std::stoull(argv[4]);
	size_t minAlignmentLength = std::stoull(argv[5]);
	std::vector<std::string> readFiles;
	for (size_t i = 6; i < argc; i++)
	{
		readFiles.emplace_back(argv[i]);
	}
	MatchIndex matchIndex { k, numWindows, windowSize };
	ReadStorage storage;
	for (auto file : readFiles)
	{
		std::mutex indexMutex;
		storage.iterateReadsFromFile(file, numThreads, false, [&matchIndex, &indexMutex](size_t readName, const std::string& sequence)
		{
			matchIndex.addMatchesFromRead(readName, indexMutex, sequence);
		});
	}
	size_t numWindowChunks = matchIndex.numWindowChunks();
	size_t numUniqueChunks = matchIndex.numUniqueChunks();
	matchIndex.clearConstructionVariablesAndCompact();
	const std::vector<size_t>& readLengths = storage.getRawReadLengths();
	std::mutex printMutex;
	std::vector<std::tuple<size_t, size_t, size_t, size_t, size_t, size_t, bool>> matches;
	auto result = matchIndex.iterateMatchChains(numThreads, storage.getRawReadLengths(), [&printMutex, &matches, minAlignmentLength, &readLengths](const size_t left, const size_t leftstart, const size_t leftend, const bool leftFw, const size_t right, const size_t rightstart, const size_t rightend, const bool rightFw)
	{
		if (leftend+1-leftstart < minAlignmentLength) return;
		if (rightend+1-rightstart < minAlignmentLength) return;
		std::lock_guard<std::mutex> lock { printMutex };
		assert(leftFw);
		matches.emplace_back(left, leftstart, leftend+1, right, rightstart, rightend+1, rightFw); // param coordinates are inclusive, switch to exclusive
		//matches.emplace_back(left, readLengths[left]-leftend-1, readLengths[left]-leftstart, right, readLengths[right]-rightend-1, readLengths[right]-rightstart, !rightFw); // param coordinates are inclusive, switch to exclusive
	});
	std::sort(matches.begin(), matches.end(), [](auto left, auto right) { return std::get<2>(left)-std::get<1>(left) + std::get<5>(left)-std::get<4>(left) > std::get<2>(right)-std::get<1>(right) + std::get<5>(right)-std::get<4>(right); });
	makeGraph(readLengths, matches, 2, "graph.gfa");
	std::cerr << result.numberReads << " reads" << std::endl;
	std::cerr << numWindowChunks << " distinct windowchunks" << std::endl;
	std::cerr << result.totalReadChunkMatches << " read-windowchunk matches (except unique)" << std::endl;
	std::cerr << numUniqueChunks << " windowchunks have only one read" << std::endl;
	std::cerr << result.readsWithMatch << " reads with a match" << std::endl;
	std::cerr << result.readPairMatches << " read-read matches" << std::endl;
	std::cerr << result.readChainMatches << " chain matches" << std::endl;
	std::cerr << result.totalMatches << " window matches" << std::endl;
	std::cerr << result.maxPerChunk << " max windowchunk size" << std::endl;
}
