//
// REViewer
// Copyright 2020 Illumina, Inc.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#include "app/Projection.hh"

#include <algorithm>
#include <set>
#include <string>

#include "graphalign/Operation.hh"

using boost::optional;
using graphtools::Alignment;
using graphtools::GraphAlignment;
using graphtools::NodeId;
using graphtools::Operation;
using graphtools::OperationType;
using graphtools::Path;
using std::string;
using std::unordered_map;
using std::vector;

int score(const GraphAlignment& alignment, int matchScore, int mismatchScore, int gapScore)
{
    int score = 0;
    for (const auto& nodeAlign : alignment.alignments())
    {
        for (const Operation& operation : nodeAlign)
        {
            switch (operation.type())
            {
            case OperationType::kMatch:
                score += matchScore * operation.referenceLength();
                break;
            case OperationType::kMismatch:
                score += mismatchScore * operation.referenceLength();
                break;
            case OperationType::kInsertionToRef:
                score += gapScore * operation.queryLength();
                break;
            case OperationType::kDeletionFromRef:
                score += gapScore * operation.referenceLength();
                break;
            default:
                break;
            }
        }
    }

    return score;
}

bool checkCommonElement(const vector<NodeId>& nodesA, int startIndexA, const vector<NodeId>& nodesB, int startIndexB)
{
    assert(startIndexA < nodesA.size() && startIndexB < nodesB.size());
    while (startIndexA != nodesA.size() && startIndexB != nodesB.size())
    {
        if (nodesA[startIndexA] < nodesB[startIndexB])
        {
            ++startIndexA;
        }
        else if (nodesB[startIndexB] < nodesA[startIndexA])
        {
            ++startIndexB;
        }
        else
        {
            assert(nodesA[startIndexA] == nodesB[startIndexB]);
            return true;
        }
    }

    return false;
}

vector<int> getStartIndexes(const Path& projPath, int initialStartIndex, const GraphAlignment& align)
{
    // Check if IRR
    std::set<NodeId> alignNodes(align.path().nodeIds().begin(), align.path().nodeIds().end());
    if (alignNodes.size() != 1)
    {
        return { initialStartIndex };
    }

    NodeId repeatNode = *alignNodes.begin();
    if (!projPath.graphRawPtr()->hasEdge(repeatNode, repeatNode))
    {
        return { initialStartIndex };
    }

    int numMotifsInPath = std::count(projPath.begin(), projPath.end(), repeatNode);
    int numMotifsInAlign = std::count(align.path().nodeIds().begin(), align.path().nodeIds().end(), repeatNode);

    if (numMotifsInPath <= numMotifsInAlign)
    {
        return { initialStartIndex };
    }

    auto repeatStartIndexIt = std::find(projPath.begin(), projPath.end(), repeatNode);
    assert(repeatStartIndexIt != projPath.end());
    int repeatStartIndex = repeatStartIndexIt - projPath.begin();

    vector<int> startIndexes;
    for (int index = 0; index != numMotifsInPath - numMotifsInAlign + 1; ++index)
    {
        startIndexes.push_back(repeatStartIndex + index);
    }

    return startIndexes;
}

struct AlignProj
{
    AlignProj(vector<int> startIndexes, GraphAlignPtr align)
        : startIndexes(std::move(startIndexes))
        , align(std::move(align))
    {
    }

    vector<int> startIndexes;
    GraphAlignPtr align;
};

static optional<AlignProj> project(const GraphAlign& align, const Path& projPath)
{
    // Find the first common node
    int alignNodeIndex = 0;
    int projNodeIndex = 0;

    const auto& alignNodes = align.path().nodeIds();
    const auto& projPathNodes = projPath.nodeIds();
    while (alignNodes[alignNodeIndex] != projPathNodes[projNodeIndex])
    {
        if (alignNodes[alignNodeIndex] < projPathNodes[projNodeIndex])
        {
            ++alignNodeIndex;
            if (alignNodeIndex == alignNodes.size())
            {
                return boost::none;
            }
        }
        else if (projPathNodes[projNodeIndex] < alignNodes[alignNodeIndex])
        {
            ++projNodeIndex;
            if (projNodeIndex == projPathNodes.size())
            {
                return boost::none;
            }
        }
        else
        {
            assert(false);
        }
    }

    // Line up the nodes
    graphtools::NodeId commonNode = alignNodes[alignNodeIndex];
    const int alignCommonNodeCount = std::count(alignNodes.begin(), alignNodes.end(), commonNode);
    const int projCommonNodeCount = std::count(projPathNodes.begin(), projPathNodes.end(), commonNode);

    if (alignCommonNodeCount < projCommonNodeCount)
    {
        projNodeIndex += projCommonNodeCount - alignCommonNodeCount;
    }
    else
    {
        alignNodeIndex += alignCommonNodeCount - projCommonNodeCount;
    }

    // Project starting from the common node
    vector<graphtools::NodeId> projNodes;
    vector<Alignment> projAligns;

    int projStartNodeIndex = projNodeIndex;
    int startingAlignIndex = alignNodeIndex;
    int endingAlignIndex = alignNodeIndex;

    while (alignNodeIndex != alignNodes.size())
    {
        if (!checkCommonElement(alignNodes, alignNodeIndex, projPathNodes, projNodeIndex))
        {
            break;
        }

        NodeId alignNode = alignNodes[alignNodeIndex];
        NodeId targetNode = projPathNodes[projNodeIndex];

        if (alignNode == targetNode)
        {
            projAligns.push_back(align.alignments()[alignNodeIndex]);
            projNodes.push_back(targetNode);
            endingAlignIndex = alignNodeIndex;

            ++alignNodeIndex;
            ++projNodeIndex;
            if (alignNodeIndex == alignNodes.size() || projNodeIndex == projPathNodes.size())
            {
                break;
            }
        }
        else if (alignNode < targetNode) // Add an insertion
        {
            auto queryLen = align.alignments()[alignNodeIndex].queryLength();
            auto refStart = projAligns.back().referenceStart();
            auto operations = projAligns.back().operations();
            operations.emplace_back(OperationType::kInsertionToRef, queryLen);
            projAligns.back() = Alignment(refStart, operations);

            ++alignNodeIndex;
            if (alignNodeIndex == alignNodes.size())
            {
                break;
            }
        }
        else if (targetNode < alignNode) // Add a deletion
        {
            auto refLen = align.path().graphRawPtr()->nodeSeq(targetNode).length();
            projAligns.emplace_back(0, std::to_string(refLen) + "D");
            projNodes.push_back(targetNode);
            ++projNodeIndex;

            if (projNodeIndex == projPathNodes.size())
            {
                break;
            }
        }
        else
        {
            assert(false);
        }
    }

    endingAlignIndex = alignNodeIndex - 1;

    // Add left soft clip if needed
    int leftSoftclipLen = 0;
    for (int alignNodeIndex = 0; alignNodeIndex != startingAlignIndex; ++alignNodeIndex)
    {
        const auto& alignNode = align.alignments()[alignNodeIndex];
        leftSoftclipLen += alignNode.queryLength();
    }

    if (leftSoftclipLen)
    {
        int refStart = projAligns.front().referenceStart();
        auto operations = projAligns.front().operations();
        operations.emplace_front(OperationType::kSoftclip, leftSoftclipLen);
        projAligns.front() = Alignment(refStart, operations);
    }

    // Add right soft clip if needed
    int rightSoftclipLen = 0;
    for (int alignNodeIndex = endingAlignIndex + 1; alignNodeIndex != alignNodes.size(); ++alignNodeIndex)
    {
        const auto& nodeAlign = align.alignments()[alignNodeIndex];
        rightSoftclipLen += nodeAlign.queryLength();
    }

    if (rightSoftclipLen)
    {
        auto refStart = projAligns.back().referenceStart();
        auto operations = projAligns.back().operations();
        operations.emplace_back(OperationType::kSoftclip, rightSoftclipLen);
        projAligns.back() = Alignment(refStart, operations);
    }

    // Initialize projected alignment
    auto projPathStart = projAligns.front().referenceStart();
    auto projPathEnd = projAligns.back().referenceStart() + projAligns.back().referenceLength();
    Path projAlignPath(projPath.graphRawPtr(), projPathStart, projNodes, projPathEnd);
    GraphAlignPtr projAlign(new GraphAlignment(projAlignPath, projAligns));

    auto startIndexes = getStartIndexes(projPath, projStartNodeIndex, *projAlign);
    return AlignProj(startIndexes, std::move(projAlign));
}

vector<ReadPathAlign> project(const GraphAlign& align, int pathIndex, const Path& path)
{
    vector<ReadPathAlign> pathAligns;
    optional<AlignProj> alignProj = project(align, path);
    if (alignProj)
    {
        for (auto startIndex : alignProj->startIndexes)
        {
            pathAligns.emplace_back(path, pathIndex, startIndex, alignProj->align);
        }
    }

    return pathAligns;
}

PairPathAlignById project(const vector<Path>& genotypePaths, const PairGraphAlignById& pairGraphAlignById)
{
    PairPathAlignById pairPathAlignById;
    for (const auto& idAndPairGraphAlign : pairGraphAlignById)
    {
        PairPathAlign pathAlign;
        int bestPairScore = std::numeric_limits<int>::lowest();
        for (int pathIndex = 0; pathIndex != genotypePaths.size(); ++pathIndex)
        {
            const auto& path = genotypePaths[pathIndex];
            vector<ReadPathAlign> readPathAligns = project(idAndPairGraphAlign.second.readAlign, pathIndex, path);
            vector<ReadPathAlign> matePathAligns = project(idAndPairGraphAlign.second.mateAlign, pathIndex, path);
            if (readPathAligns.empty() || matePathAligns.empty())
            {
                continue;
            }

            const int pairScore = score(*readPathAligns.front().align) + score(*matePathAligns.front().align);
            if (pairScore > bestPairScore)
            {
                bestPairScore = pairScore;
                pathAlign.readAligns.clear();
                pathAlign.mateAligns.clear();
            }

            if (pairScore == bestPairScore)
            {
                pathAlign.readAligns.insert(pathAlign.readAligns.end(), readPathAligns.begin(), readPathAligns.end());
                pathAlign.mateAligns.insert(pathAlign.mateAligns.end(), matePathAligns.begin(), matePathAligns.end());
            }
        }

        assert(!pathAlign.readAligns.empty() && !pathAlign.mateAligns.empty());
        pairPathAlignById.emplace(idAndPairGraphAlign.first, pathAlign);
    }

    return pairPathAlignById;
}

/*
FragOriginsById getFragOrigins(const vector<Path>& hapPaths, const FragPathAlignsById& fragPathAlignsById)
{


    FragOriginsById fragOriginsById;


    assert(!fragOriginsById.empty());
    return fragOriginsById;
} */

ReadPathAlign::ReadPathAlign(const Path& hapPath, int pathIndex, int startIndexOnPath, GraphAlignPtr align)
    : pathIndex(pathIndex)
    , startIndexOnPath(startIndexOnPath)
    , align(std::move(align))
{
    int prefixLen = 0;
    for (int nodeIndex = 0; nodeIndex != startIndexOnPath; ++nodeIndex)
    {
        NodeId node = hapPath.nodeIds()[nodeIndex];
        prefixLen += hapPath.graphRawPtr()->nodeSeq(node).length();
    }

    begin = prefixLen + this->align->path().startPosition();
    end = begin + this->align->referenceLength();
}
