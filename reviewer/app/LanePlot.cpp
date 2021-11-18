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

#include "app/LanePlot.hh"

#include <tuple>

#include <boost/optional.hpp>

#include "graphalign/GraphAlignmentOperations.hh"
#include "graphalign/LinearAlignmentOperations.hh"

using boost::optional;
using graphtools::getQuerySequencesForEachNode;
using graphtools::getSequencesForEachOperation;
using graphtools::Graph;
using graphtools::NodeId;
using graphtools::Operation;
using graphtools::OperationType;
using graphtools::Path;
using std::list;
using std::string;
using std::unordered_map;
using std::vector;

static bool overlaps(const Segment& first, const Segment& second)
{
    const int64_t leftBound = first.start > second.start ? first.start : second.start;
    const int64_t rightBound = first.end < second.end ? first.end : second.end;

    return leftBound <= rightBound;
}

struct ReadAlignOrigin
{
    ReadAlignOrigin(std::string read, GraphAlign align, GenomicRegion origin, bool consistentWithMultipleHaplotypes)
        : read(std::move(read))
        , align(std::move(align))
        , origin(origin)
        , consistentWithMultipleHaplotypes(consistentWithMultipleHaplotypes)
    {
    }

    std::string read;
    GraphAlign align;
    GenomicRegion origin;
    bool consistentWithMultipleHaplotypes;
};

static bool singlePath(const vector<FragPathAlign>& fragPathAligns)
{
    optional<int> pathIndex;
    for (const auto& fragPathAlign : fragPathAligns)
    {
        if (!pathIndex)
        {
            pathIndex = fragPathAlign.readAlign.pathIndex;
        }

        if (fragPathAlign.readAlign.pathIndex != *pathIndex)
        {
            return false;
        }
    }

    return true;
}

list<ReadAlignOrigin> extractReadInfo(
    const FragAssignment& fragAssignment, const FragById& fragById, const FragPathAlignsById& fragPathAlignsById)
{
    list<ReadAlignOrigin> readInfo;
    for (int fragIndex = 0; fragIndex != fragAssignment.fragIds.size(); ++fragIndex)
    {
        const auto& fragId = fragAssignment.fragIds[fragIndex];
        const auto& frag = fragById.at(fragId);
        const int alignIndex = fragAssignment.alignIndexByFrag[fragIndex];
        const FragPathAlign& fragAlign = fragPathAlignsById.at(fragId)[alignIndex];

        const bool consistentWithMultiplePaths = !singlePath(fragPathAlignsById.at(fragId));

        GenomicRegion readRegion(fragAlign.readAlign.pathIndex, fragAlign.readAlign.begin, fragAlign.readAlign.end);
        readInfo.emplace_back(frag.read.bases, *fragAlign.readAlign.align, readRegion, consistentWithMultiplePaths);

        GenomicRegion mateRegion(fragAlign.mateAlign.pathIndex, fragAlign.mateAlign.begin, fragAlign.mateAlign.end);
        readInfo.emplace_back(frag.mate.bases, *fragAlign.mateAlign.align, mateRegion, consistentWithMultiplePaths);
    }
    return readInfo;
}

void clipFlanks(vector<Path>& hapPaths, int flankBasesToKeep, list<ReadAlignOrigin>& infoByRead)
{
    const auto& graph = *hapPaths.front().graphRawPtr();

    auto readInfoIt = infoByRead.begin();
    while (readInfoIt != infoByRead.end())
    {
        const int hapIndex = readInfoIt->origin.contigIndex();
        const auto& hapPath = hapPaths[hapIndex];
        const int leftFlankLen = graph.nodeSeq(hapPath.firstNodeId()).length();
        const int rightFlankLen = graph.nodeSeq(hapPath.lastNodeId()).length();

        const int clipRegionStart = leftFlankLen - flankBasesToKeep;
        const int clipRegionEnd = hapPath.length() - rightFlankLen + flankBasesToKeep;

        auto& origin = readInfoIt->origin;
        auto& align = readInfoIt->align;
        if (origin.end() <= clipRegionStart || origin.start() > clipRegionEnd)
        {
            readInfoIt = infoByRead.erase(readInfoIt);
        }
        else if (origin.start() < clipRegionStart && clipRegionStart <= origin.end())
        {
            align.shrinkStart(clipRegionStart - origin.start());
            origin.setStart(clipRegionStart);
        }
        else if (origin.start() < clipRegionEnd && clipRegionEnd < origin.end())
        {
            align.shrinkEnd(origin.end() - clipRegionEnd);
            origin.setEnd(clipRegionEnd);
        }
        else // TODO: Consider the case of both ends sticking out
        {
            ++readInfoIt;
        }
    }

    for (auto& readInfo : infoByRead)
    {
        auto& origin = readInfo.origin;
        const int hapIndex = origin.contigIndex();
        const auto& hapPath = hapPaths[hapIndex];
        const int leftFlankLen = graph.nodeSeq(hapPath.firstNodeId()).length();
        const int clipRegionStart = leftFlankLen - flankBasesToKeep;

        origin.setStart(origin.start() - clipRegionStart);
        origin.setEnd(origin.end() - clipRegionStart);
    }

    for (auto& hapPath : hapPaths)
    {
        const int leftFlankLen = graph.nodeSeq(hapPath.firstNodeId()).length();
        const int rightFlankLen = graph.nodeSeq(hapPath.lastNodeId()).length();
        hapPath.shrinkStartBy(leftFlankLen - flankBasesToKeep);
        hapPath.shrinkEndBy(rightFlankLen - flankBasesToKeep);
    }
}

class ColorPicker
{
public:
    explicit ColorPicker(const Graph& graph)
    {
        const vector<string> readLoopColors = { "#fc8d62", "#66c2a5" };
        const vector<string> pathLoopColors = { "url(#OrangeWhiteOrange)", "url(#GreenWhiteGreen)" };
        int loopIndex = 0;
        for (auto node = 0; node != graph.numNodes(); ++node)
        {
            if (graph.hasEdge(node, node))
            {
                readColorByNode_[node] = readLoopColors[loopIndex % readLoopColors.size()];
                pathColorByNode_[node] = pathLoopColors[loopIndex % pathLoopColors.size()];
                loopIndex += 1;
            }
            else
            {
                readColorByNode_[node] = "#8da0cb";
                pathColorByNode_[node] = "url(#BlueWhiteBlue)";
            }
        }
    }

    string getReadColor(NodeId node) const { return readColorByNode_.at(node); }
    string getPathColor(NodeId node) const { return pathColorByNode_.at(node); }

private:
    unordered_map<NodeId, std::string> readColorByNode_;
    unordered_map<NodeId, std::string> pathColorByNode_;
};

static optional<Feature>
getFeature(const ColorPicker& colorPicker, NodeId node, const Operation& op, const string& ref, const string& query)
{
    optional<Feature> feature;
    const string fill = colorPicker.getReadColor(node);
    const int opLength = op.length();
    const string atcg = "ATCG";
    if (op.type() == OperationType::kMatch)
    {
        for (string::size_type i; i < query.size(); i++)
        {
            feature = Feature(FeatureType::kRect, 1, fill, "none");
            if (atcg.find(ref[i]) != string::npos)
            {
                feature->label = string(ref[i], 1);
            }
        }
    }
    else if (op.type() == OperationType::kMismatch)
    {
        feature = Feature(FeatureType::kRect, opLength, "#cdcdcd", "none");
        feature->label = query;
    }
    else if (op.type() == OperationType::kDeletionFromRef)
    {
        feature = Feature(FeatureType::kLine, opLength, "none", "black");
    }
    else if (op.type() == OperationType::kSoftclip)
    {
        feature = Feature(FeatureType::kRect, opLength, "#cdcdcd", "none");
        feature->label = query;
    }
    else if (op.type() == OperationType::kInsertionToRef)
    {
        feature = Feature(FeatureType::kVerticalLine, 0, "none", "black");
    }
    else
    {
        assert(op.referenceLength() == 0);
    }

    return feature;
}

// TODO: Rename segment to displaySegment?
static list<Segment> getSegments(int hapIndex, const list<ReadAlignOrigin>& infoByRead)
{
    list<Segment> segments;
    ColorPicker colorPicker(*infoByRead.front().align.path().graphRawPtr());

    for (const auto& readInfo : infoByRead)
    {
        const auto& origin = readInfo.origin;
        if (origin.contigIndex() != hapIndex)
        {
            continue;
        }

        const auto& align = readInfo.align;

        const int segmentStart = [&]
        {
            int start = origin.start();
            const auto& firstOperation = align.alignments().front().front();
            if (firstOperation.type() == OperationType::kSoftclip)
            {
                start -= firstOperation.queryLength();
            }
            return start;
        }();

        vector<Feature> features;
        const int numNodes = align.path().numNodes();
        auto readPiecesByNode = getQuerySequencesForEachNode(align, readInfo.read);
        auto nodeReadPieceIt = readPiecesByNode.begin();

        for (int nodeIndex = 0; nodeIndex != numNodes; ++nodeIndex)
        {
            const auto node = align.path().getNodeIdByIndex(nodeIndex);
            const auto& nodeSeq = align.path().graphRawPtr()->nodeSeq(node);
            const auto& nodeAlign = align.alignments()[nodeIndex];
            const auto& nodeReadPiece = *nodeReadPieceIt;
            auto opReadPieces = getSequencesForEachOperation(nodeAlign, nodeSeq, nodeReadPiece);
            auto opReadPieceIt = opReadPieces.begin();

            for (const auto& operation : nodeAlign.operations())
            {
                const auto& opRefSeq = opReadPieceIt->first;
                const auto& opQuerySeq = opReadPieceIt->second;
                auto feature = getFeature(colorPicker, node, operation, opRefSeq, opQuerySeq);
                if (feature)
                {
                    features.push_back(*feature);
                }

                ++opReadPieceIt;
            }
            ++nodeReadPieceIt;
        }

        const double opacity = readInfo.consistentWithMultipleHaplotypes ? 0.7 : 1.0;
        segments.emplace_back(segmentStart, features, opacity);
    }

    return segments;
}

static list<Segment> trimSegments(const Path& hapPath, const list<Segment>& segments)
{
    list<Segment> trimmedSegments;
    const int haplotypeLen = hapPath.length();
    for (const auto& segment : segments)
    {
        vector<Feature> trimmedFeatures;
        int segmentPos = segment.start;
        auto featureIter = segment.features.begin();
        while (segmentPos + featureIter->length <= 0)
        {
            segmentPos += featureIter->length;
            ++featureIter;
        }

        if (segmentPos < 0)
        {
            const int trimmedLen = featureIter->length + segmentPos;
            trimmedFeatures.emplace_back(featureIter->type, trimmedLen, featureIter->fill, featureIter->stroke);
            if (featureIter->label)
            {
                trimmedFeatures.back().label = featureIter->label->substr(-segmentPos, trimmedLen);
            }
            ++featureIter;
        } // TODO: Address the case where a segment needs to be trimmed on both sides

        while (featureIter != segment.features.end() && segmentPos < haplotypeLen)
        {
            if (segmentPos + featureIter->length > haplotypeLen)
            {
                const int trimmedLen = haplotypeLen - segmentPos;
                trimmedFeatures.emplace_back(featureIter->type, trimmedLen, featureIter->fill, featureIter->stroke);
                if (featureIter->label)
                {
                    trimmedFeatures.back().label = featureIter->label->substr(0, trimmedLen);
                }
            }
            else
            {
                trimmedFeatures.push_back(*featureIter);
            }

            segmentPos += featureIter->length;
            ++featureIter;
        }

        const int trimmedSegmentStart = segment.start > 0 ? segment.start : 0;
        trimmedSegments.emplace_back(trimmedSegmentStart, trimmedFeatures, segment.opacity);
    }
    return trimmedSegments;
}

static void addLabelLane(const Path& path, LanePlot& lanePlot)
{
    vector<Segment> segments;
    const Graph& graph = *path.graphRawPtr();
    for (int node = 0; node != graph.numNodes(); ++node)
    {
        if (graph.hasEdge(node, node))
        {
            const int repeatSize = std::count(path.begin(), path.end(), node);
            if (repeatSize == 0)
            {
                continue;
            }
            const int start = path.getDistanceFromPathStart(node, 0);
            const int motifSize = graph.nodeSeq(node).length();
            // const int end = start + motifSize * repeatSize;

            vector<Feature> features = { { FeatureType::kArrows, motifSize * repeatSize, "black", "black" } };
            features.back().label = std::to_string(repeatSize) + " units";
            segments.emplace_back(start, features, 1.0);
        }
    }

    const int hapHeight = 12;
    lanePlot.emplace_back(hapHeight, segments);
}

static void addHaplotypePathLane(const Path& hapPath, LanePlot& lanePlot)
{
    // Add haplotype path
    vector<Feature> pathFeatures;
    const auto& graph = *hapPath.graphRawPtr();
    ColorPicker colorPicker(graph);
    for (int nodeIndex = 0; nodeIndex != static_cast<int>(hapPath.nodeIds().size()); ++nodeIndex)
    {
        const auto& nodeSeq = hapPath.getNodeSeq(nodeIndex);
        const int featureLen = hapPath.getNodeOverlapLengthByIndex(nodeIndex);
        const auto& fill = colorPicker.getPathColor(hapPath.getNodeIdByIndex(nodeIndex));
        if (nodeIndex == 0)
        {
            pathFeatures.emplace_back(FeatureType::kRectWithLeftBreak, featureLen, fill, "black");
        }
        else if (nodeIndex == hapPath.nodeIds().size() - 1)
        {
            pathFeatures.emplace_back(FeatureType::kRectWithRightBreak, featureLen, fill, "black");
        }
        else
        {
            pathFeatures.emplace_back(FeatureType::kRect, featureLen, fill, "black");
        }
        pathFeatures.back().label = nodeSeq;
    }
    const int hapHeight = 18;
    Segment pathSegment(0, pathFeatures, 1.0);
    lanePlot.emplace_back(hapHeight, vector<Segment>({ pathSegment }));
}

static void addSegmentLanes(list<Segment>& hapSegments, LanePlot& lanePlot)
{
    using IntTuple = std::tuple<int, int>;
    hapSegments.sort([](const Segment& lhs, const Segment& rhs)
                     { return IntTuple(lhs.end, lhs.start) < IntTuple(rhs.end, rhs.start); });

    while (!hapSegments.empty())
    {
        vector<Segment> laneSegments;
        auto iter = hapSegments.begin();
        while (iter != hapSegments.end())
        {
            if (laneSegments.empty())
            {
                laneSegments.push_back(*iter);
                iter = hapSegments.erase(iter);
            }
            else
            {
                const Segment& lastSegment = laneSegments.back();
                if (!overlaps(lastSegment, *iter))
                {
                    laneSegments.push_back(*iter);
                    iter = hapSegments.erase(iter);
                }
                else
                {
                    ++iter;
                }
            }
        }
        const int readHeight = 10;
        lanePlot.emplace_back(readHeight, std::move(laneSegments));
    }
}

static LanePlot getLanePlot(const Path& path, list<Segment>& hapSegments)
{
    LanePlot lanePlot;
    addLabelLane(path, lanePlot);
    addHaplotypePathLane(path, lanePlot);
    addSegmentLanes(hapSegments, lanePlot);

    return lanePlot;
}

void removeFlankingReads(list<ReadAlignOrigin>& infoByRead)
{
    const Graph& graph = *infoByRead.front().align.path().graphRawPtr();
    const NodeId leftFlankNode = 0;
    const NodeId rightFlankNode = graph.numNodes() - 1;
    auto readInfoIt = infoByRead.begin();
    while (readInfoIt != infoByRead.end())
    {
        const NodeId firstNode = readInfoIt->align.path().firstNodeId();
        const NodeId lastNode = readInfoIt->align.path().lastNodeId();
        const bool insideLeftFlank = lastNode == leftFlankNode;
        const bool insideRightFlank = firstNode == rightFlankNode;
        if (insideLeftFlank || insideRightFlank)
        {
            readInfoIt = infoByRead.erase(readInfoIt);
        }
        else
        {
            ++readInfoIt;
        }
    }
}

vector<LanePlot> generateBlueprint(
    vector<Path> paths, const FragById& fragById, const FragAssignment& fragAssignment,
    const FragPathAlignsById& fragPathAlignsById)
{
    auto infoByRead = extractReadInfo(fragAssignment, fragById, fragPathAlignsById);
    removeFlankingReads(infoByRead);
    clipFlanks(paths, 50, infoByRead);
    vector<LanePlot> lanePlots;
    for (int pathIndex = 0; pathIndex != paths.size(); ++pathIndex)
    {
        auto segments = getSegments(pathIndex, infoByRead);
        segments = trimSegments(paths[pathIndex], segments);
        auto lanePlot = getLanePlot(paths[pathIndex], segments);
        lanePlots.push_back(lanePlot);
    }

    return lanePlots;
}
