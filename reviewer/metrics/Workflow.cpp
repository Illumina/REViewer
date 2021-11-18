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

#include "metrics/Workflow.hh"

#include <algorithm>
#include <fstream>
#include <map>
#include <sstream>
#include <vector>

using graphtools::NodeId;
using std::map;
using std::ofstream;
using std::string;
using std::vector;

static int getTotalMatchesToNode(NodeId targetNode, const vector<GraphAlignPtr>& hapAligns)
{
    int numMatches = 0;
    for (const auto& alignPtr : hapAligns)
    {
        for (auto nodeIndex = 0; nodeIndex != static_cast<int>(alignPtr->path().numNodes()); ++nodeIndex)
        {
            auto node = alignPtr->getNodeIdByIndex(nodeIndex);
            if (targetNode != node)
            {
                continue;
            }

            numMatches += static_cast<int>(alignPtr->alignments()[nodeIndex].numMatched());
        }
    }

    return numMatches;
}

static void getAlleleDepth(
    const LocusSpecification& locusSpec, const GraphPath& hapPath, const vector<GraphAlignPtr>& hapAligns,
    map<string, vector<double>>* alleleDepths)
{
    for (const auto& variantSpec : locusSpec.variantSpecs())
    {
        assert(variantSpec.nodes().size() == 1);
        const auto strNode = variantSpec.nodes().front();
        int strLen = static_cast<int>(std::count(hapPath.nodeIds().begin(), hapPath.nodeIds().end(), strNode));
        strLen *= static_cast<int>(locusSpec.regionGraph().nodeSeq(strNode).length());
        const int numMatches = getTotalMatchesToNode(strNode, hapAligns);
        const double depth = strLen > 0 ? numMatches / static_cast<double>(strLen) : 0.0;
        (*alleleDepths)[variantSpec.id()].push_back(depth);
    }
}

static map<string, vector<double>> getAlleleDepths(
    const LocusSpecification& locusSpec, const GraphPaths& paths, const FragById& fragById,
    const FragAssignment& fragAssignment, const FragPathAlignsById& fragPathAlignsById)
{
    map<string, vector<double>> alleleDepths;

    for (int hapIndex = 0; hapIndex != paths.size(); ++hapIndex)
    {
        vector<GraphAlignPtr> hapAligns;
        for (int fragIndex = 0; fragIndex != fragAssignment.fragIds.size(); ++fragIndex)
        {
            const auto& fragId = fragAssignment.fragIds[fragIndex];
            const int alignIndex = fragAssignment.alignIndexByFrag[fragIndex];
            const auto& fragPathAlign = fragPathAlignsById.at(fragId)[alignIndex];
            assert(fragPathAlign.readAlign.align && fragPathAlign.mateAlign.align);

            if (fragPathAlign.readAlign.pathIndex == hapIndex)
            {
                hapAligns.push_back(fragPathAlign.readAlign.align);
                hapAligns.push_back(fragPathAlign.mateAlign.align);
            }
        }
        getAlleleDepth(locusSpec, paths[hapIndex], hapAligns, &alleleDepths);
    }

    return alleleDepths;
}

static string encodeDepths(const vector<double>& depths)
{
    std::ostringstream encoding;
    encoding.precision(2);
    encoding << std::fixed;

    for (auto depth : depths)
    {
        if (encoding.tellp())
        {
            encoding << "/";
        }
        encoding << depth;
    }

    return encoding.str();
}

void getMetrics(
    const LocusSpecification& locusSpec, const GraphPaths& paths, const FragById& fragById,
    const FragAssignment& fragAssignment, const FragPathAlignsById& fragPathAlignsById, const string& outputPrefix)
{
    const auto alleleDepths = getAlleleDepths(locusSpec, paths, fragById, fragAssignment, fragPathAlignsById);

    const string alleleDepthsPath = outputPrefix + ".allele_depth.tsv";
    ofstream alleleDepthsFile(alleleDepthsPath);
    if (alleleDepthsFile.is_open())
    {
        alleleDepthsFile << "VariantId\tAlleleDepths" << std::endl;
        for (const auto& variantSpec : locusSpec.variantSpecs())
        {
            const auto depthsEncoding = encodeDepths(alleleDepths.at(variantSpec.id()));
            alleleDepthsFile << variantSpec.id() << "\t" << depthsEncoding << std::endl;
        }
    }
    else
    {
        throw std::runtime_error("Unable to open " + alleleDepthsPath);
    }
    alleleDepthsFile.close();
}
