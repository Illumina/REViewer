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
#include <vector>

using std::vector;

static void
getAlleleDepth(const LocusSpecification& locusSpec, const GraphPath& hapPath, const vector<GraphAlignPtr>& hapAligns)
{
    for (const auto& variantSpec : locusSpec.variantSpecs())
    {
        assert(variantSpec.nodes().size() == 1);
        const auto strNode = variantSpec.nodes().front();

        const int strSize = static_cast<int>(std::count(hapPath.nodeIds().begin(), hapPath.nodeIds().end(), strNode));

        std::cerr << variantSpec.id() << " " << strNode << " " << strSize << std::endl;
        for (const auto& alignPtr : hapAligns)
        {
            for (auto nodeIndex = 0; nodeIndex != static_cast<int>(alignPtr->path().numNodes()); ++nodeIndex)
            {
                auto node = alignPtr->getNodeIdByIndex(nodeIndex);
                if (node != strNode)
                {
                    continue;
                }

                std::cerr << alignPtr->alignments()[nodeIndex] << std::endl;
            }
        }
    }
}

void getMetrics(
    const LocusSpecification& locusSpec, const GraphPaths& paths, const FragById& fragById,
    const FragAssignment& fragAssignment, const FragPathAlignsById& fragPathAlignsById)
{
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
        getAlleleDepth(locusSpec, paths[hapIndex], hapAligns);
    }
}
