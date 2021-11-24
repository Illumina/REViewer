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

#include "snps/Workflow.hh"

#include <array>
#include <stdexcept>
#include <vector>

using graphtools::Operation;
using graphtools::OperationType;
using std::vector;

namespace snps
{

int getBaseIndex(char base)
{
    switch (base)
    {
    case 'A':
    case 'a':
        return 0;
    case 'C':
    case 'c':
        return 1;
    case 'T':
    case 't':
        return 2;
    case 'G':
    case 'g':
        return 3;
    default:
        throw std::runtime_error("Encountered unknown base " + std::to_string(base));
    }
}

struct BaseCounts
{
    BaseCounts()
        : counts({ 0, 0, 0, 0 })
    {
    }

    int getCount(char base) const { return counts[getBaseIndex(base)]; }

    std::array<int, 4> counts;
};

vector<BaseCounts> getBaseCountVector(const GraphPath& path, const vector<const ReadPathAlign*>& alignPtrs)
{
    vector<BaseCounts> baseCountVec;
    baseCountVec.resize(path.length());

    for (const auto& pathAlignPtr : alignPtrs)
    {
        const auto& pathAlign = *pathAlignPtr;
        // std::cerr << "\t" << *pathAlign.align << " (" << pathAlign.begin << "-" << pathAlign.end << ")" << std::endl;

        int posOnPath = pathAlign.begin;
        for (const auto& nodeAlign : pathAlign.align->alignments())
        {
            for (const Operation& operation : nodeAlign)
            {
                switch (operation.type())
                {
                case OperationType::kMatch:
                    for (int pos = posOnPath; pos != posOnPath + operation.referenceLength(); ++pos)
                    {
                        baseCountVec[pos].counts[getBaseIndex('A')] += 1;
                    }
                    break;
                /*
                case OperationType::kMismatch:
                    score += mismatchScore * operation.referenceLength();
                    break;
                case OperationType::kInsertionToRef:
                    score += gapScore * operation.queryLength();
                    break;
                case OperationType::kDeletionFromRef:
                    score += gapScore * operation.referenceLength();
                    break; */
                default:
                    break;
                }
            }

            // auto startNodeIndex = pathAlign.startIndexOnPath;
            // auto endNodeIndex = startNodeIndex + pathAlign.align->path().numNodes();
            // for (int nodeIndex = startNodeIndex; nodeIndex != endNodeIndex) {
            //
            // }
        }
    }
    return baseCountVec;
}

void callSnps(
    const GraphPaths& paths, const FragById& fragById, const FragAssignment& fragAssignment,
    const FragPathAlignsById& fragPathAlignsById)
{
    for (int hapIndex = 0; hapIndex != paths.size(); ++hapIndex)
    {
        vector<const ReadPathAlign*> hapAligns;
        for (int fragIndex = 0; fragIndex != fragAssignment.fragIds.size(); ++fragIndex)
        {
            const auto& fragId = fragAssignment.fragIds[fragIndex];
            const int alignIndex = fragAssignment.alignIndexByFrag[fragIndex];
            const auto& fragPathAlign = fragPathAlignsById.at(fragId)[alignIndex];
            assert(fragPathAlign.readAlign.align && fragPathAlign.mateAlign.align);

            if (fragPathAlign.readAlign.pathIndex == hapIndex)
            {
                hapAligns.push_back(&fragPathAlign.readAlign);
                hapAligns.push_back(&fragPathAlign.mateAlign);
            }
        }

        auto baseCountVec = getBaseCountVector(paths[hapIndex], hapAligns);
        for (const auto& baseCounts : baseCountVec)
        {
            auto aCount = baseCounts.getCount('A');
            auto tCount = baseCounts.getCount('T');
            auto cCount = baseCounts.getCount('C');
            auto gCount = baseCounts.getCount('G');
            std::cerr << hapIndex << "\t" << aCount << "\t" << tCount << "\t" << cCount << "\t" << gCount << std::endl;
        }
    }
}

}