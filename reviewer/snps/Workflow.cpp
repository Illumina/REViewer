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

#include "graphalign/GraphAlignmentOperations.hh"
#include "graphalign/LinearAlignmentOperations.hh"

using graphtools::Alignment;
using graphtools::getQuerySequencesForEachNode;
using graphtools::getSequencesForEachOperation;
using graphtools::Operation;
using graphtools::OperationType;
using std::string;
using std::vector;

using AlignedRead = std::pair<std::string, ReadPathAlign>;

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

void updateBaseCounts(
    const string& bases, const string& nodeSeq, const Alignment& align, int startPos, vector<BaseCounts>& baseCountVec)
{
    auto seqsByOperation = getSequencesForEachOperation(align, nodeSeq, bases);
    auto seqsByOperationIt = seqsByOperation.begin();

    for (const Operation& operation : align)
    {
        const auto& operationReferenceSeq = seqsByOperationIt->first;
        const auto& operationBaseSeq = seqsByOperationIt->second;
        switch (operation.type())
        {
        case OperationType::kMatch:
        case OperationType::kMismatch:
            for (int offset = 0; offset != operation.referenceLength(); ++offset)
            {
                auto base = operationBaseSeq[offset];
                baseCountVec[startPos + offset].counts[getBaseIndex(base)] += 1;
            }
            break;
        default:
            break;
        }

        ++seqsByOperationIt;
        startPos += static_cast<int>(operation.referenceLength());
    }
}

vector<BaseCounts> getBaseCountVector(const GraphPath& path, const vector<AlignedRead>& alignedReads)
{
    vector<BaseCounts> baseCountVec;
    baseCountVec.resize(path.length());

    for (const auto& alignedRead : alignedReads)
    {
        const auto& bases = alignedRead.first;
        const auto& pathAlign = alignedRead.second;
        const int numNodesInAlign = static_cast<int>(pathAlign.align->path().numNodes());
        auto readPiecesByNode = getQuerySequencesForEachNode(*pathAlign.align, bases);
        auto nodeReadPieceIt = readPiecesByNode.begin();

        int posOnPath = pathAlign.begin;

        for (int nodeIndex = 0; nodeIndex != numNodesInAlign; ++nodeIndex)
        {
            const auto node = pathAlign.align->path().getNodeIdByIndex(nodeIndex);
            const auto& nodeSeq = pathAlign.align->path().graphRawPtr()->nodeSeq(node);
            const auto& nodeAlign = pathAlign.align->alignments()[nodeIndex];

            updateBaseCounts(*nodeReadPieceIt, nodeSeq, nodeAlign, posOnPath, baseCountVec);

            ++nodeReadPieceIt;
            posOnPath += static_cast<int>(pathAlign.align->referenceLength());
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
        vector<AlignedRead> haplotypeReads;
        for (int fragIndex = 0; fragIndex != fragAssignment.fragIds.size(); ++fragIndex)
        {
            const auto& fragId = fragAssignment.fragIds[fragIndex];
            const int alignIndex = fragAssignment.alignIndexByFrag[fragIndex];
            const auto& fragPathAlign = fragPathAlignsById.at(fragId)[alignIndex];
            assert(fragPathAlign.readAlign.align && fragPathAlign.mateAlign.align);

            if (fragPathAlign.readAlign.pathIndex == hapIndex)
            {
                const auto& readBases = fragById.at(fragId).read.bases;
                const auto& mateBases = fragById.at(fragId).mate.bases;
                haplotypeReads.emplace_back(std::make_pair(readBases, fragPathAlign.readAlign));
                haplotypeReads.emplace_back(std::make_pair(mateBases, fragPathAlign.mateAlign));
            }
        }

        std::cerr << "Hap\tRef\tA\tT\tC\tG" << std::endl;
        auto baseCountVec = getBaseCountVector(paths[hapIndex], haplotypeReads);
        const auto& haplotypeSeq = paths[hapIndex].seq();

        for (int index = 0; index != baseCountVec.size(); ++index)
        {
            const auto& baseCounts = baseCountVec[index];
            char ref = haplotypeSeq[index];
            auto as = baseCounts.getCount('A');
            auto ts = baseCounts.getCount('T');
            auto cs = baseCounts.getCount('C');
            auto gs = baseCounts.getCount('G');
            std::cerr << hapIndex << "\t" << ref << "\t" << as << "\t" << ts << "\t" << cs << "\t" << gs << std::endl;
        }
    }
}

}