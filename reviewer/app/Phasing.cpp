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

#include "app/Phasing.hh"

#include <algorithm>
#include <fstream>
#include <set>
#include <utility>

#include "app/Projection.hh"

using boost::optional;
using std::ofstream;
using std::pair;
using std::string;
using std::vector;

int scorePath(int pathIndex, const PairPathAlignById& pairPathAlignById)
{
    int pathScore = 0;
    for (const auto& idAndPairPathAlign : pairPathAlignById)
    {
        const auto& pairAlign = idAndPairPathAlign.second;
        for (const auto& readAlign : pairAlign.readAligns)
        {
            if (readAlign.pathIndex == pathIndex)
            {
                pathScore += score(*readAlign.align);
                break;
            }
        }

        for (const auto& mateAlign : pairAlign.mateAligns)
        {
            if (mateAlign.pathIndex == pathIndex)
            {
                pathScore += score(*mateAlign.align);
                break;
            }
        }
    }

    return pathScore;
}

using ScoredGenotype = pair<DiplotypePaths, int>;

string summarizePath(const graphtools::Graph& graph, const graphtools::Path& path)
{
    string summary;
    std::set<graphtools::NodeId> observedNodes;
    for (const auto nodeId : path.nodeIds())
    {
        if (observedNodes.find(nodeId) != observedNodes.end())
        {
            continue;
        }

        const bool isLoopNode = graph.hasEdge(nodeId, nodeId);

        if (nodeId == 0)
        {
            summary += "(LF)";
        }
        else if (nodeId + 1 == graph.numNodes())
        {
            summary += "(RF)";
        }
        else
        {
            const string& nodeSeq = graph.nodeSeq(nodeId);
            summary += "(" + nodeSeq + ")";

            if (isLoopNode)
            {
                int numMotifs = std::count(path.nodeIds().begin(), path.nodeIds().end(), nodeId);
                summary += "{" + std::to_string(numMotifs) + "}";
            }
        }

        observedNodes.emplace(nodeId);
    }

    return summary;
}

void outputPhasingInfo(
    const FragById& fragById, const vector<ScoredGenotype>& scoredGenotypes, const string& phasingInfoPath)
{
    using GenotypePathSet = std::set<graphtools::Path>;
    assert(!fragById.empty());
    const auto& firstFrag = fragById.begin()->second;
    const graphtools::Graph& graph = *firstFrag.read.align.path().graphRawPtr();

    ofstream phasingInfoFile(phasingInfoPath);
    if (phasingInfoFile.is_open())
    {
        std::set<GenotypePathSet> observedGenotypes;
        for (const auto& scoredGenotype : scoredGenotypes)
        {
            const auto& genotypePaths = scoredGenotype.first;
            const int score = scoredGenotype.second;

            GenotypePathSet genotypePathSet;
            genotypePathSet.emplace(genotypePaths.front());
            genotypePathSet.emplace(genotypePaths.back());
            if (observedGenotypes.find(genotypePathSet) != observedGenotypes.end())
            {
                continue;
            }

            phasingInfoFile << summarizePath(graph, genotypePaths.front());
            if (genotypePaths.size() == 2)
            {
                phasingInfoFile << "/" << summarizePath(graph, genotypePaths.back());
            }
            phasingInfoFile << "\t" << score << std::endl;
            observedGenotypes.emplace(genotypePathSet);
        }
    }
    else
    {
        throw std::runtime_error("Unable to open " + phasingInfoPath);
    }
}

DiplotypePaths phase(
    const FragById& fragGraphAlignById, const vector<DiplotypePaths>& pathsByGenotype,
    const optional<string>& phasingInfoPath)
{
    vector<ScoredGenotype> scoredGenotypes;

    for (const auto& genotypePaths : pathsByGenotype)
    {
        assert(genotypePaths.size() == 1 || genotypePaths.size() == 2);

        auto pairPathAlignById = project(genotypePaths, fragGraphAlignById);

        int genotypeScore = scorePath(0, pairPathAlignById);
        if (genotypePaths.size() == 2)
        {
            genotypeScore += scorePath(1, pairPathAlignById);
        }

        scoredGenotypes.emplace_back(genotypePaths, genotypeScore);
    }

    std::sort(
        scoredGenotypes.begin(), scoredGenotypes.end(),
        [](const ScoredGenotype& gt1, const ScoredGenotype& gt2) { return gt1.second > gt2.second; });

    if (phasingInfoPath)
    {
        outputPhasingInfo(fragGraphAlignById, scoredGenotypes, *phasingInfoPath);
    }

    assert(!scoredGenotypes.empty());
    return scoredGenotypes.front().first;
}
