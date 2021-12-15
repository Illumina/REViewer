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

#include "Workflow.hh"

#include <fstream>
#include <set>

#include <boost/algorithm/string.hpp>

#include "spdlog/spdlog.h"
#include "thirdparty/json/json.hpp"

#include "app/Aligns.hh"
#include "app/CatalogLoading.hh"
#include "app/FragLenFilter.hh"
#include "app/GenerateSvg.hh"
#include "app/GenotypePaths.hh"
#include "app/LanePlot.hh"
#include "app/Origin.hh"
#include "app/Phasing.hh"
#include "app/Projection.hh"
#include "metrics/Metrics.hh"

using boost::optional;
using graphtools::Graph;
using graphtools::GraphAlignment;
using graphtools::Operation;
using graphtools::OperationType;
using std::map;
using std::ofstream;
using std::string;
using std::unordered_map;
using std::vector;

class LocusResults
{
public:
    LocusResults(
        const ScoredDiplotypes& scoredDiplotypes, vector<LanePlot> lanePlots, MetricsByVariant metricsByVariant)
        : scoredDiplotypes_(std::move(scoredDiplotypes))
        , lanePlots_(std::move(lanePlots))
        , metricsByVariant_(std::move(metricsByVariant))
    {
    }

    const ScoredDiplotypes& scoredDiplotypes() const { return scoredDiplotypes_; }
    const vector<LanePlot>& lanePlots() const { return lanePlots_; }
    const MetricsByVariant& metricsByVariant() const { return metricsByVariant_; }

private:
    ScoredDiplotypes scoredDiplotypes_;
    vector<LanePlot> lanePlots_;
    MetricsByVariant metricsByVariant_;
};

using ResultsByLocus = map<string, LocusResults>;

static LocusResults analyzeLocus(
    const string& referencePath, const string& readsPath, const string& vcfPath, const string& locusId,
    const LocusSpecification& locusSpec)
{
    spdlog::info("Loading specification of locus {}", locusId);

    auto fragById = getAligns(readsPath, referencePath, locusSpec);
    spdlog::info("Extracted {} frags", fragById.size());

    spdlog::info("Calculating fragment length");
    const int meanFragLen = getMeanFragLen(fragById);
    spdlog::info("Fragment length is estimated to be {}", meanFragLen);

    spdlog::info("Extracting genotype paths");
    auto pathsByDiplotype = getCandidateDiplotypePaths(meanFragLen, vcfPath, locusSpec);

    spdlog::info("Phasing");

    auto scoredDiplotypes = scoreDiplotypes(fragById, pathsByDiplotype);
    auto topDiplotype = scoredDiplotypes.front().first; // scoredDiplotypes are sorted
    spdlog::info("Found {} paths defining diplotype", topDiplotype.size());

    spdlog::info("Projecting reads onto haplotype paths");
    auto pairPathAlignById = project(topDiplotype, fragById);
    spdlog::info("Projected {} read pairs", pairPathAlignById.size());

    spdlog::info("Generating fragment alignments");
    auto fragPathAlignsById = resolveByFragLen(meanFragLen, topDiplotype, pairPathAlignById);
    spdlog::info("Generated {} fragment alignments", fragPathAlignsById.size());

    spdlog::info("Assigning fragment origins");
    auto fragAssignment = getBestFragAssignment(topDiplotype, fragPathAlignsById);
    spdlog::info("Found assignments for {} frags", fragAssignment.fragIds.size());

    spdlog::info("Generating metrics");
    auto metricsByVariant = getMetrics(locusSpec, topDiplotype, fragById, fragAssignment, fragPathAlignsById);

    spdlog::info("Generating plot blueprint");
    auto lanePlots = generateBlueprint(topDiplotype, fragById, fragAssignment, fragPathAlignsById);

    return { scoredDiplotypes, lanePlots, metricsByVariant };
}

vector<string> getLocusIds(const RegionCatalog& catalog, const string& encoding)
{
    vector<string> locusIds;
    boost::split(locusIds, encoding, boost::is_any_of(","));

    for (const auto& locusId : locusIds)
    {
        if (locusId.empty())
        {
            throw std::runtime_error("Empty locus ids are not allowed");
        }

        if (catalog.find(locusId) == catalog.end())
        {
            throw std::runtime_error(locusId + " is missing from the variant catalog");
        }
    }

    return locusIds;
}

static string summarizePath(const graphtools::Graph& graph, const graphtools::Path& path)
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
            assert(graph.numNodes() != 0);
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

void outputPhasingInfo(const Graph& graph, const ScoredDiplotypes& scoredDiplotypes, const string& phasingInfoPath)
{
    using Diplotypes = std::set<graphtools::Path>;

    ofstream phasingInfoFile(phasingInfoPath);
    if (phasingInfoFile.is_open())
    {
        std::set<Diplotypes> observedGenotypes;
        for (const auto& scoredGenotype : scoredDiplotypes)
        {
            const auto& genotypePaths = scoredGenotype.first;
            const int score = scoredGenotype.second;

            Diplotypes genotypePathSet;
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

void outputResults(
    const RegionCatalog& locusCatalog, const ResultsByLocus& resultsByLocus, const string& outputPrefix,
    bool phasingInfoNeeded)
{
    spdlog::info("Writing output to disk");
    nlohmann::json metricsReport;
    for (const auto& locusIdAndResults : resultsByLocus)
    {
        const auto& locusId = locusIdAndResults.first;
        const auto& locusResults = locusIdAndResults.second;
        const auto& locusSpec = locusCatalog.at(locusId);
        const string locusOutputPrefix = outputPrefix + "." + locusId;
        generateSvg(locusResults.lanePlots(), locusOutputPrefix + ".svg");
        if (phasingInfoNeeded)
        {
            const auto phasingInfoPath = locusOutputPrefix + ".phasing.txt";
            const Graph& graph = locusSpec.regionGraph();
            outputPhasingInfo(graph, locusResults.scoredDiplotypes(), phasingInfoPath);
        }

        for (const auto& variantMetrics : locusResults.metricsByVariant())
        {
            metricsReport[variantMetrics.variantId]["Genotype"] = variantMetrics.genotype;
            metricsReport[variantMetrics.variantId]["AlleleDepths"] = variantMetrics.alleleDepth;
        }
    }

    const string metricsFilePath = outputPrefix + ".metrics.json";
    ofstream metricsFile(metricsFilePath);
    if (metricsFile.is_open())
    {
        metricsFile << metricsReport.dump(4) << std::endl;
    }
    else
    {
        throw std::runtime_error("Unable to open " + metricsFilePath);
    }
    metricsFile.close();
}

int runWorkflow(const WorkflowArguments& args)
{
    Reference reference(args.referencePath);
    auto locusCatalog = loadLocusCatalogFromDisk(args.catalogPath, reference, args.locusExtensionLength);

    ResultsByLocus resultsByLocus;
    auto locusIds = getLocusIds(locusCatalog, args.locusId);
    for (const auto& locusId : locusIds)
    {
        auto locusSpec = locusCatalog.at(locusId);
        auto locusResults = analyzeLocus(args.referencePath, args.readsPath, args.vcfPath, locusId, locusSpec);
        resultsByLocus.emplace(locusId, locusResults);
    }

    outputResults(locusCatalog, resultsByLocus, args.outputPrefix, args.outputPhasingInfo);

    return 0;
}
