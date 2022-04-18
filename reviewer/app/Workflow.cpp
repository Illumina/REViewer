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

#include <set>
#include <stdlib.h>

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

template <typename T> static string encode(const vector<T>& depths)
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
    auto pathsByDiplotype = getCandidateDiplotypes(meanFragLen, vcfPath, locusSpec);

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

std::ofstream initPhasingFile(const std::string& prefix)
{
    const auto path = prefix + ".phasing.tsv";
    std::ofstream file(path);
    if (!file.is_open())
    {
        throw std::runtime_error("Unable to open " + path);
    }

    file << "LocusId\tDiplotype\tScore" << std::endl;

    return file;
}

std::ofstream initMetricsFile(const std::string& prefix)
{
    const auto path = prefix + ".metrics.tsv";
    std::ofstream metricsFile(path);
    if (!metricsFile.is_open())
    {
        throw std::runtime_error("Unable to open " + path);
    }

    metricsFile << "VariantId\tGenotype\tAlleleDepth" << std::endl;

    return metricsFile;
}

int runWorkflow(const WorkflowArguments& args)
{
    Reference reference(args.referencePath);
    auto locusCatalog = loadLocusCatalogFromDisk(args.catalogPath, reference, args.locusExtensionLength);
    auto locusIds = getLocusIds(locusCatalog, args.locusId);
    auto phasingFile = initPhasingFile(args.outputPrefix);
    auto metricsFile = initMetricsFile(args.outputPrefix);
    
    // For reproducibility
    srand(14345);

    for (const auto& locusId : locusIds)
    {
        auto locusSpec = locusCatalog.at(locusId);
        auto locusResults = analyzeLocus(args.referencePath, args.readsPath, args.vcfPath, locusId, locusSpec);

        const auto svgPath = args.outputPrefix + "." + locusId + ".svg";
        generateSvg(locusResults.lanePlots(), svgPath);

        for (const auto& metrics : locusResults.metricsByVariant())
        {
            const auto& genotype = encode(metrics.genotype);
            const auto& alleleDepth = encode(metrics.alleleDepth);
            metricsFile << metrics.variantId << "\t" << genotype << "\t" << alleleDepth << std::endl;
        }

        for (const auto& diplotypeAndScore : locusResults.scoredDiplotypes())
        {
            phasingFile << locusId << "\t" << diplotypeAndScore.first << "\t" << diplotypeAndScore.second << std::endl;
        }
    }

    phasingFile.close();
    metricsFile.close();

    return 0;
}
