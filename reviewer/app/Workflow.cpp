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

#include <boost/algorithm/string.hpp>

#include "spdlog/spdlog.h"

#include "app/Aligns.hh"
#include "app/CatalogLoading.hh"
#include "app/FragLenFilter.hh"
#include "app/GenerateSvg.hh"
#include "app/GenotypePaths.hh"
#include "app/LanePlot.hh"
#include "app/Origin.hh"
#include "app/Phasing.hh"
#include "app/Projection.hh"
#include "metrics/Workflow.hh"
#include "snps/Workflow.hh"

using boost::optional;
using graphtools::GraphAlignment;
using graphtools::Operation;
using graphtools::OperationType;
using std::string;
using std::unordered_map;
using std::vector;

static void analyzeLocus(
    const string& referencePath, const string& readsPath, const string& vcfPath, const string& locusId,
    const LocusSpecification& locusSpec, const string& outputPrefix, bool outputPhasingInfo)
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
    optional<string> phasingInfoPath;
    if (outputPhasingInfo)
    {
        phasingInfoPath = outputPrefix + ".phasing.txt";
    }

    auto diplotypePaths = phase(fragById, pathsByDiplotype, phasingInfoPath);
    spdlog::info("Found {} paths defining diplotype", diplotypePaths.size());

    spdlog::info("Projecting reads onto haplotype paths");
    auto pairPathAlignById = project(diplotypePaths, fragById);
    spdlog::info("Projected {} read pairs", pairPathAlignById.size());

    spdlog::info("Generating fragment alignments");
    auto fragPathAlignsById = resolveByFragLen(meanFragLen, diplotypePaths, pairPathAlignById);
    spdlog::info("Generated {} fragment alignments", fragPathAlignsById.size());

    spdlog::info("Assigning fragment origins");
    auto fragAssignment = getBestFragAssignment(diplotypePaths, fragPathAlignsById);
    spdlog::info("Found assignments for {} frags", fragAssignment.fragIds.size());

    spdlog::info("Generating metrics");
    getMetrics(locusSpec, diplotypePaths, fragById, fragAssignment, fragPathAlignsById, outputPrefix);

    spdlog::info("Calling SNPs");
    snps::callSnps(diplotypePaths, fragById, fragAssignment, fragPathAlignsById);

    spdlog::info("Generating plot blueprint");
    auto lanePlots = generateBlueprint(diplotypePaths, fragById, fragAssignment, fragPathAlignsById);

    spdlog::info("Writing SVG image to disk");
    generateSvg(lanePlots, outputPrefix + ".svg");
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

int runWorkflow(const WorkflowArguments& args)
{
    const int kFlankLength = 1000;
    Reference reference(args.referencePath);
    auto locusCatalog = loadLocusCatalogFromDisk(args.catalogPath, reference, args.locusExtensionLength);

    auto locusIds = getLocusIds(locusCatalog, args.locusId);
    for (const auto& locusId : locusIds)
    {
        const string locusOutputPrefix = args.outputPrefix + "." + locusId;
        auto locusSpec = locusCatalog.at(locusId);
        analyzeLocus(
            args.referencePath, args.readsPath, args.vcfPath, locusId, locusSpec, locusOutputPrefix,
            args.outputPhasingInfo);
    }

    return 0;
}
