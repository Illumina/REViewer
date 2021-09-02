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

using boost::optional;
using graphtools::GraphAlignment;
using graphtools::Operation;
using graphtools::OperationType;
using std::string;
using std::unordered_map;
using std::vector;

void visualizeLocus(
    const string& referencePath, const string& readsPath, const string& vcfPath, const string& locusId,
    const LocusSpecification& locusSpec, const string& svgPath, const optional<string>& phasingInfoPath)
{
    spdlog::info("Loading specification of locus {}", locusId);

    ReadPairById fragById;
    auto fragGraphAlignById = getAligns(readsPath, referencePath, locusSpec, fragById);
    spdlog::info("Extracted {} frags", fragGraphAlignById.size());

    spdlog::info("Calculating fragment length");
    const int meanFragLen = getMeanFragLen(fragGraphAlignById);
    spdlog::info("Fragment length is estimated to be {}", meanFragLen);

    spdlog::info("Extracting genotype paths");
    auto pathsByGenotype = getCandidateGenotypePaths(meanFragLen, vcfPath, locusSpec);

    spdlog::info("Phasing");
    auto genotypePaths = phase(fragGraphAlignById, pathsByGenotype, phasingInfoPath);
    spdlog::info("Found {} paths defining genotype", genotypePaths.size());

    spdlog::info("Projecting reads onto haplotype paths");
    auto pairPathAlignById = project(genotypePaths, fragGraphAlignById);
    spdlog::info("Projected {} read pairs", pairPathAlignById.size());

    spdlog::info("Generating fragment alignments");
    auto fragPathAlignsById = resolveByFragLen(meanFragLen, genotypePaths, pairPathAlignById);
    spdlog::info("Generated {} fragment alignments", fragPathAlignsById.size());

    spdlog::info("Assigning fragment origins");
    auto fragAssignment = getBestFragAssignment(genotypePaths, fragPathAlignsById);
    spdlog::info("Found assignments for {} frags", fragAssignment.fragIds.size());

    spdlog::info("Generating plot blueprint");
    auto lanePlots = generateBlueprint(genotypePaths, fragById, fragAssignment, fragPathAlignsById);

    spdlog::info("Writing SVG image to disk");
    generateSvg(lanePlots, svgPath);
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
        const string svgPath = args.outputPrefix + "." + locusId + ".svg";
        optional<string> phasingInfoPath;
        if (args.outputPhasingInfo)
        {
            phasingInfoPath = args.outputPrefix + "." + locusId + ".phasing.txt";
        }
        auto locusSpec = locusCatalog.at(locusId);
        visualizeLocus(args.referencePath, args.readsPath, args.vcfPath, locusId, locusSpec, svgPath, phasingInfoPath);
    }

    return 0;
}
