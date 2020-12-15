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

int runWorkflow(const WorkflowArguments& args)
{
    const int kFlankLength = 1000;
    Reference reference(args.referencePath);
    spdlog::info("Loading specification of locus {}", args.locusId);
    auto locusCatalog = loadLocusCatalogFromDisk(args.catalogPath, reference, kFlankLength);
    auto locusSpec = locusCatalog.at(args.locusId);

    ReadPairById fragById;
    auto fragGraphAlignById = getAligns(args.htsFilePath, args.referencePath, locusSpec, fragById);
    spdlog::info("Extracted {} frags", fragGraphAlignById.size());

    spdlog::info("Calculating fragment length");
    const int meanFragLen = getMeanFragLen(fragGraphAlignById);
    spdlog::info("Fragment length is estimated to be {}", meanFragLen);

    spdlog::info("Extracting genotype paths");
    auto pathsByGenotype = getCandidateGenotypePaths(meanFragLen, args.vcfPath, locusSpec);

    spdlog::info("Phasing");
    auto genotypePaths = phase(fragGraphAlignById, pathsByGenotype);
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
    generateSvg(lanePlots, args.outputPrefix + ".svg");

    return 0;
}
