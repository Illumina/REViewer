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

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>
#include <tests/UnitTests.hh>
#include <app/Workflow.hh>
#include <app/FragLenFilter.hh>
#include <app/Origin.hh>
#include <boost/filesystem.hpp>
using std::string;
using std::vector;
using boost::filesystem::path;



// WORKFLOW UTILS
// These functions let us run the workflow until one of several different stopping
// points, grabbing all results up to that step.

WorkflowResults runToGetAligns(const WorkflowArguments& args)
{
    WorkflowResults results;
    Reference reference(args.referencePath);
    auto locusCatalog = loadLocusCatalogFromDisk(args.catalogPath, reference, args.locusExtensionLength);
    auto locusIds = getLocusIds(locusCatalog, args.locusId);
    auto phasingFile = initPhasingFile(args.outputPrefix);
    auto metricsFile = initMetricsFile(args.outputPrefix);

    results.locusSpec = locusCatalog.at(args.locusId);
    results.fragById = getAligns(args.readsPath, args.referencePath, results.locusSpec.value());
    return results;
}

WorkflowResults runToGetMeanFragLen(const WorkflowArguments& args)
{
    WorkflowResults results = runToGetAligns(args);
    results.meanFragLen = getMeanFragLen(results.fragById.value());
    return results;
}

WorkflowResults runToGetCandidateDiplotypes(const WorkflowArguments& args)
{
    WorkflowResults results = runToGetMeanFragLen(args);
    results.pathsByDiplotype = getCandidateDiplotypes(results.meanFragLen.value(), args.vcfPath, results.locusSpec.value());
    return results;
}

WorkflowResults runToGetTopDiplotype(const WorkflowArguments& args)
{
    WorkflowResults results = runToGetCandidateDiplotypes(args);
    results.scoredDiplotypes = scoreDiplotypes(results.fragById.value(), results.pathsByDiplotype.value());
    results.topDiplotype = results.scoredDiplotypes.value().front().first;
    return results;
}

WorkflowResults runToProject(const WorkflowArguments& args)
{
    WorkflowResults results = runToGetTopDiplotype(args);
    results.pairPathAlignById = project(results.topDiplotype.value(), results.fragById.value());
    return results;
}

WorkflowResults runToResolveByFragLen(WorkflowArguments& args)
{
    WorkflowResults results = runToProject(args);
    results.fragPathAlignsById = resolveByFragLen(results.meanFragLen.value(), results.topDiplotype.value(), results.pairPathAlignById.value());
    return results;
}

WorkflowResults runToGetBestFragAssignment(WorkflowArguments& args)
{
    WorkflowResults results = runToResolveByFragLen(args);
    results.fragAssignment = getBestFragAssignment(results.topDiplotype.value(), results.fragPathAlignsById.value());
    return results;
}
//END WORKFLOW UTILS

// WORKFLOW ARGUMENTS FOR STANDARD TEST INPUTS
vector<WorkflowArguments> exampleSampleArgs()
{
    string exampleSampleInfo[3][2] = {
            {"BEAN1", "HG00684"},
            {"BEAN1", "HG01682"},
            {"JPH3", "HG00277"}};

    vector<WorkflowArguments> sampleArgs(13);
    int n = 0;
    std::generate(sampleArgs.begin(), sampleArgs.end(), [&exampleSampleInfo, &n]() {
        WorkflowArguments sampleInfo;

        path testInputPath = path("install") / "tests" / "inputs";
        path testOutputPath = path("install") / "tests" / "outputs";

        string locusName  = exampleSampleInfo[n][0];
        string sampleName = exampleSampleInfo[n][1];

        std::ostringstream bamfile;
        bamfile << locusName << "_" << sampleName << ".bam";
        sampleInfo.readsPath = (testInputPath / "bamlets" / bamfile.str()).string();

        std::ostringstream vcffile;
        vcffile << sampleName << ".vcf";
        sampleInfo.vcfPath = (testInputPath / "vcfs" / vcffile.str()).string();

        sampleInfo.catalogPath = (testInputPath / "catalogs" / "stranger_variant_catalog_hg38_chr16.json").string();
        sampleInfo.referencePath = (testInputPath / "genomes" / "HG38_chr16.fa").string();

        sampleInfo.locusId = locusName;

        std::ostringstream outfile;
        outfile << locusName << "_" << sampleName;
        sampleInfo.outputPrefix = (testOutputPath / "images" / outfile.str()).string();
        sampleInfo.locusExtensionLength = 1000;
        return sampleInfo;
    });
    return sampleArgs;
}





