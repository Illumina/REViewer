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

#include "tests/UnitTests.hh"
#include "snps/Workflow.hh"

#include <catch2/catch.hpp>
#include <boost/filesystem.hpp>

extern "C"
{
#include "htslib/faidx.h"
#include "htslib/hts.h"
#include "htslib/sam.h"
}

using namespace snps;

TEST_CASE("Initializing SNP calling workflow", "[SNP calling]")
{
    GraphPaths paths;
    FragById fragById;
    FragAssignment fragAssignment({}, {});
    FragPathAlignsById fragPathAlignsById;
    callSnps(paths, fragById, fragAssignment, fragPathAlignsById);
//    REQUIRE(snpCalls == SnpCalls());
}

TEST_CASE("Checking that example bam files can be opened", "[bam_open]")
{
    htsFile* htsFilePtr = nullptr;
    bam_hdr_t* htsHeaderPtr = nullptr;
    hts_idx_t* htsIndexPtr = nullptr;
    hts_itr_t* htsRegionPtr = nullptr;
    bam1_t* htsAlignmentPtr = nullptr;
    for (auto args : exampleSampleArgs())
    {
        std::cout << "Attempting to read bam at " << args.readsPath << std::endl;

        htsFilePtr = sam_open(args.readsPath.c_str(), "r");
        if (!htsFilePtr)
        {
        throw std::runtime_error("Failed to read BAM file " + args.readsPath);
        }

        // Required step for parsing of some CRAMs
        if (hts_set_fai_filename(htsFilePtr, args.referencePath.c_str()) != 0)
        {
        throw std::runtime_error("Failed to set index of: " + args.referencePath);
        }

        htsHeaderPtr = sam_hdr_read(htsFilePtr);
        if (!htsHeaderPtr)
        {
        throw std::runtime_error("Failed to read header of " + args.readsPath);
        }

        htsIndexPtr = sam_index_load(htsFilePtr, args.readsPath.c_str());
        if (!htsIndexPtr)
        {
        throw std::runtime_error("Failed to read index of " + args.readsPath);
        }
    }
}

TEST_CASE("Check that consensus pileups match by-eye expectations", "[pileup_consensus]")
{
    std::map<string,string> expectedPileups = {
            {"BEAN1_HG00684.bam", "TAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAATAAAATAAAA"},
            {"BEAN1_HG01682.bam", "TAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAATAAAATAAAA"},
            {"JPH3_HG00277.bam", "CTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTG"}
    };
    for (auto args : exampleSampleArgs())
    {
        string bamfile = args.readsPath;
        int n = args.readsPath.length();
        for(int i = n-1; i >= 0; i--)
        {
            if (args.readsPath[i] == '/')
            {
                bamfile = args.readsPath.substr(i+1, n - i - 1);
                break;
            }
        }
        auto results = runToGetBestFragAssignment(args);
        auto consensusString = pileupConsensus(results.topDiplotype.value(), results.fragById.value(),
                                               results.fragAssignment.value(), results.fragPathAlignsById.value());

        std::cout << "Expected from " << bamfile << ": " << expectedPileups[bamfile] << "\n\t got " << consensusString << std::endl;
        REQUIRE(consensusString == expectedPileups[bamfile]);
    }
}