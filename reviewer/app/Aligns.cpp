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

#include "app/Aligns.hh"

#include <vector>

#include <boost/algorithm/string.hpp>

#include "spdlog/spdlog.h"

#include "graphalign/GraphAlignmentOperations.hh"

extern "C"
{
#include "htslib/faidx.h"
#include "htslib/hts.h"
#include "htslib/sam.h"
}

#include "graphalign/GraphAlignment.hh"
#include "graphalign/GraphAlignmentOperations.hh"

using graphtools::decodeGraphAlignment;
using graphtools::GraphAlignment;
using std::map;
using std::string;
using std::vector;

FragById getAligns(const string& readsPath, const string& referencePath, const LocusSpecification& locusSpec)
{
    htsFile* htsFilePtr = nullptr;
    bam_hdr_t* htsHeaderPtr = nullptr;
    hts_idx_t* htsIndexPtr = nullptr;
    hts_itr_t* htsRegionPtr = nullptr;
    bam1_t* htsAlignmentPtr = nullptr;

    htsFilePtr = sam_open(readsPath.c_str(), "r");
    if (!htsFilePtr)
    {
        throw std::runtime_error("Failed to read BAM file " + readsPath);
    }

    // Required step for parsing of some CRAMs
    if (hts_set_fai_filename(htsFilePtr, referencePath.c_str()) != 0)
    {
        throw std::runtime_error("Failed to set index of: " + referencePath);
    }

    htsHeaderPtr = sam_hdr_read(htsFilePtr);
    if (!htsHeaderPtr)
    {
        throw std::runtime_error("Failed to read header of " + readsPath);
    }

    htsIndexPtr = sam_index_load(htsFilePtr, readsPath.c_str());
    if (!htsIndexPtr)
    {
        throw std::runtime_error("Failed to read index of " + readsPath);
    }

    const GenomicRegion region = locusSpec.variantSpecs().front().referenceLocus();

    // In case the reference is not fully compatible with the BAM
    faidx_t* referenceIndex = fai_load(referencePath.c_str());
    const char* contigName = faidx_iseq(referenceIndex, region.contigIndex());
    const int regionContigBamHeaderIndex = sam_hdr_name2tid(htsHeaderPtr, contigName);
    fai_destroy(referenceIndex);

    const auto& graph = locusSpec.regionGraph();
    const auto queryRegionStart = region.start() - graph.nodeSeq(0).length();
    const auto queryRegionEnd = region.end() + graph.nodeSeq(graph.numNodes() - 1).length();
    htsRegionPtr = sam_itr_queryi(htsIndexPtr, regionContigBamHeaderIndex, queryRegionStart, queryRegionEnd);
    if (htsRegionPtr == nullptr)
    {
        throw std::runtime_error("Failed to extract reads from the specified region");
    }

    FragById fragById;
    std::map<string, Read> unpairedCache;
    htsAlignmentPtr = bam_init1();
    while (sam_itr_next(htsFilePtr, htsRegionPtr, htsAlignmentPtr) >= 0)
    {
        // Decode read sequence
        string bases;
        uint8_t* htsSeqPtr = bam_get_seq(htsAlignmentPtr);
        const int readLength = htsAlignmentPtr->core.l_qseq;
        bases.resize(readLength);
        for (int index = 0; index != readLength; ++index)
        {
            bases[index] = seq_nt16_str[bam_seqi(htsSeqPtr, index)];
        }

        // Decode base call quality sequence
        string quals;
        uint8_t* htsQualsPtr = bam_get_qual(htsAlignmentPtr);
        quals.resize(readLength);

        for (int index = 0; index != readLength; ++index)
        {
            quals[index] = static_cast<char>(33 + htsQualsPtr[index]);
        }

        // Decode graph alignment
        if (!bam_get_l_aux(htsAlignmentPtr))
        {
            throw std::runtime_error("All BAM alignments are required to have \"XG\" auxiliary tag");
        }

        int auxPos = 0;
        uint8_t* aux = bam_get_aux(htsAlignmentPtr);
        const string tag(aux + auxPos, aux + auxPos + 2);
        const string tagType(aux + auxPos + 2, aux + auxPos + 3);
        auxPos += 3;
        if (tag != "XG" || tagType != "Z")
        {
            throw std::runtime_error("Unexpected auxiliary tag " + tag + ":" + tagType);
        }

        int valueLen = 0;
        while (*(aux + auxPos + valueLen) != '\0')
        {
            ++valueLen;
        }

        string fragmentId = bam_get_qname(htsAlignmentPtr);
        const string cigarEncoding(aux + auxPos, aux + auxPos + valueLen);
        vector<string> pieces;
        boost::split(pieces, cigarEncoding, boost::is_any_of(","));
        assert(pieces.size() == 3);
        const string& locusId = pieces[0];
        int pos = std::stoi(pieces[1]);
        const string& cigar = pieces[2];

        if (locusId != locusSpec.locusId())
        {
            continue;
        }

        GraphAlignment align = decodeGraphAlignment(pos, cigar, &locusSpec.regionGraph());
        if (!graphtools::checkConsistency(align, bases))
        {
            spdlog::warn("Encountered inconsistent alignment \n{}", prettyPrint(align, bases));
        }

        // Store read and alignment
        if (unpairedCache.find(fragmentId) == unpairedCache.end())
        {
            Read read(std::move(bases), std::move(quals), align);
            unpairedCache.emplace(fragmentId, read);
        }
        else
        {
            auto mate = std::move(unpairedCache.at(fragmentId));
            unpairedCache.erase(fragmentId);

            Read read(bases, quals, align);
            fragById.emplace(fragmentId, Frag(read, mate));
        }
    }

    if (!unpairedCache.empty())
    {
        spdlog::warn("Found {} unpaired reads", unpairedCache.size());
    }

    bam_destroy1(htsAlignmentPtr);
    htsAlignmentPtr = nullptr;

    hts_itr_destroy(htsRegionPtr);
    htsRegionPtr = nullptr;

    hts_idx_destroy(htsIndexPtr);
    htsIndexPtr = nullptr;

    bam_hdr_destroy(htsHeaderPtr);
    htsHeaderPtr = nullptr;

    sam_close(htsFilePtr);
    htsFilePtr = nullptr;

    return fragById;
}
