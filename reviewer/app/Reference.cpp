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

#include "Reference.hh"

#include <algorithm>
#include <memory>
#include <stdexcept>

using std::string;
using std::to_string;
using std::vector;

Reference::Reference(const string& referencePath)
    : referencePath_(referencePath)
    , contigInfo_({})
{
    htsFastaIndexPtr_ = fai_load(referencePath_.c_str());

    std::vector<std::pair<std::string, int64_t>> internalNamesAndSizes;

    for (int contigIndex = 0; contigIndex != faidx_nseq(htsFastaIndexPtr_); ++contigIndex)
    {
        const char* sequenceName = faidx_iseq(htsFastaIndexPtr_, contigIndex);
        int64_t sequenceLength = faidx_seq_len(htsFastaIndexPtr_, sequenceName);
        internalNamesAndSizes.emplace_back(sequenceName, sequenceLength);
    }

    contigInfo_ = ReferenceContigInfo(internalNamesAndSizes);
}

Reference::~Reference() { fai_destroy(htsFastaIndexPtr_); }

string Reference::getSequence(const string& contigName, int64_t start, int64_t end) const
{
    const int contigIndex = contigInfo_.getContigId(contigName);
    const char* contigNamePtr = faidx_iseq(htsFastaIndexPtr_, contigIndex);

    int extractedLength;
    // This htslib function is 0-based closed but our coordinates are half open
    char* sequencePtr = faidx_fetch_seq(htsFastaIndexPtr_, contigNamePtr, start, end - 1, &extractedLength);

    if (!sequencePtr || extractedLength < 0 || extractedLength < end - start)
    {
        const string encoding(contigName + ":" + to_string(start) + "-" + to_string(end));
        const string message = "Unable to extract " + encoding + " from " + referencePath_;
        throw std::runtime_error(message);
    }

    string sequence(sequencePtr);
    free(sequencePtr);
    std::transform(sequence.begin(), sequence.end(), sequence.begin(), ::toupper);

    return sequence;
}

string Reference::getSequence(const GenomicRegion& region) const
{
    return getSequence(contigInfo_.getContigName(region.contigIndex()), region.start(), region.end());
}
