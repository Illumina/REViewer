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

#pragma once

#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// Include the fai class from samtools
#include "htslib/faidx.h"

#include "GenomicRegion.hh"
#include "ReferenceContigInfo.hh"

class Reference
{
public:
    explicit Reference(const std::string& referencePath);
    ~Reference();

    std::string getSequence(const std::string& contigIndex, int64_t start, int64_t end) const;
    std::string getSequence(const GenomicRegion& region) const;

    const ReferenceContigInfo& contigInfo() const { return contigInfo_; }

private:
    std::string referencePath_;
    faidx_t* htsFastaIndexPtr_;
    ReferenceContigInfo contigInfo_;
};
