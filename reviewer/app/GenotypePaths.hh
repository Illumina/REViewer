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

#include <string>
#include <vector>

#include "graphcore/Path.hh"

#include "core/LocusSpecification.hh"

using Diplotype = std::vector<graphtools::Path>;

/// Computes all possible diplotype paths at the given locus
/// \param meanFragLen: Mean fragment length
/// \param vcfPath: Path to the VCF file
/// \param locusSpec: Locus specification
/// \return Vector of all possible diplotype paths
///
/// Assumption:
/// All haplotype paths start at the first base of left flank
/// (node 0) and end at the last base of the right flank
/// (last node)
///
std::vector<Diplotype>
getCandidateDiplotypePaths(int meanFragLen, const std::string& vcfPath, const LocusSpecification& locusSpec);