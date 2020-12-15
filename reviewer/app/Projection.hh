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

#include <memory>
#include <vector>

#include <boost/optional.hpp>

#include "graphalign/GraphAlignment.hh"
#include "graphcore/Path.hh"

#include "app/Aligns.hh"
#include "app/GenotypePaths.hh"

struct ReadPathAlign
{
    ReadPathAlign(const graphtools::Path& hapPath, int pathIndex, int startIndexOnPath, GraphAlignPtr align);

    int pathIndex;
    int startIndexOnPath;
    int begin;
    int end;
    GraphAlignPtr align;
};

struct PairPathAlign
{
    PairPathAlign() = default;
    PairPathAlign(std::vector<ReadPathAlign> readAligns, std::vector<ReadPathAlign> mateAligns)
        : readAligns(std::move(readAligns))
        , mateAligns(std::move(mateAligns))
    {
    }

    std::vector<ReadPathAlign> readAligns;
    std::vector<ReadPathAlign> mateAligns;
};

int score(const graphtools::GraphAlignment& alignment, int matchScore = 5, int mismatchScore = -4, int gapScore = -8);

using PairPathAlignById = std::map<std::string, PairPathAlign>;
PairPathAlignById
project(const std::vector<graphtools::Path>& genotypePaths, const PairGraphAlignById& pairGraphAlignById);
