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

#include <map>
#include <vector>

#include "app/Projection.hh"

struct FragPathAlign
{
    FragPathAlign(ReadPathAlign readAlign, ReadPathAlign mateAlign)
        : readAlign(std::move(readAlign))
        , mateAlign(std::move(mateAlign))
    {
        assert(this->readAlign.pathIndex == this->mateAlign.pathIndex);
    }

    ReadPathAlign readAlign;
    ReadPathAlign mateAlign;
};

int getMeanFragLen(const PairGraphAlignById& pairGraphAlignById);

using FragPathAlignsById = std::map<std::string, std::vector<FragPathAlign>>;

FragPathAlignsById
resolveByFragLen(int meanFragLen, const GenotypePaths& paths, const PairPathAlignById& pairPathAlignById);
