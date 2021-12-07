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

#include "graphalign/GraphAlignment.hh"

#include <cassert>

using GraphPath = graphtools::Path;
using GraphPaths = std::vector<GraphPath>;

using GraphAlign = graphtools::GraphAlignment;
using GraphAlignPtr = std::shared_ptr<GraphAlign>;

struct Read
{
    Read(std::string bases, std::string quals, GraphAlign align)
        : bases(std::move(bases))
        , quals(std::move(quals))
        , align(std::move(align))
    {
    }

    std::string bases;
    std::string quals;
    GraphAlign align;
};

struct Frag
{
    Frag(Read read, Read mate)
        : read(std::move(read))
        , mate(std::move(mate))
    {
    }

    Read read;
    Read mate;
};

using FragById = std::map<std::string, Frag>;

struct FragAssignment
{
    FragAssignment(std::vector<std::string> fragIds, std::vector<int> alignIndexByFrag)
        : fragIds(std::move(fragIds))
        , alignIndexByFrag(std::move(alignIndexByFrag))
    {
    }

    std::vector<std::string> fragIds;
    std::vector<int> alignIndexByFrag;
};

struct ReadPathAlign
{
    ReadPathAlign(const graphtools::Path& hapPath, int pathIndex, int startIndexOnPath, GraphAlignPtr align);

    int pathIndex;
    int startIndexOnPath;
    int begin;
    int end;
    GraphAlignPtr align;
};

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

using FragPathAlignsById = std::map<std::string, std::vector<FragPathAlign>>;
