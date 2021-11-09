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
#include <memory>
#include <string>

#include "graphalign/GraphAlignment.hh"

#include "app/LocusSpecification.hh"

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

FragById getAligns(const std::string& readsPath, const std::string& referencePath, const LocusSpecification& locusSpec);
