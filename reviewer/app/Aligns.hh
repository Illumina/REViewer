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

struct ReadPair
{
    ReadPair(std::string read, std::string mate)
        : read(std::move(read))
        , mate(std::move(mate))
    {
    }

    std::string read;
    std::string mate;
};

struct PairGraphAlign
{
    PairGraphAlign(GraphAlign readAlign, GraphAlign mateAlign)
        : readAlign(std::move(readAlign))
        , mateAlign(std::move(mateAlign))
    {
    }

    GraphAlign readAlign;
    GraphAlign mateAlign;
};

using ReadPairById = std::map<std::string, ReadPair>;
using PairGraphAlignById = std::map<std::string, PairGraphAlign>;

PairGraphAlignById getAligns(
    const std::string& htsFilePath, const std::string& referencePath, const LocusSpecification& locusSpec,
    ReadPairById& fragById);
