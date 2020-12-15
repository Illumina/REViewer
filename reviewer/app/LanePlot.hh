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

#include <list>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include <boost/optional.hpp>

#include "graphcore/Path.hh"

#include "app/Aligns.hh"
#include "app/GenomicRegion.hh"
#include "app/Origin.hh"

enum class FeatureType
{
    kRect,
    kRectWithLeftBreak,
    kRectWithRightBreak,
    kLine,
    kArrows,
    kVerticalLine
};

struct Feature
{
    Feature(FeatureType type, int length, std::string fill, std::string stroke)
        : type(type)
        , length(length)
        , fill(std::move(fill))
        , stroke(std::move(stroke))
    {
    }
    FeatureType type;
    int length;
    boost::optional<std::string> label;
    std::string fill;
    std::string stroke;
};

struct Segment
{
    Segment(int start, std::vector<Feature> features, double opacity)
        : start(start)
        , features(std::move(features))
        , opacity(opacity)
    {
        end = start;
        for (const auto& feature : this->features)
        {
            end += feature.length;
        }
    }

    int start;
    int end;
    std::vector<Feature> features;
    double opacity;
};

struct Lane
{
    Lane(int height, std::vector<Segment> segments)
        : height(height)
        , segments(std::move(segments))
    {
    }

    int height;
    std::vector<Segment> segments;
};

using LanePlot = std::vector<Lane>;

std::vector<LanePlot> generateBlueprint(
    std::vector<graphtools::Path> paths, const ReadPairById& fragById, const FragAssignment& fragAssignment,
    const FragPathAlignsById& fragPathAlignsById);
