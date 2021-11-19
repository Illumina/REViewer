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

#include "LocusSpecification.hh"

#include <algorithm>
#include <cassert>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "spdlog/spdlog.h"
#include "thirdparty/json/json.hpp"

#include "core/Reference.hh"

using boost::optional;
using graphtools::NodeId;
using std::map;
using std::ostream;
using std::string;
using std::to_string;
using std::vector;

using Json = nlohmann::json;

namespace spd = spdlog;

LocusSpecification::LocusSpecification(RegionId locusId, graphtools::Graph regionGraph)
    : locusId_(std::move(locusId))
    , regionGraph_(std::move(regionGraph))
{
}

void LocusSpecification::addVariantSpecification(
    std::string id, VariantClassification classification, GenomicRegion referenceLocus, vector<NodeId> nodes,
    optional<NodeId> refNode)
{
    variantSpecs_.emplace_back(std::move(id), classification, std::move(referenceLocus), std::move(nodes), refNode);
}

const VariantSpecification& LocusSpecification::getVariantSpecById(const std::string& variantSpecId) const
{
    for (const auto& variantSpec : variantSpecs_)
    {
        if (variantSpec.id() == variantSpecId)
        {
            return variantSpec;
        }
    }

    throw std::logic_error("There is no variant " + variantSpecId + " in locus " + locusId_);
}

bool LocusSpecification::requiresGenomeWideDepth() const
{
    for (const auto& variantSpec : variantSpecs_)
    {
        if (variantSpec.classification().subtype == VariantSubtype::kSMN)
        {
            return true;
        }
    }

    return false;
}
