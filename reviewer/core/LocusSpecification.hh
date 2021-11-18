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

#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <boost/optional.hpp>

#include "graphcore/Graph.hh"
#include "thirdparty/json/json.hpp"

#include "core/GenomicRegion.hh"
#include "core/Reference.hh"
#include "core/VariantSpecification.hh"

using RegionId = std::string;
using NodeToRegionAssociation = std::unordered_map<graphtools::NodeId, GenomicRegion>;

class LocusSpecification
{
public:
    LocusSpecification(RegionId locusId, graphtools::Graph regionGraph);

    const RegionId& locusId() const { return locusId_; }
    const graphtools::Graph& regionGraph() const { return regionGraph_; }
    const std::vector<VariantSpecification>& variantSpecs() const { return variantSpecs_; }
    void addVariantSpecification(
        std::string id, VariantClassification classification, GenomicRegion referenceLocus,
        std::vector<graphtools::NodeId> nodes, boost::optional<graphtools::NodeId> optionalRefNode);
    const VariantSpecification& getVariantSpecById(const std::string& variantSpecId) const;

    bool requiresGenomeWideDepth() const;

private:
    std::string locusId_;
    graphtools::Graph regionGraph_;
    std::vector<VariantSpecification> variantSpecs_;
};

using RegionCatalog = std::map<RegionId, LocusSpecification>;
