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

#include "app/GenotypePaths.hh"

#include <cassert>
#include <fstream>
#include <map>
#include <utility>

#include <boost/algorithm/string.hpp>
#include <boost/optional.hpp>

using boost::optional;
using graphtools::Graph;
using graphtools::NodeId;
using graphtools::Path;
using std::map;
using std::pair;
using std::stoi;
using std::string;
using std::vector;

static vector<int> extractRepeatLengths(const string& vcfPath, const string& repeatId)
{
    std::ifstream vcfFile(vcfPath);
    if (vcfFile.is_open())
    {
        const string query = "VARID=" + repeatId + ";";
        string line;
        while (getline(vcfFile, line))
        {
            if (boost::find_first(line, query))
            {
                vector<string> pieces;
                boost::split(pieces, line, boost::is_any_of("\t"));
                const string sampleFields = pieces[pieces.size() - 1];
                boost::split(pieces, sampleFields, boost::is_any_of(":"));
                const string genotypeEncoding = pieces[2];

                boost::split(pieces, genotypeEncoding, boost::is_any_of("/"));
                vector<int> sizes;
                for (const auto& sizeEncoding : pieces)
                {
                    sizes.push_back(stoi(sizeEncoding));
                }
                return sizes;
            }
        }
        vcfFile.close();
    }
    else
    {
        throw std::runtime_error("Unable to open file " + vcfPath);
    }

    throw std::runtime_error("No VCF record for " + repeatId);
}

static vector<int> capLengths(int upperBound, const vector<int>& lengths)
{
    vector<int> cappedLength;
    cappedLength.reserve(lengths.size());
    for (int length : lengths)
    {
        cappedLength.push_back(length <= upperBound ? length : upperBound);
    }

    return cappedLength;
}

using NodeRange = pair<NodeId, NodeId>;
using Nodes = std::vector<graphtools::NodeId>;
using GenotypeNodes = std::vector<Nodes>;

map<NodeRange, GenotypeNodes>
getVariantPathPieces(int meanFragLen, const string& vcfPath, const LocusSpecification& locusSpec)
{
    map<NodeRange, GenotypeNodes> nodesByHaplotypeByVariant;
    for (const auto& variantSpec : locusSpec.variantSpecs())
    {
        assert(variantSpec.classification().type == VariantType::kRepeat);
        assert(variantSpec.nodes().size() == 1);
        NodeId repeatNode = variantSpec.nodes().front();
        vector<vector<NodeId>> variantPaths;

        auto repeatLens = extractRepeatLengths(vcfPath, variantSpec.id());
        repeatLens = capLengths(meanFragLen, repeatLens);

        variantPaths.reserve(repeatLens.size());
        for (int repeatLen : repeatLens)
        {
            variantPaths.emplace_back(repeatLen, repeatNode);
        }

        NodeId nodeRangeFrom = std::numeric_limits<NodeId>::max();
        NodeId nodeRangeTo = std::numeric_limits<NodeId>::lowest();

        for (NodeId node : variantSpec.nodes())
        {
            if (node != -1)
            {
                nodeRangeFrom = std::min(nodeRangeFrom, node);
                nodeRangeTo = std::max(nodeRangeTo, node);
            }
        }

        assert(variantPaths.size() == 2);
        nodesByHaplotypeByVariant.emplace(std::make_pair(nodeRangeFrom, nodeRangeTo), variantPaths);
    }

    return nodesByHaplotypeByVariant;
}

static optional<pair<GenotypeNodes, NodeId>>
getVariantGenotypeNodes(const map<NodeRange, GenotypeNodes>& nodeRangeToPaths, NodeId node)
{
    for (const auto& nodeRangeAndPaths : nodeRangeToPaths)
    {
        const auto& nodeRange = nodeRangeAndPaths.first;
        const auto& paths = nodeRangeAndPaths.second;
        if (nodeRange.first <= node && node <= nodeRange.second)
        {
            assert(paths.size() == 2);
            return std::make_pair(paths, nodeRange.second);
        }
    }

    return boost::none;
}

static vector<GenotypeNodes>
extendGenotypeNodes(const vector<GenotypeNodes>& genotypes, const GenotypeNodes& genotypeExtension)
{
    vector<GenotypeNodes> extendedGenotype;
    for (auto& genotype : genotypes)
    {
        assert(genotype.size() == genotypeExtension.size());
        if (genotype.size() == 1)
        {
            const Nodes& haplotypeExtension = genotypeExtension.front();
            Nodes extendedHaplotype = genotype.front();
            extendedHaplotype.insert(extendedHaplotype.end(), haplotypeExtension.begin(), haplotypeExtension.end());
            extendedGenotype.push_back({ extendedHaplotype });
        }
        else
        {
            assert(genotype.size() == 2);
            Nodes hap1Ext1 = genotype.front();
            hap1Ext1.insert(hap1Ext1.end(), genotypeExtension.front().begin(), genotypeExtension.front().end());

            Nodes hap2Ext2 = genotype.back();
            hap2Ext2.insert(hap2Ext2.end(), genotypeExtension.back().begin(), genotypeExtension.back().end());

            extendedGenotype.push_back({ hap1Ext1, hap2Ext2 });

            Nodes hap1Ext2 = genotype.front();
            hap1Ext2.insert(hap1Ext2.end(), genotypeExtension.back().begin(), genotypeExtension.back().end());

            Nodes hap2Ext1 = genotype.back();
            hap2Ext1.insert(hap2Ext1.end(), genotypeExtension.front().begin(), genotypeExtension.front().end());

            extendedGenotype.push_back({ hap1Ext2, hap2Ext1 });
        }
    }

    return extendedGenotype;
}

vector<GenotypePaths>
getCandidateGenotypePaths(int meanFragLen, const string& vcfPath, const LocusSpecification& locusSpec)
{
    map<NodeRange, GenotypeNodes> genotypeNodesByVariant = getVariantPathPieces(meanFragLen, vcfPath, locusSpec);

    // Assume that all variants have the same number of alleles
    const int numAlleles = genotypeNodesByVariant.begin()->second.size();

    vector<GenotypeNodes> nodesByGenotype = { GenotypeNodes(numAlleles, { 0 }) };

    NodeId node = 1;
    while (node != locusSpec.regionGraph().numNodes())
    {
        auto variantPathNodesAndLastNode = getVariantGenotypeNodes(genotypeNodesByVariant, node);
        if (variantPathNodesAndLastNode)
        {
            nodesByGenotype = extendGenotypeNodes(nodesByGenotype, variantPathNodesAndLastNode->first);
            node = variantPathNodesAndLastNode->second;
        }
        else
        {
            for (auto& genotypeNodes : nodesByGenotype)
            {
                for (auto& haplotypeNodes : genotypeNodes)
                {
                    haplotypeNodes.push_back(node);
                }
            }
        }

        ++node;
    }

    vector<GenotypePaths> pathsByGenotype;
    const NodeId rightFlankNode = locusSpec.regionGraph().numNodes() - 1;
    const int rightFlankLength = locusSpec.regionGraph().nodeSeq(rightFlankNode).length();
    for (const auto& genotypeNodes : nodesByGenotype)
    {
        GenotypePaths genotypePaths;
        for (const auto& haplotypeNodes : genotypeNodes)
        {
            genotypePaths.emplace_back(&locusSpec.regionGraph(), 0, haplotypeNodes, rightFlankLength);
        }
        pathsByGenotype.push_back(genotypePaths);
    }

    return pathsByGenotype;
}
