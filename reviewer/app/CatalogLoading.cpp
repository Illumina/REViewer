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

#include "CatalogLoading.hh"

#include <algorithm>
#include <cassert>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include <boost/optional.hpp>

#include "spdlog/spdlog.h"
#include "thirdparty/json/json.hpp"

#include "app/LocusSpecDecoding.hh"
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

static bool checkIfFieldExists(const Json& record, const string& fieldName)
{
    return record.find(fieldName) != record.end();
}

static void assertFieldExists(const Json& record, const string& fieldName)
{
    if (!checkIfFieldExists(record, fieldName))
    {
        std::stringstream out;
        out << record;
        throw std::logic_error("Field " + fieldName + " must be present in " + out.str());
    }
}

static void assertRecordIsArray(const Json& record)
{
    if (!record.is_array())
    {
        std::stringstream out;
        out << record;
        throw std::logic_error("Expected array but got this instead " + out.str());
    }
}

static void makeArray(Json& record)
{
    if (record.type() != Json::value_t::array)
    {
        record = Json::array({ record });
    }
}

static VariantTypeFromUser decodeVariantTypeFromUser(const string& encoding)
{
    if (encoding == "RareRepeat")
    {
        return VariantTypeFromUser::kRareRepeat;
    }
    if (encoding == "Repeat")
    {
        return VariantTypeFromUser::kCommonRepeat;
    }
    if (encoding == "SmallVariant")
    {
        return VariantTypeFromUser::kSmallVariant;
    }
    if (encoding == "SMN")
    {
        return VariantTypeFromUser::kSMN;
    }
    else
    {
        throw std::logic_error("Encountered invalid variant type: " + encoding);
    }
}

static vector<string> generateIds(const std::string& locusId, const Json& variantRegionEncodings)
{
    if (variantRegionEncodings.size() == 1)
    {
        return { locusId };
    }

    vector<string> variantIds;
    for (const auto& regionEncoding : variantRegionEncodings)
    {
        variantIds.push_back(locusId + "_" + regionEncoding.get<string>());
    }

    return variantIds;
}

static LocusDescriptionFromUser loadUserDescription(Json& locusJson, const ReferenceContigInfo& contigInfo)
{
    LocusDescriptionFromUser userDescription;

    assertFieldExists(locusJson, "LocusId");
    userDescription.locusId = locusJson["LocusId"].get<string>();

    assertFieldExists(locusJson, "ReferenceRegion");
    makeArray(locusJson["ReferenceRegion"]);
    for (const auto& encoding : locusJson["ReferenceRegion"])
    {
        GenomicRegion region = decode(contigInfo, encoding.get<string>());
        userDescription.referenceRegions.push_back(region);
    }

    assertFieldExists(locusJson, "LocusStructure");
    userDescription.locusStructure = locusJson["LocusStructure"].get<string>();

    assertFieldExists(locusJson, "VariantType");
    makeArray(locusJson["VariantType"]);
    for (const auto& encoding : locusJson["VariantType"])
    {
        userDescription.variantTypesFromUser.push_back(decodeVariantTypeFromUser(encoding.get<string>()));
    }

    if (checkIfFieldExists(locusJson, "TargetRegion"))
    {
        makeArray(locusJson["TargetRegion"]);
        for (const auto& locusEncoding : locusJson["TargetRegion"])
        {
            GenomicRegion region = decode(contigInfo, locusEncoding.get<string>());
            userDescription.targetRegions.push_back(region);
        }
    }

    if (checkIfFieldExists(locusJson, "VariantId"))
    {
        makeArray(locusJson["VariantId"]);
        for (const auto& variantId : locusJson["VariantId"])
        {
            userDescription.variantIds.push_back(variantId.get<string>());
        }
    }
    else
    {
        userDescription.variantIds = generateIds(userDescription.locusId, locusJson["ReferenceRegion"]);
    }

    if (checkIfFieldExists(locusJson, "OfftargetRegions"))
    {
        assertRecordIsArray(locusJson["OfftargetRegions"]);
        for (const auto& locusEncoding : locusJson["OfftargetRegions"])
        {
            GenomicRegion region = decode(contigInfo, locusEncoding.get<string>());
            userDescription.offtargetRegions.push_back(region);
        }
    }

    if (checkIfFieldExists(locusJson, "ErrorRate"))
    {
        userDescription.errorRate = locusJson["ErrorRate"].get<double>();
    }
    if (checkIfFieldExists(locusJson, "LikelihoodRatioThreshold"))
    {
        userDescription.likelihoodRatioThreshold = locusJson["LikelihoodRatioThreshold"].get<double>();
    }
    if (checkIfFieldExists(locusJson, "MinimalLocusCoverage"))
    {
        userDescription.minLocusCoverage = locusJson["MinimalLocusCoverage"].get<double>();
    }

    return userDescription;
}

RegionCatalog loadLocusCatalogFromDisk(const string& catalogPath, const Reference& reference, int flankLength)
{
    std::ifstream inputStream(catalogPath.c_str());

    if (!inputStream.is_open())
    {
        throw std::runtime_error("Failed to open catalog file " + catalogPath);
    }

    Json catalogJson;
    inputStream >> catalogJson;
    makeArray(catalogJson);

    RegionCatalog catalog;
    for (auto& locusJson : catalogJson)
    {
        LocusDescriptionFromUser userDescription = loadUserDescription(locusJson, reference.contigInfo());
        try {
            LocusSpecification locusSpec = decodeLocusSpecification(userDescription, reference, flankLength);
            catalog.emplace(std::make_pair(locusSpec.locusId(), locusSpec));
        } catch (const std::exception& e) {
            std::cout << "Error on locus spec " + userDescription.locusId + ": " + e.what() << std::endl;
        }
    }

    return catalog;
}
