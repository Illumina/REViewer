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

#include "app/Phasing.hh"

#include <algorithm>

#include "app/Projection.hh"

using boost::optional;
using std::pair;
using std::string;
using std::vector;

int scorePath(int pathIndex, const PairPathAlignById& pairPathAlignById)
{
    int pathScore = 0;
    for (const auto& idAndPairPathAlign : pairPathAlignById)
    {
        const auto& pairAlign = idAndPairPathAlign.second;
        for (const auto& readAlign : pairAlign.readAligns)
        {
            if (readAlign.pathIndex == pathIndex)
            {
                pathScore += score(*readAlign.align);
                break;
            }
        }

        for (const auto& mateAlign : pairAlign.mateAligns)
        {
            if (mateAlign.pathIndex == pathIndex)
            {
                pathScore += score(*mateAlign.align);
                break;
            }
        }
    }

    return pathScore;
}

ScoredDiplotypes scoreDiplotypes(const FragById& fragById, const vector<Diplotype>& diplotypes)
{
    vector<ScoredDiplotype> scoredDiplotypes;

    for (const auto& diplotype : diplotypes)
    {
        assert(diplotype.size() == 1 || diplotype.size() == 2);

        auto pairPathAlignById = project(diplotype, fragById);

        int genotypeScore = scorePath(0, pairPathAlignById);
        if (diplotype.size() == 2)
        {
            genotypeScore += scorePath(1, pairPathAlignById);
        }

        scoredDiplotypes.emplace_back(diplotype, genotypeScore);
    }

    std::sort(
        scoredDiplotypes.begin(), scoredDiplotypes.end(),
        [](const ScoredDiplotype& gt1, const ScoredDiplotype& gt2) { return gt1.second > gt2.second; });

    assert(!scoredDiplotypes.empty());
    return scoredDiplotypes;
}

std::ostream& operator<<(std::ostream& out, const ScoredDiplotype& diplotype) { return out; }
