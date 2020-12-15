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

#include "app/Projection.hh"

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

GenotypePaths phase(const PairGraphAlignById& fragGraphAlignById, const vector<GenotypePaths>& pathsByGenotype)
{
    int maxScore = 0;
    GenotypePaths bestGenotypePaths;
    for (const auto& genotypePaths : pathsByGenotype)
    {
        assert(genotypePaths.size() == 1 || genotypePaths.size() == 2);

        auto pairPathAlignById = project(genotypePaths, fragGraphAlignById);

        int genotypeScore = scorePath(0, pairPathAlignById);
        if (genotypePaths.size() == 2)
        {
            genotypeScore += scorePath(1, pairPathAlignById);
        }

        if (genotypeScore > maxScore)
        {
            maxScore = genotypeScore;
            bestGenotypePaths = genotypePaths;
        }
    }

    return bestGenotypePaths;
}
