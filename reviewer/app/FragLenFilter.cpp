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

#include "app/FragLenFilter.hh"

static int calcFragLen(const ReadPathAlign& readAlign, const ReadPathAlign& mateAlign)
{
    if (readAlign.pathIndex != mateAlign.pathIndex)
    {
        return std::numeric_limits<int>::max();
    }

    return std::max(readAlign.end, mateAlign.end) - std::min(readAlign.begin, mateAlign.begin);
}

FragPathAlignsById
resolveByFragLen(int meanFragLen, const GenotypePaths& paths, const PairPathAlignById& pairPathAlignById)
{
    FragPathAlignsById fragPathAlignsById;

    for (const auto& idAndPairPathAlign : pairPathAlignById)
    {
        const auto& fragId = idAndPairPathAlign.first;
        const PairPathAlign& pairPathAlign = idAndPairPathAlign.second;

        int bestFragLen = std::numeric_limits<int>::max();
        for (const auto& readAlign : pairPathAlign.readAligns)
        {
            for (const auto& mateAlign : pairPathAlign.mateAligns)
            {
                if (readAlign.pathIndex != mateAlign.pathIndex)
                {
                    continue;
                }

                const auto& path = paths[readAlign.pathIndex];
                const int fragLen = calcFragLen(readAlign, mateAlign);

                if (std::abs(fragLen - meanFragLen) < std::abs(bestFragLen - meanFragLen))
                {
                    bestFragLen = meanFragLen;
                    fragPathAlignsById[fragId].clear();
                }

                if (meanFragLen == bestFragLen)
                {
                    fragPathAlignsById[fragId].emplace_back(readAlign, mateAlign);
                }
            }
        }
    }

    return fragPathAlignsById;
}

int getMeanFragLen(const PairGraphAlignById& pairGraphAlignById)
{
    if (pairGraphAlignById.empty())
    {
        throw std::runtime_error("There are no read alignments in the target region");
    }

    int leftFlankId = 0;
    int rightFlankId = pairGraphAlignById.begin()->second.readAlign.path().graphRawPtr()->numNodes() - 1;

    double fragLenSum = 0;
    int numFlankingReads = 0;
    for (const auto& idAndPairGraphAlign : pairGraphAlignById)
    {
        const PairGraphAlign& pairAlign = idAndPairGraphAlign.second;

        const auto readStartNode = pairAlign.readAlign.path().getNodeIdByIndex(0);
        const auto mateStartNode = pairAlign.mateAlign.path().getNodeIdByIndex(0);

        const bool matesStartOnLeftFlank = readStartNode == leftFlankId && mateStartNode == leftFlankId;
        const bool matesStartOnRightFlank = readStartNode == rightFlankId && mateStartNode == rightFlankId;

        if (!matesStartOnLeftFlank && !matesStartOnRightFlank)
        {
            continue;
        }

        const int readStart = pairAlign.readAlign.path().startPosition();
        const int readEnd = readStart + static_cast<int>(pairAlign.readAlign.queryLength());

        const int mateStart = pairAlign.mateAlign.path().startPosition();
        const int mateEnd = mateStart + static_cast<int>(pairAlign.mateAlign.queryLength());

        fragLenSum += readEnd <= mateEnd ? mateEnd - readStart : readEnd - mateStart;
        ++numFlankingReads;
    }

    if (numFlankingReads == 0)
    {
        throw std::runtime_error("Unable to determine fragment length due to missing flanking read fragments");
    }

    return static_cast<int>(fragLenSum / numFlankingReads);
}