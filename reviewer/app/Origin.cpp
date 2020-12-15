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

#include "app/Origin.hh"

using graphtools::NodeId;
using graphtools::Path;
using std::string;
using std::vector;

FragAssignment getBestFragAssignment(const vector<Path>& hapPaths, const FragPathAlignsById& fragPathAlignsById)
{
    vector<string> fragIds;
    fragIds.reserve(fragPathAlignsById.size());
    for (const auto& idAndFragAligns : fragPathAlignsById)
    {
        fragIds.push_back(idAndFragAligns.first);
    }

    vector<int> alignIndexByFrag(fragIds.size(), 0);
    for (int fragIndex = 0; fragIndex != fragIds.size(); ++fragIndex)
    {
        const auto& fragId = fragIds[fragIndex];
        int numOrigins = fragPathAlignsById.at(fragId).size();
        int originIndex = rand() % numOrigins;
        alignIndexByFrag[fragIndex] = originIndex;
    }

    return { fragIds, alignIndexByFrag };
}
