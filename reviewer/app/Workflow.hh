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
#include <string>

#include "app/CatalogLoading.hh"
using std::string;
using std::vector;

struct WorkflowArguments
{
    string readsPath;
    string vcfPath;
    string catalogPath;
    string referencePath;
    string locusId;
    string outputPrefix;
    int locusExtensionLength;
};

int runWorkflow(const WorkflowArguments& args);
std::ofstream initPhasingFile(const string& prefix);
std::ofstream initMetricsFile(const string& prefix);
vector<string> getLocusIds(const RegionCatalog& catalog, const string& encoding);
