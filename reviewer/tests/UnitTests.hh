//
// Created by Clamons, Samuel on 5/3/22.
//
#pragma once
#include <optional>
#include <app/Workflow.hh>
#include <app/GenotypePaths.hh>
#include <app/Phasing.hh>
#include <app/Projection.hh>
#include <metrics/Metrics.hh>

#ifndef REVIEWER_UNITTESTS_HH
#define REVIEWER_UNITTESTS_HH

int three();

struct WorkflowResults {
    boost::optional<LocusSpecification> locusSpec;
    boost::optional<FragById> fragById;
    boost::optional<int> meanFragLen;
    boost::optional<std::vector<Diplotype>> pathsByDiplotype;
    boost::optional<ScoredDiplotypes> scoredDiplotypes;
    boost::optional<Diplotype> topDiplotype;
    boost::optional<PairPathAlignById> pairPathAlignById;
    boost::optional<FragPathAlignsById> fragPathAlignsById;
    boost::optional<FragAssignment> fragAssignment;
    boost::optional<MetricsByVariant> metricsByVariant;
};

WorkflowResults runToGetAligns(const WorkflowArguments& args);
WorkflowResults runToGetMeanFragLen(const WorkflowArguments& args);
WorkflowResults runToGetCandidateDiplotypes(const WorkflowArguments& args);
WorkflowResults runToGetTopDiplotype(const WorkflowArguments& args);
WorkflowResults runToProject(const WorkflowArguments& args);
WorkflowResults runToResolveByFragLen(WorkflowArguments& args);
WorkflowResults runToGetBestFragAssignment(WorkflowArguments& args);

vector<WorkflowArguments> exampleSampleArgs();

#endif //REVIEWER_UNITTESTS_HH
