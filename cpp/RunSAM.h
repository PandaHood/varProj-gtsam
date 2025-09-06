//
// Created by nikolas on 7/1/25.
//

#ifndef RUNSAM_H
#define RUNSAM_H
#include <queue>

#include "gtsam/geometry/Pose2.h"
#include <gtsam/geometry/Pose3.h>

#include "gtsam/inference/Symbol.h"
#include "gtsam/nonlinear/BatchFixedLagSmoother.h"
#include "gtsam/nonlinear/NonlinearFactorGraph.h"
#include "gtsam/nonlinear/Values.h"
#include "gtsam/slam/BetweenFactor.h"
//#include "gtsam_unstable/nonlinear/ConcurrentBatchFilter.h"
#include "gtsam_unstable/nonlinear/ConcurrentBatchSmoother.h"
#include <gtsam_unstable/nonlinear/IncrementalFixedLagSmoother.h>
#include <gtsam_unstable/nonlinear/ConcurrentIncrementalFilter.h>
#include "ConcurrentCertifiableBatchSmoother.h"
#include <gflags/gflags.h>
#include <glog/logging.h>

//#include "gtsam_unstable/nonlinear/ConcurrentBatchSmoother.h"

using namespace std;
using namespace gtsam;
using symbol_shorthand::X;

using gtsam::Values;
using gtsam::Pose2;
using gtsam::Pose3;
using gtsam::symbol_shorthand::X;
using Scalar = double;
namespace gtsam
{
    // Structure to store pending loop closures
    template<typename VALUE>
    struct PendingLoopClosure {
        Key keyFrom;
        Key keyTo;
        VALUE measurement;
        SharedNoiseModel noise;
        PendingLoopClosure(Key from, Key to, const VALUE& meas, const SharedNoiseModel& n)
                : keyFrom(from), keyTo(to), measurement(meas), noise(n) {}
    };

    struct SAMopts
    {
        bool verbosityLM = LevenbergMarquardtParams::SUMMARY;
        double absoluteErrorTol = 1e-10;
        double relativeErrorTol = 1e-10;
        size_t max_iter = 200;
        double lag = 10.0;


        double relinearizeThreshold = 0.01;
        double relinearizeSkip = 1; // ?
    };

    template<class Problem, size_t d>
    class RunSAM
    {
    private:
        std::unique_ptr<Problem> problem_;
    public:
        SAMopts opts_;
        RunSAM(){};

        Values refineWithFilterSmoother(
                ConcurrentIncrementalFilter& concurrentFilter,
                ConcurrentCertifiableBatchSmoother& concurrentSmoother,
                const NonlinearFactorGraph& loopClosures) {
            // 1) Separator-only summary
            NonlinearFactorGraph sepFactors;
            Values sepValues;
            concurrentFilter.getSummarizedFactors(sepFactors, sepValues);

            // 2) Full filter values
            Values allValues = concurrentFilter.getLinearizationPoint();

            // 3) Other keys = allValues \ sepValues
            KeySet otherKeys;
            for (Key k : allValues.keys()) {
                if (!sepValues.exists(k)) otherKeys.insert(k);
            }

            // 4) Extract other values
            Values otherValues;
            for (Key k : otherKeys) otherValues.insert(k, allValues.at(k));

            // 5) Extract alive filter factors on otherKeys
            NonlinearFactorGraph allFactors = concurrentFilter.getFactors();
            NonlinearFactorGraph otherFactors;
            for (const auto& fptr : allFactors) {
                if (!fptr) continue;
                bool onOther = true;
                for (Key k : fptr->keys()) {
                    if (!otherKeys.count(k)) { onOther = false; break; }
                }
                if (onOther) otherFactors.push_back(fptr);
            }

            // Use otherFactors as filter summary
            NonlinearFactorGraph filterSummaryFactors = otherFactors;
            Values filterSummaryValues = otherValues;

            // 6) Extract smoother factors and values
            NonlinearFactorGraph smootherFactors = concurrentSmoother.getFactors();
            Values smootherValues = concurrentSmoother.getLinearizationPoint();

            // 7) Combine graphs: filter summary, smoother, then loop closures
            NonlinearFactorGraph combinedGraph = filterSummaryFactors;
            combinedGraph.push_back(smootherFactors);
            combinedGraph.push_back(loopClosures);

            // 8) Check for key intersections
            if (!noKeyIntersection(filterSummaryValues, smootherValues)) {
                throw std::runtime_error("Overlapping keys between filter and smoother values");
            }

            // 9) Build initial estimate
            Values initialEstimate = smootherValues;
            for (const auto& kv : filterSummaryValues) {
                initialEstimate.insert(kv.key, kv.value);
            }

            // 10) Batch optimization
            LevenbergMarquardtParams params;
            LevenbergMarquardtOptimizer optimizer(combinedGraph, initialEstimate, params);
            // Try certifiable estimator
            auto problem = std::make_unique<Problem>(d, initialEstimate.size());

            //Values optimizedValues = problem->runIncrementalBatch(combinedGraph, initialEstimate, params);

            return optimizer.optimize();
            //return optimizedValues;

        }

        Values refineWithFilterSmoother_Local(
                ConcurrentIncrementalFilter& concurrentFilter,
                ConcurrentBatchSmoother& concurrentSmoother,
                const NonlinearFactorGraph& loopClosures) {
            // 1) Separator-only summary
            NonlinearFactorGraph sepFactors;
            Values sepValues;
            concurrentFilter.getSummarizedFactors(sepFactors, sepValues);

            // 2) Full filter values
            Values allValues = concurrentFilter.getLinearizationPoint();

            // 3) Other keys = allValues \ sepValues
            KeySet otherKeys;
            for (Key k : allValues.keys()) {
                if (!sepValues.exists(k)) otherKeys.insert(k);
            }

            // 4) Extract other values
            Values otherValues;
            for (Key k : otherKeys) otherValues.insert(k, allValues.at(k));

            // 5) Extract alive filter factors on otherKeys
            NonlinearFactorGraph allFactors = concurrentFilter.getFactors();
            NonlinearFactorGraph otherFactors;
            for (const auto& fptr : allFactors) {
                if (!fptr) continue;
                bool onOther = true;
                for (Key k : fptr->keys()) {
                    if (!otherKeys.count(k)) { onOther = false; break; }
                }
                if (onOther) otherFactors.push_back(fptr);
            }

            // Use otherFactors as filter summary
            NonlinearFactorGraph filterSummaryFactors = otherFactors;
            Values filterSummaryValues = otherValues;

            // 6) Extract smoother factors and values
            NonlinearFactorGraph smootherFactors = concurrentSmoother.getFactors();
            Values smootherValues = concurrentSmoother.getLinearizationPoint();

            // 7) Combine graphs: filter summary, smoother, then loop closures
            NonlinearFactorGraph combinedGraph = filterSummaryFactors;
            combinedGraph.push_back(smootherFactors);
            combinedGraph.push_back(loopClosures);

            // 8) Check for key intersections
            if (!noKeyIntersection(filterSummaryValues, smootherValues)) {
                throw std::runtime_error("Overlapping keys between filter and smoother values");
            }

            // 9) Build initial estimate
            Values initialEstimate = smootherValues;
            for (const auto& kv : filterSummaryValues) {
                initialEstimate.insert(kv.key, kv.value);
            }

            // 10) Batch optimization
            LevenbergMarquardtParams params;
            LevenbergMarquardtOptimizer optimizer(combinedGraph, initialEstimate, params);

            //auto problem = std::make_unique<Problem>(d, initialEstimate.size());

            //Values optimizedValues = problem->runIncrementalBatch(combinedGraph, initialEstimate, params);

                return optimizer.optimize();
            //return optimizedValues;
        }

        /*Exapnd to allow 3D*/
        template<typename VALUE>
        Values runSAM(string FLAGS_filepath, string output ="/home/jason/DPGO/StiefelManifold/results/Debug/") {
            // Parameters
            const double lag = opts_.lag;
            // Create the concurrent filter and smoother
            ISAM2Params isamParams;
            isamParams.relinearizeThreshold = opts_.relinearizeThreshold;
            isamParams.relinearizeSkip = opts_.relinearizeSkip;
            //    isamParams.factorization = ISAM2Params::Factorization::QR;
            ConcurrentIncrementalFilter concurrentFilter(isamParams);
            ConcurrentIncrementalFilter concurrentFilter_local(isamParams);

            //    ConcurrentBatchSmoother concurrentSmoother;
            LevenbergMarquardtParams params;
            if (opts_.verbosityLM)
            {
                params.verbosityLM = LevenbergMarquardtParams::SUMMARY;
            }
            params.absoluteErrorTol = opts_.absoluteErrorTol;
            params.relativeErrorTol = opts_.relativeErrorTol;
            params.maxIterations = opts_.max_iter;
            ConcurrentCertifiableBatchSmoother concurrentSmoother(params);
            ConcurrentBatchSmoother concurrentSmoother_local(params);

            // Containers for new data
            NonlinearFactorGraph newFactors;
            Values newValues;
            VALUE prev;
            // gtsam::Values type container for all the Filter's newest update
            Values filterResults, filterResults_local;
            // Queue for storing pending loop closures
            queue<PendingLoopClosure<VALUE>> pendingLoopClosures;

            // Initialize with prior

            /* make functions for 2D*/
            VALUE priorPose;
            SharedNoiseModel priorNoise;
            problem_->template setPrior<VALUE>(priorPose,priorNoise);

            newFactors.addPrior<VALUE>(X(0), priorPose, priorNoise);
            newValues.insert(X(0), priorPose);
            prev = priorPose;
            // Process initial update
            concurrentFilter.update(newFactors, newValues, FastList<Key>());
            synchronize(concurrentFilter, concurrentSmoother);
            concurrentFilter_local.update(newFactors, newValues, FastList<Key>());
            synchronize(concurrentFilter_local, concurrentSmoother_local);
            // Clear containers
            newFactors.resize(0);
            newValues.clear();
            // Main loop
            size_t index = 0;
            std::vector<VALUE> poseArray;
            std::pair<size_t, size_t> keys;

            // Map original g2o IDs to contiguous indices
            std::unordered_map<size_t, size_t> poseIndex;
            size_t nextPoseIdx = 1;
            poseIndex[0] = 0;

            std::ifstream infile(FLAGS_filepath);
            std::string line;

            /* Make this somehow applicable to 3D for example*/
            int id1, id2;
            VALUE Pose;
            std::shared_ptr<noiseModel::Diagonal> noise;
            int in = 0;
            while (std::getline(infile, line)) {
                /* Make this somehow applicable to 3D for example*/
                bool skip = problem_->template read_incremental_g2o<VALUE>(line, poseIndex, nextPoseIdx, Pose, noise, id1, id2);
                if (skip) continue;

                size_t i = poseIndex[id1], j = poseIndex[id2];
                if (id1 == id2 - 1) {  // Odometry measurement
                    VALUE guess = prev.compose(Pose);
                    prev = guess;

                    newValues.insert(X(j), guess);
                    newFactors.add(gtsam::BetweenFactor<VALUE>(
                            X(i), X(j), Pose, noise));

                    // Update smoother
                    // Determine keys to marginalize
                    FastList<Key> oldKeys;
                    if (id2 > lag) {
                        oldKeys.push_back(X(id2 - lag - 1));
                    }

                    // Update filter and smoothers
                    try {
                        concurrentFilter.update(newFactors, newValues, oldKeys);
                        concurrentFilter_local.update(newFactors, newValues, oldKeys);
                        for (size_t i = 0; i < 5; i++) {
                            concurrentFilter.update();
                        }
                        // export values of both filter
                        // X(j) is actually the new key for this update
                        gtsam::Pose2 latestPoseFilter =
                                concurrentFilter.calculateEstimate<gtsam::Pose2>(X(j));
                        filterResults.insert(X(j), latestPoseFilter);

                        gtsam::Pose2 latestPoseFilter_local =
                                concurrentFilter_local.calculateEstimate<gtsam::Pose2>(X(j));
                        filterResults_local.insert(X(j), latestPoseFilter_local);

                    } catch (const std::exception& e) {
                        // Handle update failure and print the factors and keys in the filter
                        LOG(FATAL) << "Error during filter update: " << e.what() << endl;
                        LOG(FATAL) << "Current factors in the filter:" << endl;

                        const NonlinearFactorGraph& factors = concurrentFilter.getFactors();
                        factors.print("Factors:");

                        LOG(FATAL) << "Current keys in the filter:" << endl;
                        const Values& values = concurrentFilter.getLinearizationPoint();
                        values.print("Keys:");
                        LOG(FATAL) << "Exiting due to filter update failure." << endl;

                        problem_->ExportValues(values, output + "debug_poses_filter", false);

                        auto smoother_values = concurrentSmoother.getLinearizationPoint();
                        problem_->ExportValues(smoother_values, output + "debug_poses_smoother", false);

                        exit(EXIT_FAILURE);

                    }
                    // Update filter and smoothers
                    newFactors.resize(0);
                    newValues.clear();

                    // Periodic synchronization and loop closure processing
                    if (id2 % 100 == 0) {
                        // First, check if we can process any pending loop closures
                        size_t n = pendingLoopClosures.size();
                        for (size_t i = 0; i < n; ++i) {
                            PendingLoopClosure loop = pendingLoopClosures.front();
                            pendingLoopClosures.pop();

                            const Values& smootherValues = concurrentSmoother.getLinearizationPoint();
                            if (smootherValues.exists(loop.keyFrom) && smootherValues.exists(loop.keyTo)) {
                                // both endpoints present → add to new factors
                                newFactors.push_back(
                                        BetweenFactor<VALUE>(
                                                loop.keyFrom, loop.keyTo, loop.measurement, loop.noise));
                            } else {
                                // one or both endpoints missing → defer for later
                                pendingLoopClosures.push(loop);
                            }
                        }

                        // Update smoother with accumulated loop closures if any
                        if (newFactors.size() > 0) {
                            concurrentSmoother.update<VALUE,d,Problem>(newFactors, Values());
                            concurrentSmoother.update<VALUE,d,Problem>();

                            concurrentSmoother_local.update(newFactors, Values());
                            concurrentSmoother_local.update();
                        }

                        for(size_t i =0; i < 5; i++) {
                            concurrentSmoother.update<VALUE,d,Problem>();
                            concurrentSmoother_local.update();
                        }

                        concurrentSmoother.update<VALUE,d,Problem>();
                        synchronize(concurrentFilter, concurrentSmoother);
                        problem_->ExportValues(concurrentSmoother.calculateEstimate(), output + "incremental_result_" + std::to_string(in), false);

                        concurrentSmoother_local.update();
                        synchronize(concurrentFilter_local, concurrentSmoother_local);
                        problem_->ExportValues(concurrentSmoother_local.calculateEstimate(), output + "incremental_result_local_"+ std::to_string(in), false);
                        in++;

                        // Print progress
                        LOG(INFO) << "Processed pose " << id2 << ", Pending loops: "
                             << pendingLoopClosures.size() << endl;
                    }
                    // Clear containers for next iteration
                    newFactors.resize(0);
                    newValues.clear();
                    index++;

                } else {  // Loop closure
                    // Create loop closure factor but store it for later
                    pendingLoopClosures.push(PendingLoopClosure(X(i), X(j), Pose, noise));
                }
            }

            // Need a final sync after reading all the data of filter
            concurrentSmoother.update<VALUE,d,Problem>();
            synchronize(concurrentFilter, concurrentSmoother);

            concurrentSmoother_local.update();
            synchronize(concurrentFilter_local, concurrentSmoother_local);

            for (size_t i = 0; i < 2; i++){
                concurrentSmoother.update<VALUE,d,Problem>();
            }
            for (size_t i = 0; i < 2; i++){
                concurrentSmoother_local.update();
            }
            auto result_certifiable = concurrentSmoother.calculateEstimate();
            problem_->ExportValues(result_certifiable, output + "incremental_result_Certifiable", false);

            auto result_local = concurrentSmoother_local.calculateEstimate();
            problem_->ExportValues(result_local, output + "incremental_result_Local", false);


            NonlinearFactorGraph loopClosureFactors = extractAllLoopClosures<VALUE>(pendingLoopClosures);
            Values finalEstimate = refineWithFilterSmoother(concurrentFilter, concurrentSmoother, loopClosureFactors);
            problem_->ExportValues(finalEstimate, output + "incremental_result_Final", false);
            Values finalEstimate_local = refineWithFilterSmoother_Local(concurrentFilter_local,concurrentSmoother_local ,loopClosureFactors);
            problem_->ExportValues(finalEstimate_local, output + "incremental_result_final_Local", false);

            // Export filter's "real time" updates
            problem_->ExportValues(filterResults, output + "CFS_results_certifiable", false);
            problem_->ExportValues(filterResults_local, output + "CFS_results_local", false);


            //Print final statistics
            LOG(INFO) << "\nFinal Statistics:" << endl;
            LOG(INFO) << "Processed poses: " << index << endl;
            LOG(INFO) << "Remaining loop closures: " << pendingLoopClosures.size() << endl;
            LOG(INFO) << "Filter keys: " << concurrentFilter.getLinearizationPoint().size() << endl;
            LOG(INFO) << "Smoother keys: " << concurrentSmoother.getLinearizationPoint().size() << endl;
            LOG(INFO) << "Final keys:" << finalEstimate.size() << endl;
            return finalEstimate;
        }

        template <typename VALUE>
        NonlinearFactorGraph extractAllLoopClosures(
            std::queue<PendingLoopClosure<VALUE>>& pendingLoopClosures) {
            NonlinearFactorGraph loopGraph;
            // Consume all closures in the queue
            while (!pendingLoopClosures.empty()) {
                PendingLoopClosure loop = pendingLoopClosures.front();
                pendingLoopClosures.pop();
                loopGraph.add(
                        BetweenFactor<VALUE>(
                                loop.keyFrom, loop.keyTo,
                                loop.measurement, loop.noise));
            }
            return loopGraph;
        }

        /**
         * @brief Check that two Values objects share no common keys.
         * @param v1 First Values map
         * @param v2 Second Values map
         * @return true if no intersecting keys; false otherwise
         */
        bool noKeyIntersection(const Values& v1, const Values& v2) {
            for (Key key : v1.keys()) {
                if (v2.exists(key)) return false;
            }
            return true;
        }
    };
}
#endif //RUNSAM_H
