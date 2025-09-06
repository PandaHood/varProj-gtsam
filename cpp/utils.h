//
// Created by jason on 1/27/25.
//
/** This file provides a convenient set of utility functions for reading in a
set of pose-graph SLAM measurements and constructing the corresponding data
matrices used in the SE-Sync algorithm.
 *
 * Copyright (C) 2016 - 2022 by David M. Rosen (dmrosen@mit.edu)
 */
#ifndef STIEFELMANIFOLDEXAMPLE_UTILS_H
#define STIEFELMANIFOLDEXAMPLE_UTILS_H

/**
 * @file Utils.h
 * @date
 * @author Jason Xu, Nikolas Sanderson
 * @brief
 */
#pragma once
#include <string>
#include <Eigen/Sparse>
#include "RelativePoseMeasurement.h"
#include "types.h"
#include <chrono>
namespace gtsam {
    namespace DataParser {

/** Given the name of a file containing a description of a special Euclidean
 * synchronization problem expressed in the .g2o format (i.e. using "EDGE_SE2 or
 * EDGE_SE3:QUAT" measurements), this function constructs and returns the
 * corresponding vector of RelativePoseMeasurements, and reports the total
 * number of poses in the pose-graph */
        Measurement read_g2o_file(const std::string &filename, size_t &num_poses);


        Measurement read_pycfg_file(const std::string &filename);




        /** Given two matrices X, Y in SO(d)^n, this function computes and returns the
 * orbit distance d_S(X,Y) between them and (optionally) the optimal
 * registration G_S in SO(d) aligning Y to X, as described in Appendix C.1 of
 * the SE-Sync tech report.
 */
        Scalar dS(const Matrix &X, const Matrix &Y, Matrix *G_S = nullptr);

/** Given two matrices X, Y in O(d)^n, this function computes and returns the
 * orbit distance d_O(X,Y) between them and (optionally) the optimal
 * registration G_O in O(d) aligning Y to X, as described in Appendix C.1 of the
 * SE-Sync tech report.
 */
        Scalar dO(const Matrix &X, const Matrix &Y, Matrix *G_O = nullptr);

    } // namespace DataParser
} // namespace gtsam



/** This small file provides a pair of convenience functions that are useful for
 * measuring elapsed computation times.
 *
 * Copyright (C) 2017 - 2018 by David M. Rosen (dmrosen@mit.edu)
 */

namespace CFGStopwatch {

/** This function returns a chrono::time_point struct encoding the time at which
 * it is called.*/
    inline std::chrono::time_point<std::chrono::high_resolution_clock> tick() {
        return std::chrono::high_resolution_clock::now();
    }

/** When this function is called with a chrono::time_point struct returned by
 * tick(), it returns the elapsed time (in seconds) between the calls to tick()
 * and tock().*/
    inline double
    tock(const std::chrono::time_point<std::chrono::high_resolution_clock>
         &tick_time) {
        auto counter = std::chrono::high_resolution_clock::now() - tick_time;
        return std::chrono::duration_cast<std::chrono::milliseconds>(counter)
                       .count() /
               1000.0;
    }
} // namespace Stopwatch


#endif //STIEFELMANIFOLDEXAMPLE_UTILS_H
