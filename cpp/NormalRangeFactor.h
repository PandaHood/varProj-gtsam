//
// Created by jason on 8/15/25.
//

#ifndef STIEFELMANIFOLDEXAMPLE_NORMALRANGEFACTOR_H
#define STIEFELMANIFOLDEXAMPLE_NORMALRANGEFACTOR_H

/**
 * @file   NormalRangeFactor.h
 * @date
 * @author  Jason
 * @brief  Non-Lifted Range Factor
 */
#pragma once

#include <gtsam/geometry/Rot2.h>
#include <gtsam/geometry/Rot3.h>
#include <gtsam/nonlinear/NonlinearFactor.h>
#include <gtsam/inference/Symbol.h>
#include "LiftedPose.h"
#include "UnitSphere.h"
#include <Eigen/Sparse>
#include <type_traits>
namespace gtsam {

/**
 *
 */
    template <size_t d>
    class GTSAM_EXPORT NormalRangeFactor : public NoiseModelFactorN<LiftedPoseDP, Vector> {
        double measures_; ///< dimensionality constants
        size_t p_, d_;               ///< dimensionality constants
//        std::shared_ptr<Matrix> G_; ///< matrix of vectorized tangent space basis generators

        // Select Rot2 or Rot3 interface based template parameter d
        using Rot = typename std::conditional<d == 2, Rot2, Rot3>::type;
        using Trans = typename std::conditional<d == 2, Vector2 , Vector3>::type;

    public:

        // Provide access to the Matrix& version of evaluateError:
        using NoiseModelFactorN<LiftedPoseDP, Vector>::evaluateError;

        /// @name Constructor
        /// @{

        /// Constructor. Residual is scalar => use a 1D noise model.
        NormalRangeFactor(Key j1, Key j2, const double  &RangeMeasure_ , size_t p,
                          const SharedNoiseModel &model = nullptr);

        /// @}
        /// @name Testable
        /// @{

        /// print with optional string
        void
        print(const std::string &s,
              const KeyFormatter &keyFormatter = DefaultKeyFormatter) const override;

        /// assert equality up to a tolerance
        bool equals(const NonlinearFactor &expected,
                    double tol = 1e-9) const override;

        /// @}
        /// @name NoiseModelFactorN methods
        /// @{

        ///
        Vector evaluateError(const LiftedPoseDP& P1, const Vector& L1, OptionalMatrixType H1,
                             OptionalMatrixType H2) const override;
        /// @}

//    private:
        /// Calculate Jacobians if asked, Only implemented for d=2 and 3 in .cpp
        void fillJacobians(const LiftedPoseDP& P1, const Vector& L1, OptionalMatrixType H1,
                           OptionalMatrixType H2) const;

    };

// Explicit instantiation for d=2 and d=3 in .cpp file:
    using NormalRangeFactor2 = NormalRangeFactor<2>;
    using NormalRangeFactor3 = NormalRangeFactor<3>;

} // namespace gtsam





#endif //STIEFELMANIFOLDEXAMPLE_NORMALRANGEFACTOR_H
