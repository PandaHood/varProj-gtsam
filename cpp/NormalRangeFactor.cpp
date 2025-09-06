//
// Created by jason on 8/15/25.
//
/* ----------------------------------------------------------------------------

 * GTSAM Copyright 2010-2019, Georgia Tech Research Corporation,
 * Atlanta, Georgia 30332-0415
 * All Rights Reserved
 * Authors: Frank Dellaert, et al. (see THANKS for the full author list)

 * See LICENSE for the license information

 * -------------------------------------------------------------------------- */

/**
 * @file NormalRangeFactor.cpp
 * @date
 * @author Jason
 * @brief  Non-lifted range factor
 */

#include <gtsam/base/timing.h>
#include <gtsam/nonlinear/NonlinearFactorGraph.h>
#include <gtsam/slam/FrobeniusFactor.h>
#include "NormalRangeFactor.h"
#include <cmath>
#include <iostream>
#include <vector>

using namespace std;
namespace gtsam {

//******************************************************************************
    template <size_t d>
    NormalRangeFactor<d>::NormalRangeFactor(Key j1, Key j2 , const double  &RangeMeasure_ , size_t p,
                                            const SharedNoiseModel &model)
            : NoiseModelFactorN<LiftedPoseDP, Vector>(model, j1, j2),
              measures_(RangeMeasure_),
              p_(p),
              d_(d) {
        // Residual is scalar => noise model dim must be 1
        if (noiseModel()->dim() != 1)
            throw std::invalid_argument(
                    "NormalRangeFactor: noise model must be 1D for a range residual.");
    }

//******************************************************************************
    template <size_t d>
    void NormalRangeFactor<d>::print(const std::string &s,
                                     const KeyFormatter &keyFormatter) const {
        std::cout << s << "NormalRangeFactor<" << p_ << " * " << d << ">(" << keyFormatter(key<1>()) << ","
                  << keyFormatter(key<2>()) << ")\n";
        traits<double>::Print(measures_, " Range measurement ");
        noiseModel_->print("  noise model: ");
    }

//******************************************************************************
    template <size_t d>
    bool NormalRangeFactor<d>::equals(const NonlinearFactor &expected,
                                      double tol) const {
        auto e = dynamic_cast<const NormalRangeFactor *>(&expected);
        return e != nullptr && NoiseModelFactorN<LiftedPoseDP, Vector>::equals(*e, tol) &&
               p_ == e->p_ && measures_ == e->measures_;
    }

//******************************************************************************
    template <size_t d>
    void NormalRangeFactor<d>::fillJacobians(const LiftedPoseDP& P1, const Vector& L1, OptionalMatrixType H1,
                                             OptionalMatrixType H2) const {
        // Compute unit direction for proper Jacobians
        const Vector tP = P1.get_t();
        const Vector tL = L1;
        Vector diff = tP - tL;
        const double dist = diff.norm();

        // Handle near-zero range robustly
        Vector grad = Vector::Zero(p_);
        if (dist > 1e-12) grad = diff.transpose() / dist;  // 1 x p_

        if (H1) {
            // 1 x total-dim of LiftedPoseDP (rotation-block then translation-block)
            H1->resize(1, LiftedPoseDP::Dimension(d, p_));
            H1->setZero();
            // place gradient in the translation block of the pose
            const size_t transOffset = StiefelManifoldKP::Dimension(d, p_);
            H1->block(0, transOffset, 1, p_) = grad.transpose();
        }
        if (H2) {
            // 1 x p_ for the lifted vector
            H2->resize(1, p_);
            *H2 = -grad.transpose();
        }
    }

//******************************************************************************
    template <size_t d>
    Vector NormalRangeFactor<d>::evaluateError(const LiftedPoseDP& P1, const Vector& L1, OptionalMatrixType H1,
                                               OptionalMatrixType H2) const {

        const Vector T_P1 = P1.get_t();
        const Vector T_L1 = L1;

        if (T_P1.size() != static_cast<int>(p_) || T_L1.size() != static_cast<int>(p_))
            throw std::invalid_argument("Invalid dimension LiftedPoseDP or Vector values passed to "
                                        "NormalRangeFactor<d>::evaluateError");

        this->fillJacobians(P1, L1, H1, H2);

        Vector errorVector(1);
        errorVector[0] = (T_P1 - T_L1).norm() - measures_;
        return errorVector;
    }


/* ************************************************************************* */
// Explicit instantiation for d=2 and d=3
    template class NormalRangeFactor<2>;
    template class NormalRangeFactor<3>;

//******************************************************************************

} // namespace gtsam
