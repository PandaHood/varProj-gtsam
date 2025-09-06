//
// Created by jason on 2/2/25.
//
/* ----------------------------------------------------------------------------

 * GTSAM Copyright 2010-2019, Georgia Tech Research Corporation,
 * Atlanta, Georgia 30332-0415
 * All Rights Reserved
 * Authors: Frank Dellaert, et al. (see THANKS for the full author list)

 * See LICENSE for the license information

 * -------------------------------------------------------------------------- */

/**
 * @file LiftedRangeFactor.cpp
 * @date
 * @author Jason Xu, Nikolas Sanderson
 * @brief
 */

#include <gtsam/base/timing.h>
#include <gtsam/nonlinear/NonlinearFactorGraph.h>
#include <gtsam/slam/FrobeniusFactor.h>
#include "LiftedRangeFactor.h"
#include <cmath>
#include <iostream>
#include <vector>

using namespace std;
namespace gtsam {

//******************************************************************************
    template <size_t d>
    LiftedRangeFactor<d>::LiftedRangeFactor(Key j1, Key j2, Key r1 , const double  &RangeMeasure_ , size_t p,
                                            const SharedNoiseModel &model)
            : NoiseModelFactorN<LiftedPoseDP, Vector, UnitSphereD>(model, j1, j2, r1),
              measures_(RangeMeasure_),
              p_(p),
              d_(d) {
        if (noiseModel()->dim() != p_)
            throw std::invalid_argument(
                    "SEsyncFactor: model with incorrect dimension.");
    }

//******************************************************************************
    template <size_t d>
    void LiftedRangeFactor<d>::print(const std::string &s,
                                const KeyFormatter &keyFormatter) const {
        std::cout << s << "LiftedRangeFactor<" << p_ << " * " << d << ">(" << keyFormatter(key<1>()) << ","
                << keyFormatter(key<2>()) << "," << keyFormatter(key<3>()) << ")\n";
        traits<double>::Print(measures_, " Range measurement ");
        noiseModel_->print("  noise model: ");
    }

//******************************************************************************
    template <size_t d>
    bool LiftedRangeFactor<d>::equals(const NonlinearFactor &expected,
                                 double tol) const {
        auto e = dynamic_cast<const LiftedRangeFactor *>(&expected);
        return e != nullptr && NoiseModelFactorN<LiftedPoseDP, Vector, UnitSphereD>::equals(*e, tol) &&
               p_ == e->p_ && measures_ == e->measures_;
    }

//******************************************************************************
    template <size_t d>
    void LiftedRangeFactor<d>::fillJacobians(const LiftedPoseDP& P1, const Vector& L1, const UnitSphereD &R1, OptionalMatrixType H1,
                                             OptionalMatrixType H2, OptionalMatrixType H3) const {

        if (H1) {

            H1->resize(p_, LiftedPoseDP::Dimension(d, p_));
            H1->setZero();
            H1->block(0, StiefelManifoldKP::Dimension(d, p_), p_, p_) = Matrix::Identity(p_, p_);

        }
        if (H2) {
            H2->resize(p_, p_);
            *H2 = -Matrix::Identity(p_, p_);
        }
        if (H3) {
            H3->resize(p_, p_ - 1);
//            The same...
//            H3->resize(p_, LiftedPoseDP::Dimension(1, p_));
            // Need basis generator here
            *H3 = -measures_ * Matrix::Identity(p_, p_) * R1.G_;
        }

    }
//******************************************************************************
    template <size_t d>
    Vector LiftedRangeFactor<d>::evaluateError(const LiftedPoseDP& P1, const Vector& L1, const UnitSphereD &R1, OptionalMatrixType H1,
                                          OptionalMatrixType H2, OptionalMatrixType H3) const {

        const Vector T_P1 = P1.get_t();
        const Vector T_L1 = L1;

        if (T_P1.size() != static_cast<int>(p_) || T_L1.size() != static_cast<int>(p_))
            throw std::invalid_argument("Invalid dimension LiftedPoseDP or Vector values passed to "
                                        "LiftedRangeFactor<d>::evaluateError");

        this->fillJacobians(P1, L1, R1, H1, H2, H3);


        return T_P1 - T_L1 - measures_ * R1.matrix();
    }

    template <size_t d>
    void LiftedRangeFactor<d>::computeHessian(std::map<std::tuple<Key, Key,Key>, Matrix>& HMap) const
    {
        // sanity
        Eigen::Index D = 1;
        Eigen::MatrixXd I = Eigen::MatrixXd::Identity(D, D);
        Eigen::MatrixXd H = Eigen::MatrixXd::Zero(3*D, 3*D);

        // (t_i,t_i), (t_i,l_j), (l_j,t_i), (l_j,l_j)
        H.block(   0,    0, D, D) =  I;
        H.block(   0,    D, D, D) = -I;
        H.block(   D,    0, D, D) = -I;
        H.block(   D,    D, D, D) =  I;

        // cross‐terms with r
        H.block(   0, 2*D, D, D) = - measures_ * I;  // H_{t_i,r}
        H.block(   D, 2*D, D, D) = measures_ * I;  // H_{l_j,r}
        H.block(2*D,    0, D, D) =  -measures_ * I;  // H_{r,t_i}
        H.block(2*D,    D, D, D) = measures_ * I;  // H_{r,l_j}

        // (r,r) block
        H.block(2*D, 2*D, D, D) = measures_*measures_ * I;
        double sigma_cov = noiseModel_->sigmas()(0);
        double cov  = 1.0 / (1*sigma_cov*sigma_cov);

        // finally store half of the reordered Hessian
        if (HMap.count({key1(), key2(),key3()}) == 0)
            HMap[{key1(), key2(),key3()}] = cov * H;
        else
            HMap[{key1(), key2(),key3()}] += cov* H;
    }


//******************************************************************************
    template <size_t d>
    LiftedRangeFactorPosePose<d>::LiftedRangeFactorPosePose(Key j1, Key j2, Key r1 , const double  &RangeMeasure_ , size_t p,
                                            const SharedNoiseModel &model)
            : NoiseModelFactorN<LiftedPoseDP, LiftedPoseDP, UnitSphereD>(model, j1, j2, r1),
              measures_(RangeMeasure_),
              p_(p),
              d_(d) {
        if (noiseModel()->dim() != p_)
            throw std::invalid_argument(
                    "SEsyncFactor: model with incorrect dimension.");
    }

//******************************************************************************
    template <size_t d>
    void LiftedRangeFactorPosePose<d>::print(const std::string &s,
                                const KeyFormatter &keyFormatter) const {
        std::cout << s << "LiftedRangeFactor<" << p_ << " * " << d << ">(" << keyFormatter(key<1>()) << ","
                << keyFormatter(key<2>()) << "," << keyFormatter(key<3>()) << ")\n";
        traits<double>::Print(measures_, " Range measurement ");
        noiseModel_->print("  noise model: ");
    }

//******************************************************************************
    template <size_t d>
    bool LiftedRangeFactorPosePose<d>::equals(const NonlinearFactor &expected,
                                 double tol) const {
        auto e = dynamic_cast<const LiftedRangeFactorPosePose *>(&expected);
        return e != nullptr && NoiseModelFactorN<LiftedPoseDP, LiftedPoseDP, UnitSphereD>::equals(*e, tol) &&
               p_ == e->p_ && measures_ == e->measures_;
    }

//******************************************************************************
    template <size_t d>
    void LiftedRangeFactorPosePose<d>::fillJacobians(const LiftedPoseDP& P1, const LiftedPoseDP& L1, const UnitSphereD &R1, OptionalMatrixType H1,
                                             OptionalMatrixType H2, OptionalMatrixType H3) const {

        if (H1) {

            H1->resize(p_, LiftedPoseDP::Dimension(d, p_));
            H1->setZero();
            H1->block(0, StiefelManifoldKP::Dimension(d, p_), p_, p_) = Matrix::Identity(p_, p_);

        }
        if (H2) {
            H2->resize(p_, LiftedPoseDP::Dimension(d, p_));
            H2->setZero();
            // translation block: negative identity
            H2->block(0, StiefelManifoldKP::Dimension(d, p_), p_, p_) = -Matrix::Identity(p_, p_);
        }
        if (H3) {
            H3->resize(p_, p_ - 1);
//            The same...
//            H3->resize(p_, LiftedPoseDP::Dimension(1, p_));
            // Need basis generator here
            *H3 = -measures_ * Matrix::Identity(p_, p_) * R1.G_;
        }

    }
//******************************************************************************
    template <size_t d>
    Vector LiftedRangeFactorPosePose<d>::evaluateError(const LiftedPoseDP& P1, const LiftedPoseDP& L1, const UnitSphereD &R1, OptionalMatrixType H1,
                                          OptionalMatrixType H2, OptionalMatrixType H3) const {

        const Vector T_P1 = P1.get_t();
        const Vector T_L1 = L1.get_t();

        if (T_P1.size() != static_cast<int>(p_) || T_L1.size() != static_cast<int>(p_))
            throw std::invalid_argument("Invalid dimension LiftedPoseDP or Vector values passed to "
                                        "LiftedRangeFactor<d>::evaluateError");

        this->fillJacobians(P1, L1, R1, H1, H2, H3);


        return T_P1 - T_L1 - measures_ * R1.matrix();
    }
/* ************************************************************************* */

    //******************************************************************************
    template <size_t d>
    LiftedRangeFactorPointPoint<d>::LiftedRangeFactorPointPoint(Key j1, Key j2, Key r1 , const double  &RangeMeasure_ , size_t p,
                                            const SharedNoiseModel &model)
            : NoiseModelFactorN<Vector, Vector, UnitSphereD>(model, j1, j2, r1),
              measures_(RangeMeasure_),
              p_(p),
              d_(d) {
        if (noiseModel()->dim() != p_)
            throw std::invalid_argument(
                    "SEsyncFactor: model with incorrect dimension.");
    }

//******************************************************************************
    template <size_t d>
    void LiftedRangeFactorPointPoint<d>::print(const std::string &s,
                                const KeyFormatter &keyFormatter) const {
        std::cout << s << "LiftedRangeFactor<" << p_ << " * " << d << ">(" << keyFormatter(key<1>()) << ","
                << keyFormatter(key<2>()) << "," << keyFormatter(key<3>()) << ")\n";
        traits<double>::Print(measures_, " Range measurement ");
        noiseModel_->print("  noise model: ");
    }

//******************************************************************************
    template <size_t d>
    bool LiftedRangeFactorPointPoint<d>::equals(const NonlinearFactor &expected,
                                 double tol) const {
        auto e = dynamic_cast<const LiftedRangeFactorPointPoint *>(&expected);
        return e != nullptr && NoiseModelFactorN<Vector, Vector, UnitSphereD>::equals(*e, tol) &&
               p_ == e->p_ && measures_ == e->measures_;
    }

//******************************************************************************
    template <size_t d>
    void LiftedRangeFactorPointPoint<d>::fillJacobians(
     const Vector& P1, const Vector& L1, const UnitSphereD& R1,
     OptionalMatrixType H1, OptionalMatrixType H2, OptionalMatrixType H3) const {

        if (H1) {
            H1->resize(p_, p_);
            *H1 = Matrix::Identity(p_, p_);
        }
        if (H2) {
            H2->resize(p_, p_);
            *H2 = -Matrix::Identity(p_, p_);
        }
        if (H3) {
            H3->resize(p_, p_ - 1);
            *H3 = -measures_ * R1.G_;  // (R1.G_ is p x (p-1))
        }
    }
//******************************************************************************
    template <size_t d>
    Vector LiftedRangeFactorPointPoint<d>::evaluateError(const Vector& P1, const Vector& L1, const UnitSphereD &R1, OptionalMatrixType H1,
                                          OptionalMatrixType H2, OptionalMatrixType H3) const {

        const Vector T_P1 = P1;
        const Vector T_L1 = L1;

        if (T_P1.size() != static_cast<int>(p_) || T_L1.size() != static_cast<int>(p_))
            throw std::invalid_argument("Invalid dimension LiftedPoseDP or Vector values passed to "
                                        "LiftedRangeFactor<d>::evaluateError");

        this->fillJacobians(P1, L1, R1, H1, H2, H3);


        return T_P1 - T_L1 - measures_ * R1.matrix();
    }
// Explicit instantiation for d=2 and d=3
    template class LiftedRangeFactor<2>;
    template class LiftedRangeFactor<3>;
    template class LiftedRangeFactorPosePose<2>;
    template class LiftedRangeFactorPosePose<3>;
    template class LiftedRangeFactorPointPoint<2>;
    template class LiftedRangeFactorPointPoint<3>;


//******************************************************************************

} // namespace gtsam
