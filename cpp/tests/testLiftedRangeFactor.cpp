/* ----------------------------------------------------------------------------
 * Copyright 2025, Northeastern University Robust Autonomy Lab, * Boston, MA 02139
 * All Rights Reserved
 * Authors: Zhexin Xu, Nikolas Sanderson, et al. (see README for the full author list)
 * See LICENSE for the license information
 * -------------------------------------------------------------------------- */

#include "../StiefelManifold.h"
#include "../StiefelManifold-inl.h"
#include "../SEsyncFactor.h"
#include "../LiftedRangeFactor.h"
#include <iostream>
#include <gtsam/base/Matrix.h>
#include <gtsam/nonlinear/factorTesting.h>
#include <CppUnitLite/TestHarness.h>
#include <gtsam/base/testLie.h>
#include <gtsam/geometry/SO4.h>
#include <gtsam/geometry/Pose2.h>

using namespace gtsam;


// Test Jacobian EvaluateError
TEST(LiftedRangeFactor, EvaluateError) {
    int p = 3, d = 2;
    const std::default_random_engine::result_type seed = std::default_random_engine::default_seed;
    StiefelManifoldKP Q1 = StiefelManifoldKP::Random(seed, 2, 3);
    Q1.initializeTangentSpace();
    UnitSphereD R1 = UnitSphereD::Random(seed, 1, 3);
    R1.initializeTangentSpace();
    Vector t1(2), t2(2);
    t1 << 0.05, 0.05;
    t2 << 0.1, 0.2;
    double range_measurement = (t2 - t1).norm();
    auto model = noiseModel::Isotropic::Sigma(p, 1.2);
    auto t12 = Vector2(0.095, 0.195);
    LiftedRangeFactor2 factor1(1, 2, 3, range_measurement, 3, model);
    Matrix H1, H2, H3;
    LiftedPoseDP P1(Q1, LiftedPoseDP::LiftToRp(t1, p));
    Vector P2 = LiftedPoseDP::LiftToRp(t2, p);
    Vector error = factor1.evaluateError(P1, P2, R1, H1, H2, H3);
    EXPECT_LONGS_EQUAL(error.size(), p);
}

// Test Jacobian with GTSAM's numerical Jacobian TEST MACRO, d = 2
TEST(LiftedRangeFactor, NumericalJacobain2) {
    const std::default_random_engine::result_type seed = std::default_random_engine::default_seed;
    size_t d = 2;
    const Rot2 R1(0.3);
    Vector t1(2), t2(2), t12(2);
    t1 << 1, 1;
    t2 << 6, 7;
    double range_measurement = (t2 - t1).norm();
    for (const int p : {5, 4, 3, 2}) {
        auto model = noiseModel::Isotropic::Sigma(p, 1.2);
        Matrix M = Matrix::Zero(p, d);
        M.block(0, 0, 2, 2) = R1.matrix();
        StiefelManifoldKP Q1(M);
        Q1.initializeTangentSpace();
        LiftedPoseDP P1(Q1, LiftedPoseDP::LiftToRp(t1, p));
        UnitSphereD R1 = UnitSphereD::Random(seed, 1, p);
        R1.initializeTangentSpace();
        Vector P2 = LiftedPoseDP::LiftToRp(t2, p);
        auto factor = LiftedRangeFactor2 (1, 2, 3, range_measurement, p, model);
        Matrix H1, H2, H3;
        factor.evaluateError(P1, P2, R1, H1, H2);
        Values values;
        values.insert(1, P1);
        values.insert(2, P2);
        values.insert(3, R1);
        EXPECT_CORRECT_FACTOR_JACOBIANS(factor, values, 1e-7, 1e-5);
    }
}

// Just using those rotation expression from unit test of SOn
namespace so3 {
    SO3 id;
    Vector3 v1 = (Vector(3) << 0.1, 0, 0).finished();
    SO3 R1 = SO3::Expmap(v1);
    Vector3 v2 = (Vector(3) << 0.01, 0.02, 0.03).finished();
    SO3 R2 = SO3::Expmap(v2);
    SO3 R12 = R1.between(R2);
} // namespace so3

//******************************************************************************
namespace submanifold {
    SO4 id;
    Vector6 v1 = (Vector(6) << 0, 0, 0, 0.1, 0, 0).finished();
    SO3 R1 = SO3::Expmap(v1.tail<3>());
    SO4 Q1 = SO4::Expmap(v1);
    Vector6 v2 = (Vector(6) << 0, 0, 0, 0.01, 0.02, 0.03).finished();
    SO3 R2 = SO3::Expmap(v2.tail<3>());
    SO4 Q2 = SO4::Expmap(v2);
    SO3 R12 = R1.between(R2);
} // namespace submanifold

// Test Jacobian with GTSAM's numerical Jacobian TEST MACRO, d = 3
TEST(LiftedRangeFactor, NumericalJacobain3) {
    const std::default_random_engine::result_type seed = std::default_random_engine::default_seed;
    size_t d = 3;
    Vector t1(3), t2(3), t12(3);
    t1 << 0, 0, 0.1;
    t12 << 0.6, 0.7, 0.7;
    t2 = submanifold::R2.matrix().inverse() * t12;
    double range_measurement = t12.norm();
    for (const size_t p : {5, 4, 3}) {
        auto model = noiseModel::Isotropic::Sigma(p , 1.2);
        Matrix M = Matrix::Zero(p, d);
        M.topLeftCorner(3, 3) = submanifold::R1.matrix();
        StiefelManifoldKP Q1(M);
        Q1.initializeTangentSpace();
        LiftedPoseDP P1(Q1, LiftedPoseDP::LiftToRp(t1, p));
        UnitSphereD R1 = UnitSphereD::Random(seed, 1, p);
        R1.initializeTangentSpace();
        Vector P2 = LiftedPoseDP::LiftToRp(t2, p);
        auto factor = LiftedRangeFactor3 (1, 2, 3, range_measurement, p, model);
        Matrix H1, H2, H3;
        factor.evaluateError(P1, P2, R1, H1, H2, H3);
        Values values;
        values.insert(1, P1);
        values.insert(2, P2);
        values.insert(3, R1);
        EXPECT_CORRECT_FACTOR_JACOBIANS(factor, values, 1e-7, 1e-5);
    }
}

/* ************************************************************************* */
int main() {
    TestResult tr;
    return TestRegistry::runAllTests(tr);
}
/* ************************************************************************* */
