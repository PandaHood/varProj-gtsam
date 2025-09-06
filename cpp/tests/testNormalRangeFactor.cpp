/* ----------------------------------------------------------------------------
 * Copyright 2025, Northeastern University Robust Autonomy Lab
 * All Rights Reserved
 * Authors: Zhexin Xu
 * See LICENSE for the license information
 * -------------------------------------------------------------------------- */

#include "../StiefelManifold.h"
#include "../StiefelManifold-inl.h"
#include "../SEsyncFactor.h"
#include "../NormalRangeFactor.h"

#include <gtsam/nonlinear/factorTesting.h>
#include <gtsam/base/Matrix.h>
#include <gtsam/geometry/Rot2.h>
#include <gtsam/geometry/Rot3.h>
#include <CppUnitLite/TestHarness.h>

#include <random>

using namespace gtsam;

/* =============================== 2D TESTS =============================== */

// EvaluateError returns 1D residual; sizes and basic call should succeed.
TEST(NormalRangeFactor, EvaluateError_2D_Size) {
    const size_t d = 2;
    const int p = 3;

    // Simple rotation block in top-left
    Rot2 R(0.3);
    Matrix M = Matrix::Zero(p, d);
    M.block(0, 0, 2, 2) = R.matrix();

    StiefelManifoldKP Q1(M);
    Q1.initializeTangentSpace();

    // Two translations in R^p (lifted) created from R^d vectors
    Vector t1(2), t2(2);
    t1 << 0.05, 0.05;
    t2 << 0.10, 0.20;

    LiftedPoseDP P1(Q1, LiftedPoseDP::LiftToRp(t1, p));
    Vector       L1 = LiftedPoseDP::LiftToRp(t2, p);

    const double measured = (t1 - t2).norm();
    auto model = noiseModel::Isotropic::Sigma(1, 0.1);

    NormalRangeFactor<2> factor(1, 2, measured, p, model);

    Matrix H1, H2;
    Vector err = factor.evaluateError(P1, L1, H1, H2);

    EXPECT_LONGS_EQUAL(1, err.size());
    EXPECT_LONGS_EQUAL(1, H1.rows());
    EXPECT_LONGS_EQUAL(1, H2.rows());
}

// Numerical Jacobians in 2D across several p
TEST(NormalRangeFactor, NumericalJacobians_2D) {
    const size_t d = 2;
    const Rot2 R(0.3);

    Vector t1(2), t2(2);
    t1 << 1.0, 1.0;
    t2 << 6.0, 7.0;

    // any scalar is fine; Jacobians of ||tP - tL|| - c do not depend on c
    const double measured = (t2 - t1).norm();

    for (const int p : {5, 4, 3, 2}) {
        Matrix M = Matrix::Zero(p, d);
        M.block(0, 0, 2, 2) = R.matrix();

        StiefelManifoldKP Q1(M);
        Q1.initializeTangentSpace();

        LiftedPoseDP P1(Q1, LiftedPoseDP::LiftToRp(t1, p));
        Vector       L1 = LiftedPoseDP::LiftToRp(t2, p);

        auto model = noiseModel::Isotropic::Sigma(1, 1.2);
        NormalRangeFactor<2> factor(11, 22, measured, p, model);

        Values values;
        values.insert(11, P1);
        values.insert(22, L1);

        EXPECT_CORRECT_FACTOR_JACOBIANS(factor, values, 1e-7, 1e-5);
    }
}

// Zero-distance case in 2D should yield zero Jacobians (robust branch)
TEST(NormalRangeFactor, ZeroDistanceJacobian_2D) {
    const size_t d = 2; const int p = 4;

    Matrix M = Matrix::Zero(p, d);
    M.block(0, 0, 2, 2) = Rot2(0.0).matrix();

    StiefelManifoldKP Q1(M);
    Q1.initializeTangentSpace();

    Vector t(2); t << 0.0, 0.0;
    LiftedPoseDP P1(Q1, LiftedPoseDP::LiftToRp(t, p));
    Vector       L1 = LiftedPoseDP::LiftToRp(t, p);  // same location → dist = 0

    const double measured = 0.0;
    auto model = noiseModel::Isotropic::Sigma(1, 0.1);
    NormalRangeFactor<2> factor(7, 8, measured, p, model);

    Matrix H1, H2;
    Vector1 err = factor.evaluateError(P1, L1, H1, H2);

    // With identical points, gradient should be zero per your implementation
    EXPECT(assert_equal(Matrix::Zero(1, LiftedPoseDP::Dimension(d, p)), H1, 1e-12));
    EXPECT(assert_equal(Matrix::Zero(1, p), H2, 1e-12));
}

/* =============================== 3D TESTS =============================== */

// EvaluateError returns 1D residual; sizes and call should succeed in 3D.
TEST(NormalRangeFactor, EvaluateError_3D_Size) {
    const size_t d = 3;
    const int p = 5;

    // Put a 3x3 rotation in the top-left of a p x d Stiefel element
    Rot3 R = Rot3::RzRyRx(0.1, -0.05, 0.2);
    Matrix M = Matrix::Zero(p, d);
    M.topLeftCorner(3, 3) = R.matrix();

    StiefelManifoldKP Q1(M);
    Q1.initializeTangentSpace();

    Vector t1(3), t2(3);
    t1 << 0.0, 0.0, 0.1;
    t2 << 0.6, 0.7, 0.7;

    LiftedPoseDP P1(Q1, LiftedPoseDP::LiftToRp(t1, p));
    Vector       L1 = LiftedPoseDP::LiftToRp(t2, p);

    const double measured = (t1 - t2).norm();
    auto model = noiseModel::Isotropic::Sigma(1, 0.1);

    NormalRangeFactor<3> factor(1, 2, measured, p, model);

    Matrix H1, H2;
    Vector err = factor.evaluateError(P1, L1, H1, H2);

    EXPECT_LONGS_EQUAL(1, err.size());
    EXPECT_LONGS_EQUAL(1, H1.rows());
    EXPECT_LONGS_EQUAL(1, H2.rows());
}

// Numerical Jacobians in 3D across several p
TEST(NormalRangeFactor, NumericalJacobians_3D) {
    const size_t d = 3;

    Vector t1(3), t2(3);
    t1 << 0.0, 0.0, 0.1;
    t2 << 0.6, 0.7, 0.7;

    const double measured = (t2 - t1).norm();

    for (const int p : {6, 5, 4, 3}) {
        Rot3 R = Rot3::RzRyRx(0.1, 0.05, -0.15);
        Matrix M = Matrix::Zero(p, d);
        M.topLeftCorner(3, 3) = R.matrix();

        StiefelManifoldKP Q1(M);
        Q1.initializeTangentSpace();

        LiftedPoseDP P1(Q1, LiftedPoseDP::LiftToRp(t1, p));
        Vector       L1 = LiftedPoseDP::LiftToRp(t2, p);

        auto model = noiseModel::Isotropic::Sigma(1, 1.2);
        NormalRangeFactor<3> factor(101, 202, measured, p, model);

        Values values;
        values.insert(101, P1);
        values.insert(202, L1);

        EXPECT_CORRECT_FACTOR_JACOBIANS(factor, values, 1e-7, 1e-5);
    }
}

// Zero-distance case in 3D should yield zero Jacobians
TEST(NormalRangeFactor, ZeroDistanceJacobian_3D) {
    const size_t d = 3; const int p = 5;

    Matrix M = Matrix::Zero(p, d);
    M.topLeftCorner(3, 3) = Rot3().matrix();

    StiefelManifoldKP Q1(M);
    Q1.initializeTangentSpace();

    Vector t(3); t << 0.0, 0.0, 0.0;
    LiftedPoseDP P1(Q1, LiftedPoseDP::LiftToRp(t, p));
    Vector       L1 = LiftedPoseDP::LiftToRp(t, p);

    const double measured = 0.0;
    auto model = noiseModel::Isotropic::Sigma(1, 0.1);
    NormalRangeFactor<3> factor(5, 6, measured, p, model);

    Matrix H1, H2;
    Vector1 err = factor.evaluateError(P1, L1, H1, H2);

    EXPECT(assert_equal(Matrix::Zero(1, LiftedPoseDP::Dimension(d, p)), H1, 1e-12));
    EXPECT(assert_equal(Matrix::Zero(1, p), H2, 1e-12));
}

/* ************************************************************************* */
int main() {
    TestResult tr;
    return TestRegistry::runAllTests(tr);
}
/* ************************************************************************* */
