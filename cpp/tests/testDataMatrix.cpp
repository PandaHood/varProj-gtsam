//
// Created by nikolas on 7/25/25.
//


#include "../CertifiableRA.h"
#include "../CertifiablePGO.h"
#include "../CertifiableRangeSLAM.h"
#include "../CertifiableLandmark.h"
#include <gflags/gflags.h>
#include <glog/logging.h>
using namespace std;
using namespace gtsam;


#include <gtest/gtest.h>


TEST(DataMatrixTest, RotationAveragingWorks2D) {

    const int         d          = 2;
    const int         p          = 2;
    const std::string inputFile  = "/home/nikolas/StiefelManifold/data/pgo/city10000.g2o";

    size_t num_poses = 0;
    const auto measurements = DataParser::read_g2o_file(inputFile, num_poses);

    std::optional<CertificateResults> result;
    auto RA_problem = std::make_shared<CertifiableRA2>(p, measurements);

    RA_problem->setCurrentValues(RA_problem->randomInitAtLevelP(p));
    RA_problem->setCurrentGraph(RA_problem->buildGraphAtLevel(p));
    SparseMatrix M = RA_problem->recoverDataMatrixFromGraph() - RA_problem->recoverDataMatrixFromHessian();
    EXPECT_EQ(M.norm(), 0);
}

TEST(DataMatrixTest, RotationAveragingWorks3D) {

    const int         d          = 3;
    const int         p          = 3;
    const std::string inputFile  = "/home/nikolas/StiefelManifold/data/pgo/sphere2500.g2o";

    size_t num_poses = 0;
    const auto measurements = DataParser::read_g2o_file(inputFile, num_poses);

    std::optional<CertificateResults> result;
    auto RA_problem = std::make_shared<CertifiableRA3>(p, measurements);

    RA_problem->setCurrentValues(RA_problem->randomInitAtLevelP(p));
    RA_problem->setCurrentGraph(RA_problem->buildGraphAtLevel(p));
    SparseMatrix M = RA_problem->recoverDataMatrixFromGraph() - RA_problem->recoverDataMatrixFromHessian();
    EXPECT_EQ(M.norm(), 0);
}

TEST(DataMatrixTest, PGOWorks2D) {

    const int         d          = 2;
    const int         p          = 2;
    const std::string inputFile  = "/home/nikolas/StiefelManifold/data/pgo/city10000.g2o";

    size_t num_poses = 0;
    const auto measurements = DataParser::read_g2o_file(inputFile, num_poses);

    std::optional<CertificateResults> result;
    auto PGO_problem = std::make_shared<CertifiablePGO2>(p, measurements);

    PGO_problem->setCurrentValues(PGO_problem->randomInitAtLevelP(p));
    PGO_problem->setCurrentGraph(PGO_problem->buildGraphAtLevel(p));
    SparseMatrix M = PGO_problem->recoverDataMatrixFromGraph() - PGO_problem->recoverDataMatrixFromHessian();
    EXPECT_NEAR(M.norm(), 0, 1e-4);
}

TEST(DataMatrixTest, PGOWorks3D) {

    const int         d          = 3;
    const int         p          = 3;
    const std::string inputFile  = "/home/nikolas/StiefelManifold/data/pgo/sphere2500.g2o";

    size_t num_poses = 0;
    const auto measurements = DataParser::read_g2o_file(inputFile, num_poses);

    std::optional<CertificateResults> result;
    auto PGO_problem = std::make_shared<CertifiablePGO3>(p, measurements);

    PGO_problem->setCurrentValues(PGO_problem->randomInitAtLevelP(p));
    PGO_problem->setCurrentGraph(PGO_problem->buildGraphAtLevel(p));
    SparseMatrix M = PGO_problem->recoverDataMatrixFromGraph() - PGO_problem->recoverDataMatrixFromHessian();
    EXPECT_NEAR(M.norm(), 0,1e-4);
}

TEST(DataMatrixTest, RangeSLAMWorks2D) {

    const int         d          = 2;
    const int         p          = 2;
    const std::string inputFile  = "/home/nikolas/StiefelManifold/data/range/plaza1.pyfg";

    const auto measurements = DataParser::read_pycfg_file(inputFile);

    std::optional<CertificateResults> result;
    auto RA_problem = std::make_shared<CertifiableRangeSLAM2>(p, measurements);

    RA_problem->setCurrentValues(RA_problem->randomInitAtLevelP(p));
    RA_problem->setCurrentGraph(RA_problem->buildGraphAtLevel(p));
    SparseMatrix M = RA_problem->recoverDataMatrixFromGraph() - RA_problem->recoverDataMatrixFromHessian();
    EXPECT_NEAR(M.norm(), 0, 1e-4);
}

TEST(DataMatrixTest, RangeSLAMWorks3D) {

    const int         d          = 3;
    const int         p          = 3;
    const std::string inputFile  = "/home/nikolas/StiefelManifold/data/range/single_drone.pyfg";

    size_t num_poses = 0;
    const auto measurements = DataParser::read_g2o_file(inputFile, num_poses);

    std::optional<CertificateResults> result;
    auto RA_problem = std::make_shared<CertifiableRangeSLAM2>(p, measurements);

    RA_problem->setCurrentValues(RA_problem->randomInitAtLevelP(p));
    RA_problem->setCurrentGraph(RA_problem->buildGraphAtLevel(p));
    SparseMatrix M = RA_problem->recoverDataMatrixFromGraph() - RA_problem->recoverDataMatrixFromHessian();
    EXPECT_NEAR(M.norm(), 0,1e-4);
}

TEST(DataMatrixTest, LandmarkWorks2D) {

    const int         d          = 2;
    const int         p          = 2;
    const std::string inputFile  = "/home/nikolas/StiefelManifold/data/range/plaza1.pyfg";
    size_t num_poses = 0;
    const auto measurements = DataParser::read_g2o_file(inputFile, num_poses);

    std::optional<CertificateResults> result;
    auto RA_problem = std::make_shared<CertifiableLandmark2>(p, measurements);

    RA_problem->setCurrentValues(RA_problem->randomInitAtLevelP(p));
    RA_problem->setCurrentGraph(RA_problem->buildGraphAtLevel(p));
    SparseMatrix M = RA_problem->recoverDataMatrixFromGraph() - RA_problem->recoverDataMatrixFromHessian();
    EXPECT_NEAR(M.norm(), 0, 1e-4);
}


// Entry point
int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
