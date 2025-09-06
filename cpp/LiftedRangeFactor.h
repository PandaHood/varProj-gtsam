//
// Created by jason on 2/2/25.
//

#ifndef STIEFELMANIFOLDEXAMPLE_LIFTEDRANGEFACTOR_H
#define STIEFELMANIFOLDEXAMPLE_LIFTEDRANGEFACTOR_H

/**
 * @file   LiftedRangeFactor.h
 * @date
 * @author  Jason Xu, Nikolas Sanderson
 * @brief  LiftedRangeFactor
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
    class GTSAM_EXPORT LiftedRangeFactorPosePose : public NoiseModelFactorN<LiftedPoseDP, LiftedPoseDP, UnitSphereD>
     {
         double measures_; ///< dimensionality constants
         size_t p_, d_;               ///< dimensionality constants
         //        std::shared_ptr<Matrix> G_; ///< matrix of vectorized tangent space basis generators

         // Select Rot2 or Rot3 interface based template parameter d
         using Rot = typename std::conditional<d == 2, Rot2, Rot3>::type;
         using Trans = typename std::conditional<d == 2, Vector2 , Vector3>::type;

     public:

         // Provide access to the Matrix& version of evaluateError:
         using NoiseModelFactorN<LiftedPoseDP, LiftedPoseDP, UnitSphereD>::evaluateError;

         /// @name Constructor
         /// @{

         /// Constructor. Note we convert to d*p-dimensional noise model.
         LiftedRangeFactorPosePose(Key j1, Key j2, Key r1 , const double  &RangeMeasure_ , size_t p,
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
         //      ///
         Vector evaluateError(const LiftedPoseDP& P1, const LiftedPoseDP& L1, const UnitSphereD &R1, OptionalMatrixType H1,
                              OptionalMatrixType H2, OptionalMatrixType H3) const override;
         /// @}

         //    private:
         /// Calculate Jacobians if asked, Only implemented for d=2 and 3 in .cpp
         void fillJacobians(const LiftedPoseDP& P1, const LiftedPoseDP& L1, const UnitSphereD &R1, OptionalMatrixType H1,
                            OptionalMatrixType H2, OptionalMatrixType H3) const;
     };

    template <size_t d>
    class GTSAM_EXPORT LiftedRangeFactorPointPoint : public NoiseModelFactorN<Vector, Vector, UnitSphereD>
     {
         double measures_; ///< dimensionality constants
         size_t p_, d_;               ///< dimensionality constants
         //        std::shared_ptr<Matrix> G_; ///< matrix of vectorized tangent space basis generators

         // Select Rot2 or Rot3 interface based template parameter d
         using Rot = typename std::conditional<d == 2, Rot2, Rot3>::type;
         using Trans = typename std::conditional<d == 2, Vector2 , Vector3>::type;

     public:

         // Provide access to the Matrix& version of evaluateError:
         using NoiseModelFactorN<Vector, Vector, UnitSphereD>::evaluateError;

         /// @name Constructor
         /// @{

         /// Constructor. Note we convert to d*p-dimensional noise model.
         LiftedRangeFactorPointPoint(Key j1, Key j2, Key r1 , const double  &RangeMeasure_ , size_t p,
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
         //      ///
         Vector evaluateError(const Vector& P1, const Vector& L1, const UnitSphereD &R1, OptionalMatrixType H1,
                              OptionalMatrixType H2, OptionalMatrixType H3) const override;
         /// @}

         //    private:
         /// Calculate Jacobians if asked, Only implemented for d=2 and 3 in .cpp
         void fillJacobians(const Vector& P1, const Vector& L1, const UnitSphereD &R1, OptionalMatrixType H1,
                            OptionalMatrixType H2, OptionalMatrixType H3) const;
     };

    template <size_t d>
    class GTSAM_EXPORT LiftedRangeFactor : public NoiseModelFactorN<LiftedPoseDP, Vector, UnitSphereD> {
        double measures_; ///< dimensionality constants
        size_t p_, d_;               ///< dimensionality constants
//        std::shared_ptr<Matrix> G_; ///< matrix of vectorized tangent space basis generators

        // Select Rot2 or Rot3 interface based template parameter d
        using Rot = typename std::conditional<d == 2, Rot2, Rot3>::type;
        using Trans = typename std::conditional<d == 2, Vector2 , Vector3>::type;

    public:

        // Provide access to the Matrix& version of evaluateError:
        using NoiseModelFactorN<LiftedPoseDP, Vector, UnitSphereD>::evaluateError;

        /// @name Constructor
        /// @{

        /// Constructor. Note we convert to d*p-dimensional noise model.
        LiftedRangeFactor(Key j1, Key j2, Key r1 , const double  &RangeMeasure_ , size_t p,
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
//      ///
        Vector evaluateError(const LiftedPoseDP& P1, const Vector& L1, const UnitSphereD &R1, OptionalMatrixType H1,
                             OptionalMatrixType H2, OptionalMatrixType H3) const override;
        /// @}

//    private:
        /// Calculate Jacobians if asked, Only implemented for d=2 and 3 in .cpp
        void fillJacobians(const LiftedPoseDP& P1, const Vector& L1, const UnitSphereD &R1, OptionalMatrixType H1,
                           OptionalMatrixType H2, OptionalMatrixType H3) const;

        void computeHessian(std::map<std::tuple<Key,Key,Key>,Matrix> &HMap) const;


        // For building data matrix block
/************************************OLD**************************************/
//        std::vector<Eigen::Triplet<Scalar>> getRangeDataBlock(size_t num_poses, size_t num_landmarks) {
//            std::vector<Eigen::Triplet<Scalar>> triplets;
//            size_t measurement_stride = 1;
//            triplets.reserve(measurement_stride);
//            size_t i = gtsam::symbolIndex(this->key1());
//            size_t j = gtsam::symbolIndex(this->key2());
//            size_t k = gtsam::symbolIndex(this->key3());
//
//            auto diag = std::dynamic_pointer_cast<noiseModel::Diagonal>(this->noiseModel());
//            if (!diag) {
//                throw std::runtime_error("Expected a Diagonal noise model");
//            }
//            const Vector &sigmas = diag->sigmas();
//            double precision = 1 /( sigmas(0) *  sigmas(0));
//
//                triplets.emplace_back(num_poses * (d+1) + num_landmarks + k, num_poses * (d+1) + num_landmarks + k,
//                                      precision * measures_ * measures_);
//
//            return triplets;
//        }
//
//        std::vector<Eigen::Triplet<Scalar>> getTransToRangeBlock(size_t num_poses, size_t num_landmarks) {
//            std::vector<Eigen::Triplet<Scalar>> triplets;
//            size_t measurement_stride = 2;
//            triplets.reserve(measurement_stride);
//            size_t i = gtsam::symbolIndex(this->key1());
//            size_t j = gtsam::symbolIndex(this->key2());
//            size_t k = gtsam::symbolIndex(this->key3());
//
//            auto diag = std::dynamic_pointer_cast<noiseModel::Diagonal>(this->noiseModel());
//            if (!diag) {
//                throw std::runtime_error("Expected a Diagonal noise model");
//            }
//            const Vector &sigmas = diag->sigmas();
//            double precision = 1 /( sigmas(0) *  sigmas(0));
//            double w_k     = precision * measures_;
//            triplets.emplace_back(num_poses * (d+1) + num_landmarks + k, i, w_k);
//            triplets.emplace_back(num_poses * (d+1) + num_landmarks + k, j + num_poses * (d+1),  -w_k);
//
//            return triplets;
//        }
//
//        std::vector<Eigen::Triplet<Scalar>> getTransToRangeBlockTranspose(size_t num_poses, size_t num_landmarks) {
//            std::vector<Eigen::Triplet<Scalar>> triplets;
//            triplets.reserve(2);                       // ← we will insert exactly two
//            const std::size_t i = gtsam::symbolIndex(this->key1()); // first pose
//            const std::size_t j = gtsam::symbolIndex(this->key2()); // second pose
//            const std::size_t k = gtsam::symbolIndex(this->key3()); // measurement id
//            auto diag = std::dynamic_pointer_cast<noiseModel::Diagonal>(this->noiseModel());
//            if (!diag) throw std::runtime_error("Expected a Diagonal noise model");
//            const double sigma = diag->sigmas()(0);
//            const double precision = 1.0 / (sigma * sigma);         // ρ_k = 1/σ²
//            const double w_k = precision * measures_;               // ρ_k · r_k
//            const std::size_t measRow =
//                    num_poses * (d + 1)  + k + num_landmarks;
//            triplets.emplace_back(i,                       measRow, w_k);
//            triplets.emplace_back(j + num_poses * (d + 1), measRow,  -w_k);
//
//            return triplets;
//        }
//
//        std::vector<Eigen::Triplet<Scalar>> getTranslationLaplacianBlock(size_t num_poses, size_t num_landmarks) {
//            std::vector<Eigen::Triplet<Scalar>> triplets;
//            triplets.reserve(4 * d);                       // ← we will insert exactly two
//            const std::size_t i = gtsam::symbolIndex(this->key1()); // first pose
//            const std::size_t j = gtsam::symbolIndex(this->key2()); // second pose
//            const std::size_t k = gtsam::symbolIndex(this->key3()); // measurement id
//            auto diag = std::dynamic_pointer_cast<noiseModel::Diagonal>(this->noiseModel());
//            if (!diag) throw std::runtime_error("Expected a Diagonal noise model");
//            const double sigma = diag->sigmas()(0);
//            const double precision = 1.0 / (sigma * sigma);         // ρ_k = 1/σ²
//            const double w_k = precision * measures_;               // ρ_k · r_k
//            // fill the d×d identity pattern for each coordinate
//                int r1 = i;
//                int r2 = (num_poses * (d + 1)) + j;
//
//                triplets.emplace_back(r1, r1,  precision);         // L_ii
//                triplets.emplace_back(r1, r2, -precision);         // L_ij
//                triplets.emplace_back(r2, r1, -precision);         // L_ji
//                triplets.emplace_back(r2, r2,  precision);         // L_jj
//
//            return triplets;
//        }
/************************************OLD**************************************/

        // B6^T*B6
        std::vector<Eigen::Triplet<Scalar>> getRangeDataBlock(size_t num_poses, size_t num_landmarks) {
            std::vector<Eigen::Triplet<Scalar>> triplets;
            size_t measurement_stride = 1;
            triplets.reserve(measurement_stride);
            size_t k = gtsam::symbolIndex(this->key3()); // Unit sphere id
            size_t offset_range = num_poses * (d + 1) + num_landmarks;
            auto diag = std::dynamic_pointer_cast<noiseModel::Diagonal>(this->noiseModel());
            if (!diag) {
                throw std::runtime_error("Expected a Diagonal noise model");
            }
            const Vector &sigmas = diag->sigmas();
            double precision = 1 /(1 * sigmas(0) *  sigmas(0));

            triplets.emplace_back(offset_range + k, offset_range + k,
                                  precision * measures_ * measures_);

            return triplets;
        }

        // B5^T * B5
        std::vector<Eigen::Triplet<Scalar>> getLandmarkPrecisionBlock(size_t num_poses) {
            std::vector<Eigen::Triplet<Scalar>> triplets;
            triplets.reserve(1);
            const std::size_t j = gtsam::symbolIndex(this->key2()); // Landmark id
            size_t offset_pose = num_poses * (d + 1);
            auto diag = std::dynamic_pointer_cast<noiseModel::Diagonal>(this->noiseModel());
            if (!diag) throw std::runtime_error("Expected a Diagonal noise model");
            const double sigma = diag->sigmas()(0);
            const double precision = 1.0 / (1 * sigma * sigma);         // ρ_k = 1/σ²

            triplets.emplace_back(offset_pose + j, offset_pose + j,  precision);         // L_ii

            return triplets;
        }

        // B4^T * B4
        std::vector<Eigen::Triplet<Scalar>> getTranslationPrecisionBlock() {
            std::vector<Eigen::Triplet<Scalar>> triplets;
            triplets.reserve(1);
            const std::size_t i = gtsam::symbolIndex(this->key1()); // Pose id
            auto diag = std::dynamic_pointer_cast<noiseModel::Diagonal>(this->noiseModel());
            if (!diag) throw std::runtime_error("Expected a Diagonal noise model");
            const double sigma = diag->sigmas()(0);
            const double precision = 1.0 / (1 * sigma * sigma);         // ρ_k = 1/σ²
            int r1 = i;
            triplets.emplace_back(r1, r1,  precision);         // L_ii

            return triplets;
        }

        // B4^T*B6
        std::vector<Eigen::Triplet<Scalar>> getTransToRangeBlock(size_t num_poses, size_t num_landmarks) {
            std::vector<Eigen::Triplet<Scalar>> triplets;
            size_t measurement_stride = 1;
            triplets.reserve(measurement_stride);
            size_t i = gtsam::symbolIndex(this->key1()); // Pose id
            size_t k = gtsam::symbolIndex(this->key3()); // Unit sphere id

            size_t offset_range = num_poses * (d + 1) + num_landmarks;
            auto diag = std::dynamic_pointer_cast<noiseModel::Diagonal>(this->noiseModel());
            if (!diag) {
                throw std::runtime_error("Expected a Diagonal noise model");
            }
            const Vector &sigmas = diag->sigmas();
            double precision = 1 /(1 * sigmas(0) *  sigmas(0));
//            double w_k     = precision * measures_;
            double w_k     = - precision * measures_;

            triplets.emplace_back(i, offset_range + k,w_k);

            return triplets;
        }

        // B6^T*B4
        std::vector<Eigen::Triplet<Scalar>> getTransToRangeBlockTranspose(size_t num_poses, size_t num_landmarks) {
            std::vector<Eigen::Triplet<Scalar>> triplets;
            size_t measurement_stride = 1;
            triplets.reserve(measurement_stride);
            size_t i = gtsam::symbolIndex(this->key1()); // Pose id
            size_t k = gtsam::symbolIndex(this->key3()); // Unit sphere id

            size_t offset_range = num_poses * (d + 1) + num_landmarks;
            auto diag = std::dynamic_pointer_cast<noiseModel::Diagonal>(this->noiseModel());
            if (!diag) {
                throw std::runtime_error("Expected a Diagonal noise model");
            }
            const Vector &sigmas = diag->sigmas();
            double precision = 1 /(1 * sigmas(0) *  sigmas(0));
//            double w_k     = precision * measures_;
            double w_k     = - precision * measures_;

            triplets.emplace_back(offset_range + k, i, w_k);

            return triplets;
        }

        // B4^T*B5
        std::vector<Eigen::Triplet<Scalar>> getTransToLandmarkBlock(size_t num_poses) {
            std::vector<Eigen::Triplet<Scalar>> triplets;
            size_t measurement_stride = 1;
            triplets.reserve(measurement_stride);
            size_t i = gtsam::symbolIndex(this->key1()); // Pose id
            size_t j = gtsam::symbolIndex(this->key2()); // Landmark id

            size_t offset_landmark = num_poses * (d + 1);
            auto diag = std::dynamic_pointer_cast<noiseModel::Diagonal>(this->noiseModel());
            if (!diag) {
                throw std::runtime_error("Expected a Diagonal noise model");
            }
            const Vector &sigmas = diag->sigmas();
            double precision = 1 /(1 * sigmas(0) *  sigmas(0));
            triplets.emplace_back(i, offset_landmark + j,-precision);

            return triplets;
        }

        // B5^T*B4
        std::vector<Eigen::Triplet<Scalar>> getTransToLandmarkBlockTranspose(size_t num_poses) {
            std::vector<Eigen::Triplet<Scalar>> triplets;
            size_t measurement_stride = 1;
            triplets.reserve(measurement_stride);
            size_t i = gtsam::symbolIndex(this->key1()); // Pose id
            size_t j = gtsam::symbolIndex(this->key2()); // Landmark id

            size_t offset_landmark = num_poses * (d + 1);
            auto diag = std::dynamic_pointer_cast<noiseModel::Diagonal>(this->noiseModel());
            if (!diag) {
                throw std::runtime_error("Expected a Diagonal noise model");
            }
            const Vector &sigmas = diag->sigmas();
            double precision = 1 /(1 * sigmas(0) *  sigmas(0));
            triplets.emplace_back(offset_landmark + j, i,-precision);

            return triplets;
        }

        // B5^T*B6
        std::vector<Eigen::Triplet<Scalar>> getLandmarkToRange(size_t num_poses, size_t num_landmarks) {
            std::vector<Eigen::Triplet<Scalar>> triplets;
            size_t measurement_stride = 1;
            triplets.reserve(measurement_stride);
            size_t j = gtsam::symbolIndex(this->key2()); // Landmark id
            size_t k = gtsam::symbolIndex(this->key3()); // Unit sphere id

            size_t offset_landmark = num_poses * (d + 1);
            size_t offset_range = offset_landmark + num_landmarks;

            auto diag = std::dynamic_pointer_cast<noiseModel::Diagonal>(this->noiseModel());
            if (!diag) {
                throw std::runtime_error("Expected a Diagonal noise model");
            }
            const Vector &sigmas = diag->sigmas();
            double precision = 1 /(1 * sigmas(0) *  sigmas(0));
            triplets.emplace_back(offset_landmark + j, offset_range + k,precision * measures_);

            return triplets;
        }

        // B6^T*B5
        std::vector<Eigen::Triplet<Scalar>> getLandmarkToRangeTranspose(size_t num_poses, size_t num_landmarks) {
            std::vector<Eigen::Triplet<Scalar>> triplets;
            size_t measurement_stride = 1;
            triplets.reserve(measurement_stride);
            size_t j = gtsam::symbolIndex(this->key2()); // Landmark id
            size_t k = gtsam::symbolIndex(this->key3()); // Unit sphere id

            size_t offset_landmark = num_poses * (d + 1);
            size_t offset_range = offset_landmark + num_landmarks;

            auto diag = std::dynamic_pointer_cast<noiseModel::Diagonal>(this->noiseModel());
            if (!diag) {
                throw std::runtime_error("Expected a Diagonal noise model");
            }
            const Vector &sigmas = diag->sigmas();
            double precision = 1 /(1 * sigmas(0) *  sigmas(0));
            triplets.emplace_back(offset_range + k, offset_landmark + j, precision * measures_);

            return triplets;
        }

        size_t countTriplets() const {
            return 9;
        }

        void appendBlocksFromFactor(std::size_t num_poses, size_t num_landmarks, std::vector<Eigen::Triplet<Scalar>>& triplets)
        {
            const auto& t1 = this->getRangeDataBlock(num_poses, num_landmarks);
            const auto& t2 = this->getTransToRangeBlock(num_poses, num_landmarks);
            const auto& t3 = this->getTransToRangeBlockTranspose(num_poses, num_landmarks);
            const auto& t4 = this->getTranslationPrecisionBlock();
            const auto& t5 = this->getLandmarkPrecisionBlock(num_poses);
            const auto& t6 = this->getLandmarkToRange(num_poses, num_landmarks);
            const auto& t7 = this->getLandmarkToRangeTranspose(num_poses, num_landmarks);
            const auto& t8 = this->getTransToLandmarkBlock(num_poses);
            const auto& t9 = this->getTransToLandmarkBlockTranspose(num_poses);


            for (const auto* block : { &t1, &t2, &t3, &t4, &t5, &t6, &t7, &t8, &t9})
                triplets.insert(triplets.end(), block->begin(), block->end());
        }

    };

// Explicit instantiation for d=2 and d=3 in .cpp file:
    using LiftedRangeFactor2 = LiftedRangeFactor<2>;
    using LiftedRangeFactor3 = LiftedRangeFactor<3>;

} // namespace gtsam





#endif //STIEFELMANIFOLDEXAMPLE_LIFTEDRANGEFACTOR_H
