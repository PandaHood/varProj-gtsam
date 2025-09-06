//
// Created by nikolas on 5/2/25.
//

#ifndef CERTIFIABLELANDMARK_H
#define CERTIFIABLELANDMARK_H

#pragma once
#include "Certifiable_problem.h"
#include "LandmarkFactor.h"

/**
 * @brief Certifiable landmark estimation problem.
 *
 * Implements certifiable estimation for simultaneous pose‐and‐landmark SLAM
 * problems via relaxation levels p.
 *
 * @tparam d  Ambient dimension (must be 2 or 3).
 */
template <size_t d>
class CertifiableLandmark : public CertifiableProblem
{
    static_assert(d == 2 || d == 3, "CertifiableLandmark only supports d = 2 or 3.");

    /// Number of landmarks in the problem.
    size_t num_landmark_;

    /// Number of landmark measurements.
    size_t num_lmk_measurements_;

    /// Number of pose‐to‐pose measurements.
    size_t num_pose_measurements_;

public:

    /**
     * @brief Constructor.
     * @param p             Initial relaxation rank.
     * @param measurements  Parsed measurements struct (num_poses, landmarks, etc.).
     */
    CertifiableLandmark(size_t p, const DataParser::Measurement& measurements, const std::string g2oPath)
        : CertifiableProblem(d, p, measurements)
    {
        num_landmark_ = measurements.num_landmarks;
        dataPath = g2oPath;
    }

    /**
     * @brief Initialize the problem: build graph, random init, and recover data matrix.
     */
    void init() {
        num_lmk_measurements_  = measurements_.landmarkMeasurements.size();
        num_pose_measurements_ = measurements_.poseMeasurements.size();
        auto t0 = CFGStopwatch::tick();
        if (opts_.initType == CertifiableProblemOpts::InitType::Random){
            currentValues_  = randomInitAtLevelP(currentRank_);
        }
        else if (opts_.initType == CertifiableProblemOpts::InitType::Odom){
            currentValues_  = InitializationFromOdom(dataPath, currentRank_);
        }
        else {
            throw std::runtime_error("Unknown initialization type in CertifiableProblemOpts::initType");
        }
        currentGraph_   = buildGraphAtLevel(currentRank_);
        M_              = recoverDataMatrixFromHessian();
        auto t1 = CFGStopwatch::tock(t0);
        certificateResults_.initialization_time.push_back(t1);
    }

    /**
     * @brief Get the number of landmarks.
     * @return Number of landmarks.
     */
    inline size_t getNumLandmark() const { return num_landmark_; }

    /**
     * @brief Build the factor graph at relaxation level p.
     * @param p  Relaxation rank.
     * @return   A GTSAM NonlinearFactorGraph containing both pose and landmark factors.
     */
    NonlinearFactorGraph buildGraphAtLevel(size_t p) override
    {
        NonlinearFactorGraph inputGraph;

        // Pose factors
        for (const auto& meas : measurements_.poseMeasurements)
        {
            const double Kappa = meas.kappa;
            const double tau   = meas.tau;
            Vector sigmas = Vector::Zero(p * d + p);
            sigmas.head(p * d).setConstant(std::sqrt(1.0 / (1 * Kappa)));
            sigmas.tail(p).setConstant(std::sqrt(1.0 / (1 * tau)));
            auto noise = noiseModel::Diagonal::Sigmas(sigmas);

            if constexpr (d == 2) {
                inputGraph.emplace_shared<SEsyncFactor2>(
                  meas.i, meas.j, meas.R, meas.t, p, noise);
            } else {
                inputGraph.emplace_shared<SEsyncFactor3>(
                  meas.i, meas.j, meas.R, meas.t, p, noise);
            }
        }

        // Landmark factors
        for (const auto& meas : measurements_.landmarkMeasurements)
        {
            const double nu = meas.nu;
            Vector sigmas = Vector::Zero(p);
            sigmas.tail(p).setConstant(std::sqrt(1.0 / (1 * nu)));
            auto noise = noiseModel::Diagonal::Sigmas(sigmas);

            if constexpr (d == 2) {
                inputGraph.emplace_shared<LiftedLandmarkFactor2>(
                  meas.i, Symbol('L', meas.j), meas.l, p, noise);
            } else {
                inputGraph.emplace_shared<LiftedLandmarkFactor3>(
                  meas.i, Symbol('L', meas.j), meas.l, p, noise);
            }
        }

        return inputGraph;
    }

    // need to create incremental version
    NonlinearFactorGraph buildGraphAtLevelUsingGraph(size_t p) override {
        throw("Not implemented");
    }

    /**
     * @brief Compute block‐diagonal lambda matrix blocks from Y.
     * @param Y  Variable matrix.
     * @return   Dense block matrix of size (d × d·num_poses).
     */
    Matrix computeLambdaBlocks(const Matrix& Y) override {
        Matrix SY   = M_ * Y;
        Matrix Yt   = Y.transpose();
        Matrix LambdaBlocks(d, num_pose_ * d);
        size_t offset = num_pose_;

        for (size_t i = 0; i < num_pose_; ++i) {
            Matrix P = SY.block(offset + i * d, 0, d, Y.cols()) *
                       Yt.block(0, offset + i * d, Y.cols(), d);
            LambdaBlocks.block(0, i * d, d, d) = 0.5 * (P + P.transpose());
        }
        return LambdaBlocks;
    }

    /**
     * @brief Convert dense lambda blocks into a sparse lambda matrix.
     * @param LambdaBlocks  Dense block matrix.
     * @return              Sparse certificate matrix Λ.
     */
    SparseMatrix computeLambdaFromLambdaBlocks(const Matrix& LambdaBlocks) override {
        std::vector<Eigen::Triplet<Scalar>> elements;
        elements.reserve(d * d * num_pose_);
        size_t offset = num_pose_;

        for (size_t i = 0; i < num_pose_; ++i)
            for (size_t r = 0; r < d; ++r)
                for (size_t c = 0; c < d; ++c)
                    elements.emplace_back(
                        offset + i * d + r,
                        offset + i * d + c,
                        LambdaBlocks(r, i * d + c));

        SparseMatrix Lambda(
            offset + d * num_pose_ + num_landmark_,
            offset + d * num_pose_ + num_landmark_);
        Lambda.setFromTriplets(elements.begin(), elements.end(),
                               std::plus<Scalar>());
        return Lambda;
    }

    /**
     * @brief Build the sparse element matrix S from Values.
     * @param values  GTSAM Values containing LiftedPoseDP and raw vectors.
     * @return        Sparse matrix S of size.
     */
    SparseMatrix elementMatrix(const Values& values) override {
        const size_t N = num_pose_;
        const size_t p = currentRank_;
        const size_t L = num_landmark_;
        const size_t offset_pose = N * (d + 1);

        std::vector<Eigen::Triplet<double>> triplets;
        triplets.reserve(N * p * (d + 1) + L * p);

        // Pose entries
        for (const auto& kv : values.extract<LiftedPoseDP>()) {
            size_t i = kv.first;
            const auto& pose = kv.second;
            const auto& t   = pose.get_TranslationVector(); // p × 1
            const auto& mat = pose.matrix();                // p × (d+1), but we only need the left-upper block(p x d)

            for (size_t row = 0; row < p; ++row)
                triplets.emplace_back(i, row, t(row));

            size_t row0 = N + i * d;
            for (size_t row = 0; row < p; ++row)
                for (size_t col = 0; col < d; ++col)
                    triplets.emplace_back(row0 + col, row, mat(row, col));
        }

        // Landmark entries
        for (const auto& kv : values.extract<Vector>()) {
            size_t i = kv.first;
            const auto& t = kv.second; // p × 1
            for (size_t row = 0; row < p; ++row)
                triplets.emplace_back(offset_pose + i, row, t(row));
        }

        SparseMatrix S(N + N * d + L, p);
        S.setFromTriplets(triplets.begin(), triplets.end(),
                          std::plus<Scalar>());
        return S;
    }

    /**
     * @brief Randomly initialize values at relaxation level Pmin.
     * @param Pmin  relaxation level.
     * @return      GTSAM Values with random LiftedPoseDP and landmark vectors.
     */
    Values randomInitAtLevelP(const size_t Pmin) override
    {
        Values initial;
        // Poses
        for (const auto& meas : measurements_.poseMeasurements) {
            if (!initial.exists(meas.i)) {
                StiefelManifoldKP Y =
                    StiefelManifoldKP::Random(
                        std::default_random_engine::default_seed, d, Pmin);
                Vector trans = Vector::Random(Pmin);
                initial.insert(meas.i, LiftedPoseDP(Y, trans));
            }
            if (!initial.exists(meas.j)) {
                StiefelManifoldKP Y =
                    StiefelManifoldKP::Random(
                        std::default_random_engine::default_seed, d, Pmin);
                Vector trans = Vector::Random(Pmin);
                initial.insert(meas.j, LiftedPoseDP(Y, trans));
            }
        }
        // Landmarks
        for (const auto& meas : measurements_.landmarkMeasurements) {
            Key key = Symbol('L', meas.j);
            if (!initial.exists(key)) {
                Vector trans = Vector::Random(Pmin);
                initial.insert(key, trans);
            }
        }
//        initial.print("Test initials construction:");
        return initial;
    }

    /**
     * @brief Convert a flat eigenvector into per-variable tangent VectorValues.
     * @param p      Relaxation level.
     * @param v      Flattened descent vector size.
     * @param values Lifted Values at level p.
     * @return       VectorValues for each variable on its tangent space.
     */
    VectorValues TangentVectorValues(size_t p,
                                     const Vector v,
                                     const Values values) override {
        VectorValues delta;
        Matrix Ydot = Matrix::Zero(v.size(), p);
        Ydot.rightCols<1>() = v;

        const auto poses     = values.extract<LiftedPoseDP>();
        const auto landmarks = values.extract<Vector>();
        size_t offset_pose = num_pose_ * (d + 1);

        // Poses tangent blocks
        for (const auto& kv : poses) {
            size_t idx = gtsam::symbolIndex(kv.first);
            const auto& pose = kv.second;
            StiefelManifoldKP Y = pose.get_Y();

            Matrix tangM = Ydot.block(num_pose_ + idx * d, 0, d, p).transpose();
            Vector transV = Ydot.block(idx, 0, 1, p).transpose();

            Vector xi = StiefelManifoldKP::Vectorize(tangM);
            Vector tangentVec = Y.G_.transpose() * xi;

            Vector combined(tangentVec.size() + transV.size());
            combined << tangentVec, transV;
            delta.insert(idx, combined);
        }

        // Landmarks tangent blocks
        for (const auto& kv : landmarks) {
            size_t idx = gtsam::symbolIndex(kv.first);
            Vector lmkV = Ydot.block(offset_pose + idx, 0, 1, p).transpose();
            delta.insert(num_pose_ + idx, lmkV);
        }
        return delta;
    }

    /**
     * @brief Project ambient gradient Ydot onto tangent space at Y.
     * @param p     Relaxation level.
     * @param Y     Basepoint matrix.
     * @param Ydot  Ambient gradient matrix.
     * @return      Tangent‐space gradient at Y.
     */
    Matrix tangent_space_projection(const size_t p,
                                    const Matrix& Y,
                                    const Matrix& Ydot) override {
        Matrix result = Ydot;
        size_t offset_t = num_pose_;
        size_t rot_sz   = d * num_pose_;

        result.block(offset_t, 0, rot_sz, p) =
            StiefelManifoldKP::Proj(
                Y.block(0, 0, rot_sz, p).transpose(),
                result.block(0, 0, rot_sz, p).transpose()
            ).transpose();
        return result;
    }

    /**
     * @brief Recover the data/certificate matrix from the current factor graph.
     * @return Sparse data matrix L.
     */
    SparseMatrix recoverDataMatrixFromGraph() override {
        std::vector<Eigen::Triplet<Scalar>> triplets;
        triplets.reserve(
            num_pose_measurements_ * (4 + 2*d*d + 2*d + 2) +
            num_lmk_measurements_  * (d*d + 4*d + 4)
        );

        if (d == 2) {
            for (const auto& f_ptr : currentGraph_) {
                if (auto factor = std::dynamic_pointer_cast<SEsyncFactor2>(f_ptr)) {
                    factor->appendBlocksFromFactor(num_pose_, triplets);
                }
                if (auto factor = std::dynamic_pointer_cast<LiftedLandmarkFactor2>(f_ptr)) {
                    factor->appendBlocksFromFactor(num_pose_, triplets);
                }
            }
        } else {
            for (const auto& f_ptr : currentGraph_) {
                if (auto factor = std::dynamic_pointer_cast<SEsyncFactor3>(f_ptr)) {
                    factor->appendBlocksFromFactor(num_pose_, triplets);
                }
                if (auto factor = std::dynamic_pointer_cast<LiftedLandmarkFactor3>(f_ptr)) {
                    factor->appendBlocksFromFactor(num_pose_, triplets);
                }
            }
        }

        SparseMatrix L((d + 1) * num_pose_ + num_landmark_,
                       (d + 1) * num_pose_ + num_landmark_);
//        L.setFromTriplets(triplets.begin(), triplets.end(),
//                          std::plus<Scalar>());
        L.setFromTriplets(triplets.begin(), triplets.end());
        return L;
    }

    SparseMatrix recoverDataMatrixFromHessian()
    {

        std::map<std::pair<Key,Key>,Matrix> HMap;
        std::map<std::pair<Key,Key>,Matrix> HMap_landmark;
        recoverHessian(HMap, HMap_landmark);
        Ordering order = Ordering::Natural(currentGraph_);
        // Build the index in one shot:
        const size_t N = order.size();
        std::vector<size_t> dims(N), offsets(N);

        std::map<Key,size_t> Key2index;
        for (size_t i = 0; i < N; ++i) {
            const Key k = order[i];
            // look at the letter in the key:
            if (symbolChr(k)=='L') {
                dims[i] = 1;
            }
            else {
                // assume it's your LiftedPoseDP
                auto lp = currentValues_.at<LiftedPoseDP>(k);
                dims[i] = lp.get_Cols();
            }
            offsets[i] = (i==0 ? 0u : offsets[i-1] + dims[i-1]);
            Key2index[k] = i;
        }
        size_t totalDim = offsets.back() + dims.back();
        // 2) Gather triplets
        std::vector<Eigen::Triplet<Scalar>> triplets;
        triplets.reserve(num_pose_measurements_ * (4 + 2*d*d + 2*d + 2) +
            num_lmk_measurements_  * (d*d + 4*d + 4));

        appendPGOHessian(triplets,HMap, offsets, order, Key2index);
        appendLandmarkHessian(triplets,HMap_landmark, offsets, order, Key2index);

        SparseMatrix L(totalDim, totalDim);
        L.setFromTriplets(triplets.begin(), triplets.end());
        // 1) make sure the sparse structure is finalized
        L.makeCompressed();

        return L;
    }

    void recoverHessian(std::map<std::pair<Key,Key>,Matrix> &HMap, std::map<std::pair<Key,Key>,Matrix> &HMap_landmark)
    {
        if (d == 2) {
            for (auto& f_ptr : currentGraph_) {
                if (auto factor = std::dynamic_pointer_cast<SEsyncFactor2>(f_ptr))
                {
                    factor->computeHessian(HMap);
                }
                if (auto factor = std::dynamic_pointer_cast<LiftedLandmarkFactor2>(f_ptr))
                {
                    factor->computeHessian(HMap_landmark);
                }
            }
        } else
        {
            for (auto& f_ptr : currentGraph_)
            {
                if (auto factor = std::dynamic_pointer_cast<SEsyncFactor3>(f_ptr))
                {
                    factor->computeHessian(HMap);
                }
                if (auto factor = std::dynamic_pointer_cast<LiftedLandmarkFactor3>(f_ptr))
                {
                    factor->computeHessian(HMap_landmark);
                }
            }
        }
    }

    void appendPGOHessian(std::vector<Eigen::Triplet<Scalar>> &triplets, std::map<std::pair<Key,Key>,Matrix> &HMap, std::vector<size_t> offsets, Ordering order, std::map<Key,size_t> Key2index)
    {
        for (auto const& kv : HMap) {
            Key const& k1 = kv.first.first;
            Key const& k2 = kv.first.second;
            size_t i = symbolIndex(k1), j = symbolIndex(k2);

            Matrix Hfull = kv.second;  // 4d×4d block in [Yi,ti,Yj,tj] order


            for (int bi = 0; bi < 4; ++bi) {
                const bool rowRot = (bi == 0 || bi == 2);
                const int  rowDim = rowRot ? d : 1;
                const int  rowBase = rowRot
                  ? int(num_pose_ + d * (bi < 2 ? i : j))
                  :     int(    (bi < 2 ? i : j));


                for (int bj = 0; bj < 4; ++bj) {
                    const bool colRot = (bj == 0 || bj == 2);
                    const int  colDim = colRot ? d : 1;
                    const int  colBase = colRot
                      ? int(num_pose_ + d * (bj < 2 ? i : j))
                      :     int(    (bj < 2 ? i : j));

                    for (int r = 0; r < rowDim; ++r) {
                        for (int c = 0; c < colDim; ++c) {
                            // always index Hfull in block‐row‐block‐col order:
                            int Hrow = bi*d + r;
                            int Hcol = bj*d + c;
                            const Scalar v = Hfull(Hrow, Hcol);
                            if (v != Scalar(0)) {
                                triplets.emplace_back(rowBase + r,
                                                      colBase + c,
                                                      v);
                            }
                        }
                    }
                }
            }
        }
    }

    void appendLandmarkHessian(std::vector<Eigen::Triplet<Scalar>> &triplets, std::map<std::pair<Key,Key>,Matrix> &HMap_landmark, std::vector<size_t> offsets, Ordering order, std::map<Key,size_t> Key2index)
    {
        for (auto const& kv : HMap_landmark) {
            const Key&   k1    = kv.first.first;
            const Key&   k2    = kv.first.second;
            size_t       i     = symbolIndex(k1);
            size_t       j     = Key2index[k2];
            const size_t ti_base = i;                   // tᵢ starts at row/col i
            const size_t li_base = offsets[j];          // ℓᵢ starts at row/col offsets[j]
            const size_t ri_base = num_pose_ + d * i;   // rᵢ (rotation) starts here
            const Matrix  Hfull = kv.second;            // (3d × 3d) Hessian block

            // — tᵢ–tᵢ, tᵢ–ℓᵢ, ℓᵢ–tᵢ, ℓᵢ–ℓᵢ —
            triplets.emplace_back(ti_base,    ti_base,    Hfull(    0,     0));
            triplets.emplace_back(ti_base,    li_base,    Hfull(    0,     d));
            triplets.emplace_back(li_base,    ti_base,    Hfull(    d,     0));
            triplets.emplace_back(li_base,    li_base,    Hfull(    d,     d));

            // — tᵢ–rᵢ —
            for (size_t r = 0; r < d; ++r) {
                triplets.emplace_back(
                    ti_base,
                    ri_base + r,
                    Hfull(/* row = */ r,       /* col = */ 2*d)
                );
            }

            // — rᵢ–tᵢ —
            for (size_t r = 0; r < d; ++r) {
                triplets.emplace_back(
                    ri_base + r,
                    ti_base,
                    Hfull(/* row = */ 2*d,    /* col = */ r)
                );
            }

            // — ℓᵢ–rᵢ —
            for (size_t r = 0; r < d; ++r) {
                triplets.emplace_back(
                    li_base,
                    ri_base + r,
                    Hfull(/* row = */ d + r,  /* col = */ 2*d)
                );
            }

            // — rᵢ–ℓᵢ —
            for (size_t r = 0; r < d; ++r) {
                triplets.emplace_back(
                    ri_base + r,
                    li_base,
                    Hfull(/* row = */ 2*d,    /* col = */ d + r)
                );
            }

            // — rᵢ–rᵢ —
            for (size_t rr = 0; rr < d; ++rr) {
                for (size_t cc = 0; cc < d; ++cc) {
                    triplets.emplace_back(
                        ri_base + rr,
                        ri_base + cc,
                        Hfull(/* row = */ 2*d + rr, /* col = */ 2*d + cc)
                    );
                }
            }
        }
    }



    /**
     * @brief Clamp a value into [lower_bound, upper_bound].
     * @param val          Input value.
     * @param lower_bound  Minimum allowed.
     * @param upper_bound  Maximum allowed.
     * @return             Clamped value.
     */
    Scalar thresholdVal(Scalar val,
                        Scalar lower_bound,
                        Scalar upper_bound) {
        if (val < lower_bound) return lower_bound;
        if (val > upper_bound) return upper_bound;
        return val;
    }

    Matrix dataMatrixProduct(const Matrix &Y, const SparseMatrix &M) {
        return M * Y;
    }

    Scalar evaluateObjective(const Matrix &Y, const SparseMatrix &M) {
        return  0.5 * (Y.transpose() * dataMatrixProduct(Y, M)).trace();
    }

    /**
     * @brief Attempt to find a certifiable solution between pMin and pMax.
     *
     * Runs optimization, gradient checks, and fast_verification at each level;
     * performs level‐lifting/saddle-escape upon failure until success or exhaustion.
     *
     * @param pMin  Minimum relaxation rank.
     * @param pMax  Maximum relaxation rank.
     * @return      The smallest p for which the solution is certified, or std::nullopt.
     */
    std::optional<CertificateResults> Solve(size_t pMin, size_t pMax) {
        Values Qstar;
        auto t6 = std::chrono::high_resolution_clock::now();
        for (size_t p = pMin; p <= pMax; ++p) {
            std::cout << "Starting optimization at rank = " << p << std::endl;
            setCurrentRank(p);
            Qstar = tryOptimizingAtLevel(p);
            setCurrentValues(Qstar);

//            auto t2 = CFGStopwatch::tick();
            auto nonlinear_graph = buildGraphAtLevel(p);
            auto linear_graph = nonlinear_graph.linearize(Qstar);
//            auto grad_norm = linear_graph->gradientAtZero();
//            std::cout << "Gradient norm at level p = " << p << " is : " << grad_norm.norm() << std::endl;
//            std::cout << "Objective value is : " << std::scientific << std::setprecision(6) << nonlinear_graph.error(Qstar) << std::endl;
//            certificateResults_.gradnorm.push_back(grad_norm.norm());
//            auto t3 = CFGStopwatch::tock(t2);
//            certificateResults_.initialization_time.push_back(t3);

            auto t4 = CFGStopwatch::tick();
            SparseMatrix S = elementMatrix(Qstar);
            Matrix lambdaBlocks = computeLambdaBlocks(S);
            SparseMatrix Lambda = computeLambdaFromLambdaBlocks(lambdaBlocks);
            SparseMatrix M = getDataMatrix();

            Scalar  obj = evaluateObjective(S, M);

            // For test
            std::cout << "GTSAM's error : " << nonlinear_graph.error(Qstar) <<std::endl;
            std::cout << "SDP's error : " << obj  <<std::endl;

            bool success = false;
            Scalar eta;
            if (opts_.useAbsoluteEta == true) {
                eta = opts_.eta;
            } else {
                eta = thresholdVal(obj * opts_.REL_CERT_ETA, opts_.MIN_CERT_ETA, opts_.MAX_CERT_ETA);
            }

            size_t nx = opts_.nx;
            Vector v;
            Scalar theta;
            size_t num_lobpcg_iters;
            size_t max_iters = opts_.max_iters;
            Scalar max_fill_factor =  opts_.max_fill_factor;
            Scalar drop_tol = opts_.drop_tol;

            success = fast_verification(M - Lambda, eta, nx, theta, v,
                                        num_lobpcg_iters, max_iters, max_fill_factor, drop_tol);
            auto t5 = CFGStopwatch::tock(t4);
            certificateResults_.verification_times.push_back(t5);
            if (!success) {
                std::cout << "Using eta = " << eta << ", Theta : " << theta << std::endl;
                increaseCurrentRank();
                currentValues_ = initializeWithDescentDirection(Qstar, M, v, theta, 1e-2);
            } else {
                std::cout << "Solution verified at level p = " << p << std::endl;
                std::cout << "Final eta = " << eta << ", Theta : " << theta << std::endl;
                certificateResults_.Yopt = S;
                certificateResults_.Lambda = Lambda;
                certificateResults_.xhat = RoundSolutionS();
                certificateResults_.total_computation_time = CFGStopwatch::tock(t6);
                certificateResults_.endingRank = p;
                return certificateResults_;
            }
        }

        std::cout << "No certifiable solution found in p ∈ [" << pMin << ", " << pMax << "]" << std::endl;
        return std::nullopt;
    }

    /**
     * @brief Round the relaxed solution back to problem dimension.
     * @return Matrix R.
     */
    Matrix RoundSolutionS() override
    {
        Matrix S = elementMatrix(currentValues_).transpose();


        // First, compute a thin SVD of Y
        Eigen::JacobiSVD<Matrix> svd(S, Eigen::ComputeFullV);

        Vector sigmas = svd.singularValues();

        // Construct a diagonal matrix comprised of the first d singular values
        DiagonalMatrix Sigma_d(d);
        DiagonalMatrix::DiagonalVectorType &diagonal = Sigma_d.diagonal();
        for (size_t i = 0; i < d; ++i)
            diagonal(i) = sigmas(i);

        // First, construct a rank-d truncated singular value decomposition for Y
        Matrix R = Sigma_d * svd.matrixV().leftCols(d).transpose();
        Vector determinants(num_pose_);

        // Compute the offset at which the rotation matrix blocks begin
        size_t rot_offset;
        rot_offset = num_pose_;


        size_t ng0 = 0; // This will count the number of blocks whose
        // determinants have positive sign
        for (size_t i = 0; i < num_pose_; i++) {
            // Compute the determinant of the ith dxd block of R
            determinants(i) = R.block(0, rot_offset + i * d, d, d).determinant();
            if (determinants(i) > 0)
                ++ng0;
        }

        if (ng0 < num_pose_ / 2) {
            // Less than half of the total number of blocks have the correct sign, so
            // reverse their orientations

            // Get a reflection matrix that we can use to reverse the signs of those
            // blocks of R that have the wrong determinant
            Matrix reflector = Matrix::Identity(d, d);
            reflector(d - 1, d - 1) = -1;

            R = reflector * R;
        }

        // Finally, project each dxd rotation block to SO(d)
        for (size_t i = 0; i < num_pose_; i++) {
            R.block(0, rot_offset + i * d, d, d) = project_to_SOd(R.block(0, rot_offset + i * d, d, d));
        }
        return R;
    }

    /**
     * @brief Export the solution in G2O or TUM format.
     * @param path  Base filename (without extension).
     * @param R     Rotation/translation solution matrix.
     * @param g2o   If true, write G2O; otherwise, write TUM.
     */
    void ExportData(const string& path, const Eigen::MatrixXd& R, bool g2o) override
    {
        if (g2o) {
            std::ofstream file(path + ".g2o");
            if (d == 2) {
                for (auto i = 0; i < num_pose_; ++i) {
                    Eigen::Matrix3d R3 = Eigen::Matrix3d::Identity();
                    R3.topLeftCorner<2,2>() =
                            R.block(0, num_pose_ + i*2, 2, 2);
                    auto theta = std::atan2(R3(1,0), R3(0,0));
                    Vector t = R.block(0, i, d, 1);
                    file << "VERTEX_SE2" << " "<< i << " " << t(0) << " " << t(1)
                         << " " << theta << "\n";
                }

                for (auto i = 0; i < num_landmark_; ++i)
                {
                    Vector t = R.block(0, (d+1)*num_pose_ + i, d, 1);
                    file << "VERTEX_XY" <<" " << i <<" "<< t(0) << " " << t(1)
                          << "\n";
                }
            } else
            {
                for (auto i = 0; i < num_pose_; ++i) {
                    Quaternion q(R.block<3,3>(0, num_pose_ + i*3));
                    Vector t = R.block(0, i, d, 1);
                    file << "VERTEX_SE3:QUAT" << " "<< i <<" " << t(0) << " " << t(1)
                         << " " << t(2) << " " << q.x() << " "
                         << q.y() << " " << q.z() << " "
                         << q.w() << "\n";
                }
                for (auto i = 0; i < num_landmark_; ++i)
                {
                    Vector t = R.block(0, (d+1)*num_pose_ + i, d, 1);
                    file << "VERTEX_XYZ" <<" " << i <<" "<< t(0) << " " << t(1) << " " << t(2)
                          << "\n";
                }
            }
        } else {
            std::ofstream file(path + ".txt");
            if (d == 2) {
                for (auto i = 0; i < num_pose_; ++i) {
                    Eigen::Matrix3d R3 = Eigen::Matrix3d::Identity();
                    R3.topLeftCorner<2,2>() =
                        R.block(0, num_pose_ + i*2, 2, 2);
                    auto q = Eigen::Quaterniond(R3);
                    Vector t = R.block(0, i, d, 1);
                    file << i << " " << t(0) << " " << t(1)
                         << " " << q.x() << " " << q.y()
                         << " " << q.z() << " " << q.w() << "\n";
                }
            } else {
                for (auto i = 0; i < num_pose_; ++i) {
                    Quaternion q(R.block<3,3>(0, num_pose_ + i*3));
                    Vector t = R.block(0, i, d, 1);
                    file << i << " " << t(0) << " " << t(1)
                         << " " << t(2) << " " << q.x() << " "
                         << q.y() << " " << q.z() << " "
                         << q.w() << "\n";
                }
            }
            file.close();
        }
    }


//    gtsam::Values InitializationFromOdom(const gtsam::DataParser::Measurement& M,
//                                         size_t P,
//                                         bool useDatasetLandmarks = true,
//                                         bool landmarkMeasIsRelativeToPose = true) {
//        using namespace Eigen;
//
//        if (P < 2)
//            throw std::invalid_argument("InitializationFromOdom: P must satisfy P >= 2 for 2D.");
//
//        // ---- Gather 2D pose edges (already reindexed by DataParser) ----
//        struct Edge2D { size_t i, j; MatrixXd R; VectorXd t; };
//        std::vector<Edge2D> odomEdges;
//        std::set<size_t> allPoseIds;
//
//        for (const auto& m : M.poseMeasurements) {
//            if (m.R.rows() == 2 && m.R.cols() == 2 && m.t.size() == 2) {
//                allPoseIds.insert(m.i);
//                allPoseIds.insert(m.j);
//                // Keep only odometry-like edges: j == i+1 (after DataParser reindexing)
//                if (m.j == m.i + 1) {
//                    odomEdges.push_back(Edge2D{m.i, m.j, m.R, m.t});
//                }
//            }
//        }
//
//        if (odomEdges.empty() && !allPoseIds.empty()) {
//            // No odometry edges found (e.g., only loop closures) — we will still return identity poses.
//        }
//
//        std::sort(odomEdges.begin(), odomEdges.end(),
//                  [](const Edge2D& a, const Edge2D& b){ return a.i < b.i; });
//
//        // ---- Chain odometry: T_j = T_i * T_ij ----
//        struct PoseW { MatrixXd R; VectorXd t; };
//        std::map<size_t, PoseW> poseMap;
//
//        if (!odomEdges.empty()) {
//            const size_t startId = odomEdges.front().i;
//            poseMap[startId] = PoseW{ MatrixXd::Identity(2,2), VectorXd::Zero(2) };
//
//            for (const auto& e : odomEdges) {
//                if (!poseMap.count(e.i)) {
//                    // If chain is broken, start a new one at identity for that i
//                    poseMap[e.i] = PoseW{ MatrixXd::Identity(2,2), VectorXd::Zero(2) };
//                }
//                const PoseW& Pi = poseMap[e.i];
//                PoseW Pj;
//                Pj.R = Pi.R * e.R;
//                Pj.t = Pi.t + Pi.R * e.t;
//                poseMap[e.j] = std::move(Pj);
//            }
//        }
//
//        // Ensure every pose id seen in measurements has an entry (identity if not chained)
//        for (size_t pid : allPoseIds) {
//            if (!poseMap.count(pid)) {
//                poseMap[pid] = PoseW{ MatrixXd::Identity(2,2), VectorXd::Zero(2) };
//            }
//        }
//
//        gtsam::Values initial;
//
//        // ---- Insert lifted poses: key = pose id (matches your randomInitAtLevelP) ----
//        for (const auto& kv : poseMap) {
//            const size_t poseKey = kv.first;
//            const PoseW& Pw = kv.second;
//
//            MatrixXd Ymat = MatrixXd::Zero(static_cast<Index>(P), 2);
//            Ymat.topRows(2) = Pw.R; // Y = [R; 0]
//            auto Y = StiefelManifoldKP::projectToManifold(Ymat);
//
//            VectorXd tP = VectorXd::Zero(static_cast<Index>(P));
//            tP.head(2) = Pw.t;
//
//            initial.insert(poseKey, LiftedPoseDP(Y, tP));
//        }
//
//        // ---- Landmarks: DataParser already mapped landmark IDs to contiguous indices (meas.j) ----
//        std::default_random_engine rng(std::default_random_engine::default_seed);
//        std::normal_distribution<double> gauss(0.0, 1.0);
//
//        // Keep first observation per contiguous landmark id
//        std::set<size_t> seenLmkIdx;
//        for (const auto& lm : M.landmarkMeasurements) {
//            const size_t lidx = lm.j; // contiguous index from DataParser
//            if (seenLmkIdx.count(lidx)) continue;
//            seenLmkIdx.insert(lidx);
//
//            VectorXd lP(static_cast<Index>(P));
//            if (useDatasetLandmarks) {
//                // Compute world point from first observation
//                VectorXd xy(2);
//                if (landmarkMeasIsRelativeToPose) {
//                    // lm.i is the (already reindexed) pose index
//                    auto it = poseMap.find(lm.i);
//                    if (it != poseMap.end()) {
//                        xy = it->second.R * lm.l + it->second.t; // relative -> world
//                    } else {
//                        // Fallback: if the observing pose was not present (shouldn't happen),
//                        // use the measurement vector directly.
//                        xy = lm.l;
//                    }
//                } else {
//                    // Dataset already in world frame
//                    xy = lm.l;
//                }
//                lP.setZero();
//                lP.head(2) = xy;
//            } else {
//                // Random landmark
//                for (Index r = 0; r < lP.size(); ++r) lP[r] = gauss(rng);
//            }
//
//            initial.insert(gtsam::Symbol('L', lidx), lP); // contiguous L0,L1,L2,...
//        }
////        initial.print("Test initials construction:");
//        return initial;
//    }




    /**
     * @brief Initialize lifted poses (and optionally landmarks) from a 2D odometry dataset file.
     *
     * This method parses a `.g2o` dataset containing `VERTEX_SE2`, `EDGE_SE2`, and optionally
     * `VERTEX_XY` entries, reconstructs absolute poses by chaining odometry edges, remaps IDs
     * to contiguous indices, and creates a `gtsam::Values` containing `LiftedPoseDP` variables
     * for poses and vectors for landmarks.
     *
     * @param dataPath             Path to the `.g2o` dataset file.
     * @param P                    Ambient dimension of the lifted representation
     *                              (must satisfy P ≥ 2 for 2D datasets).
     * @param useDatasetLandmarks  If true, initializes landmarks from dataset coordinates;
     *                              if false, initializes them with random Gaussian vectors.
     * @return                     `gtsam::Values` containing the initialized variables.
     *
     * @throws std::invalid_argument if P < 2.
     * @throws std::runtime_error if the dataset file cannot be opened.
     *
     * @note Pose IDs from the dataset are remapped to contiguous indices starting from zero.
     * @note Odometry chaining assumes edges are `EDGE_SE2` measurements in either forward or reverse order.
     */
    gtsam::Values InitializationFromOdom(const std::string& dataPath,
                                         size_t P,
                                         bool useDatasetLandmarks = true)
    {
        if (P < 2)
            throw std::invalid_argument("InitializationFromOdom: P must satisfy P >= 2 for 2D.");

        struct PoseW { Eigen::MatrixXd R; Eigen::VectorXd t; };
        struct Edge2D { size_t i, j; Eigen::MatrixXd R; Eigen::VectorXd t; };
        struct LmkRec { size_t origId; double x; double y; };

        std::map<size_t, PoseW> poseAbs;
        std::vector<Edge2D> odomEdges;
        std::set<size_t> allPoseIds;
        std::vector<LmkRec> lmksInOrder;
        std::unordered_map<size_t, size_t> lmkIdRemap;

        std::ifstream fin(dataPath);
        if (!fin)
            throw std::runtime_error("InitializationFromOdom: cannot open file: " + dataPath);

        std::string line;
        while (std::getline(fin, line)) {
            if (line.empty() || line[0] == '#') continue;

            std::istringstream iss(line);
            std::string tag;
            iss >> tag;

            if (tag == "VERTEX_SE2") {
                size_t id; double x, y, th;
                if (!(iss >> id >> x >> y >> th)) continue;

                Eigen::MatrixXd R = Eigen::MatrixXd::Identity(2, 2);
                const double c = std::cos(th), s = std::sin(th);
                R(0,0) =  c; R(0,1) = -s;
                R(1,0) =  s; R(1,1) =  c;

                Eigen::VectorXd t = Eigen::VectorXd::Zero(2);
                t(0) = x; t(1) = y;

                poseAbs[id] = PoseW{R, t};
                allPoseIds.insert(id);
            }
            else if (tag == "EDGE_SE2") {
                size_t i, j; double dx, dy, dth;
                double I11, I12, I13, I22, I23, I33;
                if (!(iss >> i >> j >> dx >> dy >> dth >> I11 >> I12 >> I13 >> I22 >> I23 >> I33))
                    continue;

                Eigen::MatrixXd Rij = Eigen::MatrixXd::Identity(2,2);
                const double c = std::cos(dth), s = std::sin(dth);
                Rij(0,0)= c; Rij(0,1)= -s; Rij(1,0)= s; Rij(1,1)= c;

                Eigen::VectorXd tij(2);
                tij << dx, dy;

                odomEdges.push_back(Edge2D{i, j, Rij, tij});
                allPoseIds.insert(i);
                allPoseIds.insert(j);
            }
            else if (tag == "VERTEX_XY") {
                size_t oid; double x, y;
                if (!(iss >> oid >> x >> y)) continue;

                if (!lmkIdRemap.count(oid)) {
                    lmkIdRemap[oid] = lmksInOrder.size();
                    lmksInOrder.push_back({oid, x, y});
                }
            }
        }
        fin.close();

        // --- Build world poses from odometry ---
        std::map<size_t, PoseW> poseW;

        for (const auto& kv : poseAbs) poseW[kv.first] = kv.second;

        if (poseW.empty() && !allPoseIds.empty()) {
            const size_t seed = *allPoseIds.begin();
            poseW[seed] = PoseW{Eigen::MatrixXd::Identity(2,2), Eigen::VectorXd::Zero(2)};
        }

        bool progressed = true;
        while (progressed) {
            progressed = false;

            for (const auto& e : odomEdges) {
                if (poseW.count(e.i) && !poseW.count(e.j)) {
                    const PoseW& Pi = poseW[e.i];
                    PoseW Pj;
                    Pj.R = Pi.R * e.R;
                    Pj.t = Pi.t + Pi.R * e.t;
                    poseW[e.j] = std::move(Pj);
                    progressed = true;
                } else if (!poseW.count(e.i) && poseW.count(e.j)) {
                    const PoseW& Pj = poseW[e.j];
                    Eigen::MatrixXd Rij_inv = e.R.transpose();
                    Eigen::VectorXd tij_inv = -Rij_inv * e.t;
                    PoseW Pi;
                    Pi.R = Pj.R * Rij_inv;
                    Pi.t = Pj.t + Pj.R * tij_inv;
                    poseW[e.i] = std::move(Pi);
                    progressed = true;
                }
            }

            if (!progressed) {
                for (size_t pid : allPoseIds) {
                    if (!poseW.count(pid)) {
                        poseW[pid] = PoseW{Eigen::MatrixXd::Identity(2,2), Eigen::VectorXd::Zero(2)};
                        progressed = true;
                        break;
                    }
                }
            }
        }

        for (size_t pid : allPoseIds) {
            if (!poseW.count(pid)) {
                poseW[pid] = PoseW{Eigen::MatrixXd::Identity(2,2), Eigen::VectorXd::Zero(2)};
            }
        }

        // --- Remap pose IDs to contiguous ---
        std::vector<size_t> poseOrigIds;
        poseOrigIds.reserve(poseW.size());
        for (const auto& kv : poseW) poseOrigIds.push_back(kv.first);
        std::sort(poseOrigIds.begin(), poseOrigIds.end());

        std::unordered_map<size_t, size_t> poseIdRemap;
        for (size_t k = 0; k < poseOrigIds.size(); ++k) {
            poseIdRemap[poseOrigIds[k]] = k;
        }

        // --- Build initial Values ---
        gtsam::Values initial;

        for (const auto& kv : poseW) {
            const size_t poseKeyContig = poseIdRemap.at(kv.first);
            const PoseW& Pw = kv.second;

            Eigen::MatrixXd Ymat = Eigen::MatrixXd::Zero(static_cast<Eigen::Index>(P), 2);
            Ymat.topRows(2) = Pw.R;
            auto Y = StiefelManifoldKP::projectToManifold(Ymat); // keep original usage

            Eigen::VectorXd tP = Eigen::VectorXd::Zero(static_cast<Eigen::Index>(P));
            tP.head(2) = Pw.t;

            initial.insert(poseKeyContig, LiftedPoseDP(Y, tP));
        }

        std::default_random_engine rng(std::default_random_engine::default_seed);
        std::normal_distribution<double> gauss(0.0, 1.0);

        for (size_t contiguousIdx = 0; contiguousIdx < lmksInOrder.size(); ++contiguousIdx) {
            const auto& rec = lmksInOrder[contiguousIdx];

            Eigen::VectorXd lP(static_cast<Eigen::Index>(P));
            if (useDatasetLandmarks) {
                lP.setZero();
                lP.head(2) << rec.x, rec.y;
            } else {
                for (Eigen::Index r = 0; r < lP.size(); ++r) lP[r] = gauss(rng);
            }

            initial.insert(gtsam::Symbol('L', contiguousIdx), lP);
        }

        return initial;
    }


    /**
     * @brief Perform local search-based optimization (using odometry-based initialization).
     *
     * This method builds a factor graph at the given rank, initializes it with
     * the current odometry-based values, and runs a Levenberg–Marquardt optimizer.
     * Optimization results (SDP value, elapsed time) are recorded in
     * `certificateResults_`.
     *
     * @param p  Rank (ambient dimension) of the lifted representation.
     * @return   Always returns std::nullopt; optimization results are stored in
     *           `certificateResults_` as a side effect.
     */
    std::optional<CertificateResults> LocalSearchWithOdomInitials(size_t p)
    {

        Values Qstar;
        auto t6 = CFGStopwatch::tick();
        std::cout << "Starting Local search-based optimization "<< std::endl;
        NonlinearFactorGraph graph = buildGraphAtLevel(p);
        auto lmParams = opts_.lmParams;
        lmParams.maxIterations    = opts_.maxIterations;
        lmParams.relativeErrorTol = opts_.relativeErrorTol;
        lmParams.absoluteErrorTol = opts_.absoluteErrorTol;
        lmParams.verbosityLM      = opts_.verbosityLM;
        auto lm = std::make_shared<LevenbergMarquardtOptimizer>(graph, currentValues_, lmParams);
        auto t0 = CFGStopwatch::tick();
        auto result = lm->optimize();
        certificateResults_.SDPval.push_back(lm->error());
        auto t1 = CFGStopwatch::tock(t0);
        certificateResults_.elapsed_optimization_times.push_back(t1);
        certificateResults_.total_computation_time = CFGStopwatch::tock(t6);

        std::cout << "Finish Local search-based optimization. " << std::endl;
        return certificateResults_;

    }






};


using CertifiableLandmark2 = CertifiableLandmark<2>;
using CertifiableLandmark3 = CertifiableLandmark<3>;


#endif //CERTIFIABLELANDMARK_H
