//
// Created by jason on 5/31/25.
//

#ifndef STIEFELMANIFOLDEXAMPLE_CERTIFIABLERANGESLAM_H
#define STIEFELMANIFOLDEXAMPLE_CERTIFIABLERANGESLAM_H

#include "Certifiable_problem.h"
#include "LiftedRangeFactor.h"
#include "NormalRangeFactor.h"
/**
 * @brief
 *
 *
 * @tparam d  Ambient pose dimension (must be 2 or 3).
 */
template <size_t d>
class CertifiableRangeSLAM : public CertifiableProblem {
    static_assert(d == 2 || d == 3, "CertifiableRangeSLAM only supports d = 2 or 3.");

    /// Number of landmarks in the problem.
    size_t num_landmark_;

    /// Number of pose‐to‐pose measurements.
    size_t num_pose_measurements_;

    // Number of range measurement.
    size_t num_range_;

public:
    /**
     * @brief Construct with initial rank and measurement data.
     * @param p             Initial relaxation rank.
     * @param measurements  Parsed measurement struct (num_poses, etc.).
     */
    CertifiableRangeSLAM(size_t p, const DataParser::Measurement& measurements, const std::string & g2oPath)
            : CertifiableProblem(d, p, measurements)
    {
        num_landmark_ = measurements.num_landmarks;
        num_range_ = measurements.num_ranges;
        certificateResults_.startingRank = p;
        dataPath = g2oPath;
    }

    /**
     * @brief Initialize graph, data matrix, and random values; record init time.
     */
    void init() {
        num_pose_measurements_ = measurements_.poseMeasurements.size();
        auto t0 = CFGStopwatch::tick();
        if (opts_.initType == CertifiableProblemOpts::InitType::Random){
            currentValues_  = randomInitAtLevelP(currentRank_);
        }
        else if (opts_.initType == CertifiableProblemOpts::InitType::Odom){
            currentValues_  = InitializationFromOdom(dataPath, currentRank_);
        }
        else if (opts_.initType == CertifiableProblemOpts::InitType::LocalSearch){
            currentValues_  = InitializationFromOdom_NonLifted(dataPath, currentRank_);
        }
        else {
            throw std::runtime_error("Unknown initialization type in CertifiableProblemOpts::initType");
        }

        if (opts_.initType == CertifiableProblemOpts::InitType::LocalSearch){
            currentGraph_   = buildGraphAtLevel_NonLifted(currentRank_);
        }
        else {
            currentGraph_   = buildGraphAtLevel(currentRank_);
        }
        //M_              = recoverDataMatrixFromHessian();
        //std::cout<<(recoverDataMatrixFromGraph()-recoverDataMatrixFromHessian()).norm()<<std::endl;
        M_ = recoverDataMatrixFromGraph();
        auto t1 = CFGStopwatch::tock(t0);
        certificateResults_.initialization_time.push_back(t1);
    }

    /**
     * @brief Assemble the factor graph at relaxation level p.
     * @param p  Relaxation rank.
     * @return   Factor graph at level-p
     */
    NonlinearFactorGraph buildGraphAtLevel(size_t p) override {
        NonlinearFactorGraph inputGraph;
        size_t k = 0;
        for (const auto meas : measurements_.rangeMeasurements)
        {
            Vector sigmas = Vector::Zero(p);
            Vector sigmas_MLE = Vector::Zero(d);
            sigmas.head(p).setConstant(sqrt(1 * meas.sigma));
            sigmas_MLE.setConstant(sqrt(1 * meas.sigma));
            auto noise = noiseModel::Diagonal::Sigmas(sigmas);
            auto noise_MLE = noiseModel::Diagonal::Sigmas(sigmas_MLE);
            if (d == 2) {
                inputGraph.emplace_shared<LiftedRangeFactor2>( meas.i, Symbol('L', meas.j), Symbol('R', k), meas.range, p , noise);
            }
            else if (d == 3)
            {
                inputGraph.emplace_shared<LiftedRangeFactor3>( meas.i, Symbol('L', meas.j), Symbol('R', k), meas.range, p , noise);
            }
            k++;
        } // range factor
        for (const auto meas : measurements_.poseMeasurements)
        {
            Vector sigmas = Vector::Zero(p * d + p);
            sigmas.head(p * d).setConstant(std::sqrt(1.0 / (1 * meas.kappa)));
            sigmas.tail(p).setConstant(std::sqrt(1.0 / (1 * meas.tau)));
            auto noise = noiseModel::Diagonal::Sigmas(sigmas);

            if (d == 2) {
                double theta = std::atan2(meas.R(1, 0), meas.R(0, 0)); // Extract angle
                auto odom = gtsam::Pose2(meas.t(0), meas.t(1), theta);
                inputGraph.emplace_shared<SEsyncFactor2>(meas.i,  meas.j, meas.R, meas.t, p, noise);
            }
            else if (d ==3 ) {
                gtsam::Rot3 rot(meas.R);
                gtsam::Point3 trans(meas.t(0), meas.t(1), meas.t(2));
                auto odom = gtsam::Pose3(rot, trans);
                inputGraph.emplace_shared<SEsyncFactor3>(meas.i, meas.j, meas.R, meas.t, p, noise);
            }
            else {
                std::cerr << "Un" << std::endl;
            }
        } // PGO Factor
        return inputGraph;
    }

    // Not used..., but need to instantiate it.
    Matrix computeLambdaBlocks(const Matrix& Y) override {
        return Matrix();
    }

    // need to create incremental version
    NonlinearFactorGraph buildGraphAtLevelUsingGraph(size_t p) override {
        NonlinearFactorGraph graph;
        return graph;
    }
    /**
     * @brief Compute the dense block‐diagonal certificate matrix.
     * @param Y  Variable matrix.
     * @return   Dense block matrix of size
     */
    std::pair<Matrix, Vector> computeLambdaBlocksRangeSLAM(const Matrix& Y) {
        Matrix QY = M_ * Y;
        Matrix Yt = Y.transpose();
        // Preallocate storage for diagonal blocks of Lambda
        Matrix stiefel_Lambda_blocks(d, num_pose_ * d);

        // Index of the row/column at which the rotational blocks begin in matrix X
        size_t offset = num_pose_;

        for (size_t i = 0; i < num_pose_; ++i) {
            Matrix P = QY.block(offset + i * d, 0, d, Y.cols()) *
                       Yt.block(0, offset + i * d, Y.cols(), d);
            stiefel_Lambda_blocks.block(0, i * d, d, d) = .5 * (P + P.transpose());
        }


        Vector oblique_Lambda_blocks(num_range_);
        auto offset_range = (d+1) * num_pose_ + num_landmark_;
        for (size_t i = 0; i < num_range_; ++i) {
            Matrix P = QY.block(offset_range + i, 0, 1, Y.cols()) *
                       Yt.block(0, offset_range + i, Y.cols(), 1);
            oblique_Lambda_blocks.block(i, 0, 1, 1) = .5 * (P + P.transpose());
        }

        return std::make_pair(stiefel_Lambda_blocks, oblique_Lambda_blocks);
    }


    // Not used..., but need to instantiate it.
    SparseMatrix computeLambdaFromLambdaBlocks(const Matrix& LambdaBlocks) override {
        return SparseMatrix();
    }

    /**
     * @brief Convert dense Λ_blocks into a sparse Λ matrix for certification.
     * @param LambdaBlocks  Dense block matrix.
     * @return              Sparse certificate matrix Λ.
     */
    SparseMatrix computeLambdaFromLambdaBlocksRangeSLAM
            (const std::pair<Matrix, Vector> &Lambda_blocks) {
        std::vector<Eigen::Triplet<Scalar>> elements;
        elements.reserve(d * num_pose_ + num_range_);

        // add the symmetric diagonal blocks for the Stiefel constraints
        for (auto i = 0; i < num_pose_; ++i) { // block index
            for (auto r = 0; r < d; ++r) {     // block row index
                for (auto c = 0; c < d; ++c) {   // block column index
                    elements.emplace_back(num_pose_ + i * d + r, num_pose_ + i * d + c,
                                          Lambda_blocks.first(r, i * d + c));
                }
            }
        }

        auto rot_mat_sz = num_pose_ * (d+1) + num_landmark_;
        // add the diagonal block for the Oblique constraints
        for (auto i = 0; i < num_range_; ++i) {
            elements.emplace_back(rot_mat_sz + i, rot_mat_sz + i,
                                  Lambda_blocks.second(i));
        }

        // add additional zeros if we're using the explicit formulation
        int Lambda_size = (d+1) * num_pose_ + num_range_ + num_landmark_;

        SparseMatrix Lambda(Lambda_size, Lambda_size);
        Lambda.setFromTriplets(elements.begin(), elements.end());
        return Lambda;
    }

    /**
     * @brief Build the element matrix S from current Values.
     * @param values  GTSAM Values containing LiftedPoseDP variables.
     * @return        Sparse matrix S of size.
     */
    SparseMatrix elementMatrix(const Values& values) override {
        using Index     = long;
        using Triplet64 = Eigen::Triplet<double, Index>;

        // ---------- 1. counts ----------------------------------------------------
        const std::size_t N_pose = values.count<LiftedPoseDP>();
        const std::size_t N_r    = values.count<UnitSphereD>();
        const std::size_t N_lmk  = values.count<Vector>();

        if (N_pose == 0) throw std::runtime_error("At least one LiftedPoseDP required.");
        const auto poses = values.extract<LiftedPoseDP>();
        if (poses.empty()) throw std::runtime_error("No LiftedPoseDP in Values.");

        const std::size_t p = poses.begin()->second.get_Rows();

        // ---------- 2. row offsets (new order) -----------------------------------
        const std::size_t off_t_pose = 0;
        const std::size_t off_r_pose = off_t_pose + N_pose;
        const std::size_t off_t_lmk  = off_r_pose + N_pose * d;
        const std::size_t off_r      = off_t_lmk + N_lmk;
        const std::size_t nrows      = off_r + N_r;

        // ---------- 3. triplets --------------------------------------------------
        std::vector<Triplet64> T;
        T.reserve(static_cast<std::size_t>(p) * nrows);

        // ---------- 4-a. t_pose (translation part of LiftedPoseDP) --------------
        for (const auto& kv : poses) {
            const auto& pose = kv.second;
            const size_t idx = gtsam::symbolIndex(kv.first);
            const Eigen::VectorXd& t = pose.get_TranslationVector();

            const Index rowVar = static_cast<Index>(off_t_pose + idx);
            for (std::size_t r = 0; r < p; ++r)
                T.emplace_back(rowVar, static_cast<Index>(r), t(static_cast<int>(r)));
        }

        // ---------- 4-b. r_pose (Stiefel part of LiftedPoseDP) ------------------
        for (const auto& kv : poses) {
            const auto& pose = kv.second;
            const size_t idx = gtsam::symbolIndex(kv.first);
            Eigen::MatrixXd R = pose.matrix(); // p × d

            const std::size_t row0 = off_r_pose + idx * d;
            for (std::size_t c = 0; c < d; ++c) {
                const Index rowVar = static_cast<Index>(row0 + c);
                for (std::size_t r = 0; r < p; ++r)
                    T.emplace_back(rowVar, static_cast<Index>(r), R(static_cast<int>(r), static_cast<int>(c)));
            }

        }

        // ---------- 4-c. t_lmk (landmark translation vectors) --------------------
        for (const auto& kv : values.extract<Vector>()) {
            const auto& lmk = kv.second;
            const size_t idx = gtsam::symbolIndex(kv.first);
            const Eigen::VectorXd& t = lmk;

            const Index rowVar = static_cast<Index>(off_t_lmk + idx);
            for (std::size_t r = 0; r < p; ++r)
                T.emplace_back(rowVar, static_cast<Index>(r), t(static_cast<int>(r)));
        }

        // ---------- 4-d. r (UnitSphereD) ----------------------------------------
        for (const auto& kv : values.extract<UnitSphereD>()) {
            const auto& unitSphere = kv.second;
            const size_t idx = gtsam::symbolIndex(kv.first);
            const Eigen::VectorXd& v = unitSphere.matrix();

            const Index rowVar = static_cast<Index>(off_r + idx);
            for (std::size_t r = 0; r < p; ++r)
                T.emplace_back(rowVar, static_cast<Index>(r), v(static_cast<int>(r)));
        }

        // ---------- 5. build sparse ---------------------------------------------
        SparseMatrix S(static_cast<Index>(nrows), static_cast<Index>(p));
        S.setFromTriplets(T.begin(), T.end());
        return S;
    }

    /**
     * @brief Randomly initialize all poses at level Pmin using uniform Stiefel and random translation.
     * @param Pmin  Target relaxation rank.
     * @return      GTSAM Values container with random LiftedPoseDP entries.
     */
    Values randomInitAtLevelP(const size_t Pmin) override {
        Values initial;
        // Add poses
        for (size_t j = 0; j < num_pose_; j++) {
            StiefelManifoldKP Y = StiefelManifoldKP::Random(std::default_random_engine::default_seed, d, Pmin);
            Vector trans = Vector::Random(Pmin);
            initial.insert( j, LiftedPoseDP(Y, trans));
        }



        // Add landmarks
        for (size_t j = 0; j < num_landmark_; j++)
        {
            Vector trans1 = Vector::Random(Pmin);
            initial.insert(Symbol('L', j), trans1);
        }

        // Add landmarks range spheres
        for (size_t j = 0; j < num_range_; j++)
        {
            Matrix X0 = Matrix::Random(Pmin, 1).normalized();
            UnitSphereD Y = UnitSphereD(X0);
            initial.insert(Symbol('R', j), Y);
        }

        return initial;
    }

    /**
     * @brief Convert a flat eigenvector into tangent VectorValues.
     * @param p  Relaxation rank.
     * @param v  Flattened direction of size.
     * @param values  Lifted Values at level p.
     * @return  VectorValues on each local tangent space.
     */
    VectorValues TangentVectorValues(
            size_t p, const Vector v, const Values values) override
    {
        VectorValues delta;
        // 1. Construct full tangent matrix using descent direction
        Matrix Ydot = Matrix::Zero(v.size(), p);  // Each row: p-dimensional vector
        Ydot.rightCols<1>() = v;

        // 2. Count variables
        const std::size_t N_pose = values.count<LiftedPoseDP>();
        const std::size_t N_lmk  = values.count<Vector>();
        const std::size_t off_t_pose = 0;
        const std::size_t off_r_pose = off_t_pose + N_pose;
        const std::size_t off_t_lmk  = off_r_pose + N_pose * d;
        const std::size_t off_r      = off_t_lmk + N_lmk;

        const auto poses = values.extract<LiftedPoseDP>();
        for (const auto& kv : poses) {
            const auto& pose = kv.second;
            StiefelManifoldKP Y = pose.get_Y();
            const size_t idx = gtsam::symbolIndex(kv.first);

            Matrix tangentMatrix = Ydot.block(N_pose + idx * d, 0, d, p).transpose();
            Vector transV = Ydot.block(idx, 0, 1, p).transpose();

            Vector xi = StiefelManifoldKP::Vectorize(tangentMatrix);
            Vector tangentVector = Y.G_.transpose() * xi;

            Vector combined(tangentVector.size() + transV.size());
            combined << tangentVector, transV;

            delta.insert(idx, combined);
        }

        // --- c. t_lmk (landmarks)
        for (const auto& kv : values.extract<Vector>()) {
            const size_t idx = gtsam::symbolIndex(kv.first);
            Vector xi = Ydot.block(off_t_lmk + idx, 0, 1, p).transpose();
            delta.insert(N_pose + idx, xi);
        }

        // --- d. r (unit spheres)
        for (const auto& kv : values.extract<UnitSphereD>()) {
            const size_t idx = gtsam::symbolIndex(kv.first);
            UnitSphereD Y = kv.second;
            Vector xi = Ydot.block(off_r + idx, 0, 1, p).transpose();
            Vector tangentVector = Y.G_.transpose() * xi;
            delta.insert(N_pose + N_lmk + idx, tangentVector);
        }
        return delta;
    }

    /**
     * @brief Project ambient‐space variation Ydot onto the tangent space at Y.
     * @param p     Relaxation rank.
     * @param Y     Basepoint matrix.
     * @param Ydot  Ambient variation.
     * @return      Tangent‐space projection.
     */
    Matrix tangent_space_projection(
            const size_t p, const Matrix& Y, const Matrix& Ydot) override
    {
        const std::size_t N_pose = num_pose_;
        const std::size_t N_r    = num_range_;
        const std::size_t N_lmk  = num_landmark_;

        Matrix result = Ydot;
        size_t offset_t = N_pose;
        size_t offset_lmk = N_pose * (d + 1) + N_lmk;
        // Stiefel component
        auto rot_mat_sz = d * N_pose;
        result.block(offset_t, 0, rot_mat_sz, p) =
                StiefelManifoldKP::Proj(Y.block(offset_t, 0, rot_mat_sz, p).transpose(),
                                        result.block(offset_t, 0, rot_mat_sz, p).transpose()).transpose();

        // Oblique component
        int r = N_r;
        result.block(offset_lmk, 0, r, p) =
                StiefelManifoldKP::Proj(Y.block(offset_lmk, 0, r, p).transpose(),
                                        result.block(offset_lmk, 0, r, p).transpose()).transpose();
        // remaining component is untouched
        return result;
    }

    /**
     * @brief Recover the data matrix from the current factor graph.
     * @return Sparse data matrix L of size.
     */
    SparseMatrix recoverDataMatrixFromGraph() override {
        std::vector<Eigen::Triplet<Scalar>> triplets;
        triplets.reserve(num_pose_measurements_ * (4 + 2*d*d + 2*d + 2) + 9 * num_range_);
        if (d ==2)
        {
            for (const auto& f_ptr : currentGraph_) {
                if (auto factor = std::dynamic_pointer_cast<SEsyncFactor2>(f_ptr)) {
                    factor->appendBlocksFromFactor(num_pose_, triplets);
                    continue;
                }

                if (auto factor = std::dynamic_pointer_cast<LiftedRangeFactor2>(f_ptr)) {
                    factor->appendBlocksFromFactor(num_pose_, num_landmark_, triplets);
                    continue;
                }
            }
        }
        else
        {
            for (const auto& f_ptr : currentGraph_) {
                if (auto factor = std::dynamic_pointer_cast<SEsyncFactor3>(f_ptr)) {
                    factor->appendBlocksFromFactor(num_pose_, triplets);
                    continue;
                }

                if (auto factor = std::dynamic_pointer_cast<LiftedRangeFactor3>(f_ptr)) {
                    factor->appendBlocksFromFactor(num_pose_, num_landmark_, triplets);
                    continue;
                }
            }
        }
        SparseMatrix L((d+1) * num_pose_ + num_landmark_ + num_range_, (d+1) * num_pose_ + num_landmark_ + num_range_);
        L.setFromTriplets(triplets.begin(), triplets.end());
        L.makeCompressed();
        return L;
    }
     SparseMatrix recoverDataMatrixFromHessian() {
        std::map<std::pair<Key,Key>,Matrix> HMap;
        std::map<std::tuple<Key,Key,Key>,Matrix> HMap_range;
        if (d == 2) {
            for (auto& f_ptr : currentGraph_) {
                if (auto factor = std::dynamic_pointer_cast<SEsyncFactor2>(f_ptr))
                {
                    factor->computeHessian(HMap);
                }
                if (auto factor = std::dynamic_pointer_cast<LiftedRangeFactor2>(f_ptr))
                {
                    factor->computeHessian(HMap_range);
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
                if (auto factor = std::dynamic_pointer_cast<LiftedRangeFactor3>(f_ptr))
                {
                    factor->computeHessian(HMap_range);
                }
            }
        }
        Ordering order = Ordering::Natural(currentGraph_);
        // Build the index in one shot:
        const size_t N = order.size();
        std::vector<size_t> dims(N), offsets(N);
        for (size_t i = 0; i < N; ++i) {
            const Key k = order[i];
            // look at the letter in the key:
            if (symbolChr(k)=='L') {
                dims[i] = 1;
            } else if (symbolChr(k)=='R'){
                dims[i] = 1;
            }
            else {
                auto lp = currentValues_.at<LiftedPoseDP>(k);
                dims[i] = lp.get_Cols();
            }
            offsets[i] = (i==0 ? 0u : offsets[i-1] + dims[i-1]);
        }
        size_t totalDim = offsets.back() + dims.back();
        // 2) Gather triplets
        std::vector<Eigen::Triplet<Scalar>> triplets;
        triplets.reserve(num_pose_measurements_ * (4 + 2*d*d + 2*d + 2) + 9 * num_range_);
        appendPGOHessian(triplets,HMap);
        appendRangeHessian(triplets,HMap_range,offsets,order);
        SparseMatrix L(totalDim, totalDim);
        L.setFromTriplets(triplets.begin(), triplets.end());
        L.makeCompressed();
        return L;
    }

    void appendPGOHessian(std::vector<Eigen::Triplet<Scalar>> &triplets, std::map<std::pair<Key,Key>,Matrix> &HMap)
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
    void appendRangeHessian(std::vector<Eigen::Triplet<Scalar>> &triplets, std::map<std::tuple<Key,Key,Key>,Matrix> &HMap_range, std::vector<size_t> offsets, Ordering order)
    {
        for (auto const& kv : HMap_range) {
            const Key& k1 = std::get<0>(kv.first); // Pose
            const Key& k2 = std::get<1>(kv.first); // Landmark
            const Key& k3 = std::get<2>(kv.first); // Unit Sphere
            size_t j;
            size_t k;
            for (size_t ii = 0; ii < order.size(); ++ii)
            {
                Key key = order[ii];
                if (symbolChr(key)=='L' && symbolIndex(key)==symbolIndex(k2))
                {
                    j = ii;
                    break;
                }
            }
            for (size_t ii = 0; ii < order.size(); ++ii)
            {
                Key key = order[ii];
                if (symbolChr(key)=='R' && symbolIndex(key)==symbolIndex(k3))
                {
                    k = ii;
                    break;
                }
            }
            size_t i = symbolIndex(k1);
            Matrix Hfull = kv.second;
            // Htiti
            triplets.emplace_back(i,i, Hfull(0,0));
            // Htili
            triplets.emplace_back(i,offsets[j], Hfull(0,1));
            // Hliti
            triplets.emplace_back(offsets[j],i, Hfull(1,0));
            //Hlili
            triplets.emplace_back(offsets[j],offsets[j], Hfull(1,1));
            //Hti r
            triplets.emplace_back(i,offsets[k], Hfull(0,2));
            //Hr ti
            triplets.emplace_back(offsets[k],i, Hfull(2,0));
            //Hli r
            triplets.emplace_back(offsets[j],offsets[k], Hfull(1,2));
            //Hr li
            triplets.emplace_back(offsets[k],offsets[j], Hfull(2,1));
            //Hr r
            triplets.emplace_back(offsets[k],offsets[k], Hfull(2,2));
        }
    }

    Values exportInitialsFromSDPSolution(const Matrix& R)
    {
        Values initialsRefine;
        if (d == 3) {
            for (auto i = 0; i < num_pose_; ++i) {
                Rot3 rot3 = Rot3(R.block(0, num_pose_ + i * d, d, d));
                Vector3 trans = R.block(0, i, d, 1);
                initialsRefine.insert(i, LiftedPoseDP(StiefelManifoldKP(rot3.matrix()), trans));
            }

            for (auto j = 0; j < num_landmark_; ++j) {
                initialsRefine.insert(
                        Symbol('L', j),
                        Vector(R.block(0, num_pose_ * (d+1) + j, d, 1))
                );
            }

            for (int k = 0; k < num_range_; ++k) {
                initialsRefine.insert(Symbol('R', k), UnitSphereD(R.block(0, num_pose_ * (d+1) + num_landmark_ + k, d, 1)));
            }
        } else {
            for (auto i = 0; i < num_pose_; ++i) {
                Rot2 rot2 = Rot2::fromCosSin(
                        R.block(0, num_pose_ + i * d, d, d)(0, 0),
                        R.block(0, num_pose_ + i * d, d, d)(1, 0)
                );
                Vector2 trans = R.block(0, i, d, 1);
                initialsRefine.insert(i, LiftedPoseDP(StiefelManifoldKP(rot2.matrix()), trans));
            }

            for (auto j = 0; j < num_landmark_; ++j) {
                initialsRefine.insert(
                        Symbol('L', j),
                        Vector(R.block(0, num_pose_ * (d+1) + j, d, 1))
                );
            }

            for (int k = 0; k < num_range_; ++k) {
                initialsRefine.insert(Symbol('R', k), UnitSphereD(R.block(0, num_pose_ * (d+1) + num_landmark_ + k, d, 1)));
            }
        }

        return initialsRefine;
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
        return 0.5 * (Y.transpose() * dataMatrixProduct(Y, M)).trace();
    }

    /**
     * @brief Search for a certifiable solution between pMin and pMax.
     *
     * Performs LM optimization, gradient checks, and fast_verification at each level.
     * Lifts along descent if verification fails. Returns certificate results on success.
     *
     * @param pMin  Minimum relaxation rank.
     * @param pMax  Maximum relaxation rank.
     * @return      CertificateResults on success, or std::nullopt if none found.
     */
    std::optional<CertificateResults> Solve(size_t pMin, size_t pMax) {
        Values Qstar;
        auto t6 = CFGStopwatch::tick();
        for (size_t p = pMin; p <= pMax; ++p) {
            std::cout << "Starting optimization at rank = " << p << std::endl;
            setCurrentRank(p);
            Qstar = tryOptimizingAtLevel(p);
            setCurrentValues(Qstar);
            auto t2 = std::chrono::high_resolution_clock::now();
            auto nonlinear_graph = buildGraphAtLevel(p);
            auto linear_graph = nonlinear_graph.linearize(Qstar);
            auto grad_norm = linear_graph->gradientAtZero();
            std::cout << "Gradient norm at level p = " << p << " is : " << grad_norm.norm() << std::endl;
            std::cout << "Objective value is : " << std::scientific << std::setprecision(6) << nonlinear_graph.error(Qstar) << std::endl;

            certificateResults_.gradnorm.push_back(grad_norm.norm());

            auto t4 = CFGStopwatch::tick();
            SparseMatrix S = elementMatrix(Qstar);
            std::pair<Matrix, Vector> lambdaBlocks = computeLambdaBlocksRangeSLAM(S);
            SparseMatrix Lambda = computeLambdaFromLambdaBlocksRangeSLAM(lambdaBlocks);
            SparseMatrix M = getDataMatrix();
            Scalar obj = evaluateObjective(S, M);

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
                std::cout << "Solution verified at level p = " << p << std::endl;
                certificateResults_.Yopt = S;
                certificateResults_.Lambda = Lambda;
                Matrix R = RoundSolutionS();
                certificateResults_.xhat = R;
                auto t7 = CFGStopwatch::tock(t6);
                certificateResults_.total_computation_time = t7;
                certificateResults_.endingRank = p;
                NonlinearFactorGraph graphRefine = buildGraphAtLevel(d);
                currentValues_ = exportInitialsFromSDPSolution(R);
                Values solutionRefine = tryOptimizingAtLevel(d);
                auto linearGraphRefine = graphRefine.linearize(solutionRefine);
                auto gradNormRefine = linearGraphRefine->gradientAtZero();
                std::cout << "Final Gradient norm is : " << std::scientific << std::setprecision(6) << gradNormRefine.norm() << std::endl;
                std::cout << "Final objective error : " << std::scientific << std::setprecision(6) << graphRefine.error(solutionRefine) <<std::endl;
                certificateResults_.gradnorm.push_back(gradNormRefine.norm());
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
    Matrix RoundSolutionS() override {
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

        // Project each spherical variable to the unit sphere by normalizing
        // the respective rows from (rot_mat_sz + 1) to (rot_mat_sz + r)

//    R.block(0, (d +1)* num_pose_ + num_landmarks, d, d).colwise().normalize();
        R.block(0, (d +1)* num_pose_ + num_landmark_, d, num_range_).colwise().normalize();

        return R;
    }

    /**
     * @brief Export the solution in G2O or TUM format.
     * @param path  Base filename (without extension).
     * @param R     Rotation/translation solution matrix.
     * @param g2o   If true, write G2O; otherwise, write TUM.
     */
    void ExportData(const string &path, const Eigen::MatrixXd &R, bool g2o) override {
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
            } else {
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
            file.close();
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


    /**
     * @brief Initialize a set of lifted pose and landmark variables from an odometry-and-range dataset.
     *
     * This function parses a raw dataset file containing odometry edges (`EDGE_SE2`),
     * range measurements (`EDGE_RANGE`), and ground-truth landmark positions
     * (`VERTEX_XY` or `VERTEX_XYZ`). It:
     *  - Reads and stores odometry edges for each robot (poses identified by `<robotChar><id>`).
     *  - Reindexes robot poses to global contiguous IDs.
     *  - Selects only consecutive odometry edges (j = i+1) to chain an initial trajectory.
     *  - Builds an initial guess for each robot pose as a lifted representation
     *    `LiftedPoseDP(Y, t)` where `Y` is the rotation block projected to the Stiefel manifold
     *    and `t` is the translation block.
     *  - Reads ground-truth landmark coordinates from `VERTEX_XY` (2D) and `VERTEX_XYZ` (3D)
     *    entries, merges them, sorts by landmark ID, and initializes contiguous landmark keys
     *    `L0, L1, ...` with the given positions in the first `d` coordinates of the lifted vector.
     *  - Initializes range-measurement "unit sphere" variables (`R0`, `R1`, ...) as
     *    random unit vectors (one per EDGE_RANGE measurement).
     *
     * @param dataPath Path to the dataset file.
     * @param P Lifting dimension (must satisfy `P >= d`).
     * @return gtsam::Values Initialized variable assignments for all poses, landmarks, and spheres.
     *
     * @throws std::invalid_argument If `P < d`.
     * @throws std::runtime_error If the dataset file cannot be opened.
     *
     * @note
     *   - `d` is the ambient spatial dimension (typically 2 or 3).
     *   - Odometry edges are assumed to be between consecutive poses only (`j = i + 1`).
     *   - Landmark coordinates beyond `d` dimensions are truncated.
     *   - For `VERTEX_XYZ` entries, only the first `min(d, 3)` coordinates are used.
     *   - Requires `StiefelManifoldKP::projectToManifold()` and `LiftedPoseDP` types to be defined.
     */
    gtsam::Values InitializationFromOdom(const std::string& dataPath, size_t P) {
        using namespace Eigen;

        if (P < static_cast<size_t>(d))
            throw std::invalid_argument("InitializationFromOdomRange: P must satisfy P >= d.");

        // --------- Parse the file (EDGE_SE2, EDGE_RANGE, VERTEX_XY, VERTEX_XYZ) ---------
        struct RawOdom { char ca, cb; size_t ia, ib; double dx, dy, dtheta; };
        std::vector<RawOdom> rawOdom;
        size_t numRanges = 0;

        // For reindexing robot poses to global contiguous ids:
        std::set<char> robots;
        std::map<char, std::set<size_t>> robotPoseIds;

        // Ground-truth landmarks
        std::map<size_t, Vector2d> vertexXY_GT;   // raw local lmk id -> (x,y)
        std::map<size_t, Vector3d> vertexXYZ_GT;  // raw local lmk id -> (x,y,z)

        std::ifstream in(dataPath);
        if (!in.is_open())
            throw std::runtime_error("InitializationFromOdomRange: cannot open file: " + dataPath);

        std::string line;
        while (std::getline(in, line)) {
            if (line.empty() || line[0] == '#') continue;
            std::istringstream ss(line);
            std::string tag; ss >> tag;

            if (tag == "EDGE_SE2") {
                // Format: EDGE_SE2 <ts> <ida> <idb> dx dy dtheta I11 I12 I13 I22 I23 I33
                double ts, dx, dy, dth, I11, I12, I13, I22, I23, I33;
                std::string ida, idb;
                if (!(ss >> ts >> ida >> idb >> dx >> dy >> dth >> I11 >> I12 >> I13 >> I22 >> I23 >> I33))
                    continue;

                char ca = ida[0], cb = idb[0];
                size_t ia = std::stoull(ida.substr(1)), ib = std::stoull(idb.substr(1));
                robots.insert(ca); robots.insert(cb);
                robotPoseIds[ca].insert(ia);
                robotPoseIds[cb].insert(ib);

                rawOdom.push_back(RawOdom{ca, cb, ia, ib, dx, dy, dth});
            }
            else if (tag == "EDGE_RANGE") {
                // Format: EDGE_RANGE <ts> <ida> <idl> range sigma
                double ts, range, sigma;
                std::string ida, idl;
                if (!(ss >> ts >> ida >> idl >> range >> sigma)) continue;

                char ca = ida[0];
                size_t ia = std::stoull(ida.substr(1));
                robots.insert(ca);
                robotPoseIds[ca].insert(ia);

                ++numRanges; // spheres are 1:1 with range measurements
            }
            else if (tag == "VERTEX_XY") {
                // Format: VERTEX_XY L<j> x y
                std::string lid; double x, y;
                if (!(ss >> lid >> x >> y)) continue;
                if (lid.empty() || lid[0] != 'L') continue;
                size_t jl = std::stoull(lid.substr(1));
                vertexXY_GT[jl] = Vector2d(x, y);
            }
            else if (tag == "VERTEX_XYZ") {
                // Format: VERTEX_XYZ L<j> x y z
                std::string lid; double x, y, z;
                if (!(ss >> lid >> x >> y >> z)) continue;
                if (lid.empty() || lid[0] != 'L') continue;
                size_t jl = std::stoull(lid.substr(1));
                vertexXYZ_GT[jl] = Vector3d(x, y, z);
            }
            // VERTEX_SE2 / VERTEX_SE3:QUAT lines can be ignored for initialization
        }
        in.close();

        // --------- Reindex robot poses to global contiguous ids ---------
        std::map<std::pair<char,size_t>, size_t> robotLocal2Global;
        size_t nextGlobal = 0;
        for (char r : robots) {
            for (size_t localIdx : robotPoseIds[r]) {
                robotLocal2Global[{r, localIdx}] = nextGlobal++;
            }
        }

        // --------- Convert odometry edges to global ids, keep only consecutive ---------
        struct OdomEdge { size_t i, j; MatrixXd R; VectorXd t; };
        std::vector<OdomEdge> odomEdges;

        auto rot2 = [](double th) {
            MatrixXd R(2,2);
            double c = std::cos(th), s = std::sin(th);
            R << c, -s, s, c;
            return R;
        };

        for (const auto& e : rawOdom) {
            size_t ig = robotLocal2Global.at({e.ca, e.ia});
            size_t jg = robotLocal2Global.at({e.cb, e.ib});
            if (jg != ig + 1) continue; // odometry = consecutive only

            OdomEdge oe;
            oe.i = ig; oe.j = jg;
            oe.R = rot2(e.dtheta);
            oe.t = (Vector2d() << e.dx, e.dy).finished();
            odomEdges.push_back(std::move(oe));
        }

        std::sort(odomEdges.begin(), odomEdges.end(),
                  [](const OdomEdge& a, const OdomEdge& b){ return a.i < b.i; });

        // --------- Chain odometry ---------
        struct PoseW { MatrixXd R; VectorXd t; };
        std::map<size_t, PoseW> poseMap;

        if (!odomEdges.empty()) {
            const size_t startId = odomEdges.front().i;
            poseMap[startId] = PoseW{ MatrixXd::Identity(d, d), VectorXd::Zero(d) };

            for (const auto& e : odomEdges) {
                if (!poseMap.count(e.i))
                    poseMap[e.i] = PoseW{ MatrixXd::Identity(d, d), VectorXd::Zero(d) };

                const PoseW& Pi = poseMap[e.i];
                PoseW Pj;
                Pj.R = Pi.R * e.R;
                Pj.t = Pi.t + Pi.R * e.t;
                poseMap[e.j] = std::move(Pj);
            }
        }

        // Ensure every pose key exists (identity if it wasn’t in a chained edge)
        for (const auto& kv : robotLocal2Global) {
            size_t g = kv.second;
            if (!poseMap.count(g)) {
                poseMap[g] = PoseW{ MatrixXd::Identity(d, d), VectorXd::Zero(d) };
            }
        }

        // --------- Build Values: poses ---------
        gtsam::Values initial;

        for (const auto& kv : poseMap) {
            const size_t key = kv.first;
            const PoseW& Pw  = kv.second;

            MatrixXd Ymat = MatrixXd::Zero(static_cast<Index>(P), d);
            Ymat.topRows(d) = Pw.R;                  // Y = [R; 0]
            auto Y = StiefelManifoldKP::projectToManifold(Ymat);

            VectorXd tP = VectorXd::Zero(static_cast<Index>(P));
            tP.head(d) = Pw.t;

            initial.insert(key, LiftedPoseDP(Y, tP));
        }

        // --------- Landmarks from GT: VERTEX_XY and VERTEX_XYZ (contiguous L0,L1,...) ---------
        // We combine both sets and assign L-keys in ascending order of the raw id (jl).
        std::vector<std::pair<size_t, VectorXd>> combinedGT;  // (raw id, first-d coords)

        for (const auto& kv : vertexXY_GT) {
            VectorXd v(d); v.setZero();
            v.head<2>() = kv.second;
            combinedGT.emplace_back(kv.first, v);
        }
        for (const auto& kv : vertexXYZ_GT) {
            VectorXd v(d); v.setZero();
            const int use = std::min<int>(d, 3);
            v.head(use) = kv.second.head(use);
            combinedGT.emplace_back(kv.first, v);
        }

        std::sort(combinedGT.begin(), combinedGT.end(),
                  [](const auto& a, const auto& b){ return a.first < b.first; });

        size_t lidx = 0;
        for (const auto& kv : combinedGT) {
            const VectorXd& dvec = kv.second;            // size d, filled from XY/XYZ
            VectorXd lP = VectorXd::Zero(static_cast<Index>(P));
            lP.head(d) = dvec;                            // put GT coords in first d entries
            initial.insert(gtsam::Symbol('L', lidx++), lP);
        }

        // --------- Range spheres: random unit vectors (one per EDGE_RANGE) ---------
        std::default_random_engine rng(std::default_random_engine::default_seed);
        std::normal_distribution<double> gauss(0.0, 1.0);

        for (size_t r = 0; r < numRanges; ++r) {
            MatrixXd X0(static_cast<Index>(P), 1);
            for (Index i = 0; i < static_cast<Index>(P); ++i) X0(i,0) = gauss(rng);
            X0.col(0).normalize();
            UnitSphereD Y(X0);
            initial.insert(gtsam::Symbol('R', r), Y);
        }

        return initial;
    }




    /**
     * @brief Initialize lifted robot poses and landmarks from a range-aided 2D odometry dataset
     *        (non-lifted range factor version, no unit-sphere variables).
     *
     * This method parses a `.g2o`-style dataset containing:
     *  - `EDGE_SE2` (odometry edges)
     *  - `EDGE_RANGE` (range measurements, only used to ensure pose IDs are included)
     *  - `VERTEX_XY` (ground-truth landmark positions)
     *
     * It reconstructs absolute robot poses by chaining consecutive odometry edges, remaps
     * robot-specific local pose IDs to contiguous global indices, and initializes:
     *  - `LiftedPoseDP` variables for all robot poses (rotation and translation).
     *  - Landmark position vectors (initialized from dataset coordinates).
     *
     * @param dataPath  Path to the `.g2o` dataset file.
     * @param P         Ambient dimension of the lifted representation (must satisfy P ≥ d for 2D).
     * @return          `gtsam::Values` containing initialized poses and landmarks.
     *
     * @throws std::invalid_argument if P < d.
     * @throws std::runtime_error if the dataset file cannot be opened.
     *
     * @note Unlike the lifted range-aided initialization, this version does not create
     *       unit-sphere variables for range factors; it is intended for use with normal range factors.
     * @note Pose IDs from each robot are remapped to contiguous global IDs in ascending order.
     * @note Landmarks are assigned contiguous IDs starting from 'L0' in order of their raw file IDs.
     */
    gtsam::Values InitializationFromOdom_NonLifted(const std::string& dataPath, size_t P) {
        using namespace Eigen;

        if (P < static_cast<size_t>(d))
            throw std::invalid_argument("InitializationFromOdomRange: P must satisfy P >= d.");

        // --------- Parse the file (EDGE_SE2, EDGE_RANGE, VERTEX_XY, VERTEX_XYZ) ---------
        struct RawOdom { char ca, cb; size_t ia, ib; double dx, dy, dtheta; };
        std::vector<RawOdom> rawOdom;
        size_t numRanges = 0;

        // For reindexing robot poses to global contiguous ids:
        std::set<char> robots;
        std::map<char, std::set<size_t>> robotPoseIds;

        // Ground-truth landmarks
        std::map<size_t, Vector2d> vertexXY_GT;   // raw local lmk id -> (x,y)
        std::map<size_t, Vector3d> vertexXYZ_GT;  // raw local lmk id -> (x,y,z)

        std::ifstream in(dataPath);
        if (!in.is_open())
            throw std::runtime_error("InitializationFromOdomRange: cannot open file: " + dataPath);

        std::string line;
        while (std::getline(in, line)) {
            if (line.empty() || line[0] == '#') continue;
            std::istringstream ss(line);
            std::string tag; ss >> tag;

            if (tag == "EDGE_SE2") {
                // Format: EDGE_SE2 <ts> <ida> <idb> dx dy dtheta I11 I12 I13 I22 I23 I33
                double ts, dx, dy, dth, I11, I12, I13, I22, I23, I33;
                std::string ida, idb;
                if (!(ss >> ts >> ida >> idb >> dx >> dy >> dth >> I11 >> I12 >> I13 >> I22 >> I23 >> I33))
                    continue;

                char ca = ida[0], cb = idb[0];
                size_t ia = std::stoull(ida.substr(1)), ib = std::stoull(idb.substr(1));
                robots.insert(ca); robots.insert(cb);
                robotPoseIds[ca].insert(ia);
                robotPoseIds[cb].insert(ib);

                rawOdom.push_back(RawOdom{ca, cb, ia, ib, dx, dy, dth});
            }
            else if (tag == "EDGE_RANGE") {
                // Format: EDGE_RANGE <ts> <ida> <idl> range sigma
                double ts, range, sigma;
                std::string ida, idl;
                if (!(ss >> ts >> ida >> idl >> range >> sigma)) continue;

                char ca = ida[0];
                size_t ia = std::stoull(ida.substr(1));
                robots.insert(ca);
                robotPoseIds[ca].insert(ia);

                ++numRanges; // spheres are 1:1 with range measurements
            }
            else if (tag == "VERTEX_XY") {
                // Format: VERTEX_XY L<j> x y
                std::string lid; double x, y;
                if (!(ss >> lid >> x >> y)) continue;
                if (lid.empty() || lid[0] != 'L') continue;
                size_t jl = std::stoull(lid.substr(1));
                vertexXY_GT[jl] = Vector2d(x, y);
            }
            else if (tag == "VERTEX_XYZ") {
                // Format: VERTEX_XYZ L<j> x y z
                std::string lid; double x, y, z;
                if (!(ss >> lid >> x >> y >> z)) continue;
                if (lid.empty() || lid[0] != 'L') continue;
                size_t jl = std::stoull(lid.substr(1));
                vertexXYZ_GT[jl] = Vector3d(x, y, z);
            }
            // VERTEX_SE2 / VERTEX_SE3:QUAT lines can be ignored for initialization
        }
        in.close();

        // --------- Reindex robot poses to global contiguous ids ---------
        std::map<std::pair<char,size_t>, size_t> robotLocal2Global;
        size_t nextGlobal = 0;
        for (char r : robots) {
            for (size_t localIdx : robotPoseIds[r]) {
                robotLocal2Global[{r, localIdx}] = nextGlobal++;
            }
        }

        // --------- Convert odometry edges to global ids, keep only consecutive ---------
        struct OdomEdge { size_t i, j; MatrixXd R; VectorXd t; };
        std::vector<OdomEdge> odomEdges;

        auto rot2 = [](double th) {
            MatrixXd R(2,2);
            double c = std::cos(th), s = std::sin(th);
            R << c, -s, s, c;
            return R;
        };

        for (const auto& e : rawOdom) {
            size_t ig = robotLocal2Global.at({e.ca, e.ia});
            size_t jg = robotLocal2Global.at({e.cb, e.ib});
            if (jg != ig + 1) continue; // odometry = consecutive only

            OdomEdge oe;
            oe.i = ig; oe.j = jg;
            oe.R = rot2(e.dtheta);
            oe.t = (Vector2d() << e.dx, e.dy).finished();
            odomEdges.push_back(std::move(oe));
        }

        std::sort(odomEdges.begin(), odomEdges.end(),
                  [](const OdomEdge& a, const OdomEdge& b){ return a.i < b.i; });

        // --------- Chain odometry ---------
        struct PoseW { MatrixXd R; VectorXd t; };
        std::map<size_t, PoseW> poseMap;

        if (!odomEdges.empty()) {
            const size_t startId = odomEdges.front().i;
            poseMap[startId] = PoseW{ MatrixXd::Identity(d, d), VectorXd::Zero(d) };

            for (const auto& e : odomEdges) {
                if (!poseMap.count(e.i))
                    poseMap[e.i] = PoseW{ MatrixXd::Identity(d, d), VectorXd::Zero(d) };

                const PoseW& Pi = poseMap[e.i];
                PoseW Pj;
                Pj.R = Pi.R * e.R;
                Pj.t = Pi.t + Pi.R * e.t;
                poseMap[e.j] = std::move(Pj);
            }
        }

        // Ensure every pose key exists (identity if it wasn’t in a chained edge)
        for (const auto& kv : robotLocal2Global) {
            size_t g = kv.second;
            if (!poseMap.count(g)) {
                poseMap[g] = PoseW{ MatrixXd::Identity(d, d), VectorXd::Zero(d) };
            }
        }

        // --------- Build Values: poses ---------
        gtsam::Values initial;

        for (const auto& kv : poseMap) {
            const size_t key = kv.first;
            const PoseW& Pw  = kv.second;

            MatrixXd Ymat = MatrixXd::Zero(static_cast<Index>(P), d);
            Ymat.topRows(d) = Pw.R;                  // Y = [R; 0]
            auto Y = StiefelManifoldKP::projectToManifold(Ymat);

            VectorXd tP = VectorXd::Zero(static_cast<Index>(P));
            tP.head(d) = Pw.t;

            initial.insert(key, LiftedPoseDP(Y, tP));
        }

        // --------- Landmarks from GT: VERTEX_XY and VERTEX_XYZ (contiguous L0,L1,...) ---------
        // We combine both sets and assign L-keys in ascending order of the raw id (jl).
        std::vector<std::pair<size_t, VectorXd>> combinedGT;  // (raw id, first-d coords)

        for (const auto& kv : vertexXY_GT) {
            VectorXd v(d); v.setZero();
            v.head<2>() = kv.second;
            combinedGT.emplace_back(kv.first, v);
        }
        for (const auto& kv : vertexXYZ_GT) {
            VectorXd v(d); v.setZero();
            const int use = std::min<int>(d, 3);
            v.head(use) = kv.second.head(use);
            combinedGT.emplace_back(kv.first, v);
        }

        std::sort(combinedGT.begin(), combinedGT.end(),
                  [](const auto& a, const auto& b){ return a.first < b.first; });

        size_t lidx = 0;
        for (const auto& kv : combinedGT) {
            const VectorXd& dvec = kv.second;            // size d, filled from XY/XYZ
            VectorXd lP = VectorXd::Zero(static_cast<Index>(P));
            lP.head(d) = dvec;                            // put GT coords in first d entries
            initial.insert(gtsam::Symbol('L', lidx++), lP);
        }
//        initial.print("Test initials construction:");
        return initial;
    }

    /**
     * @brief Build a non-lifted factor graph for range-aided pose graph optimization.
     *
     * Constructs a `NonlinearFactorGraph` containing:
     *  - `NormalRangeFactor{2,3}` factors for range measurements.
     *  - `SEsyncFactor{2,3}` factors for odometry/pose measurements.
     *
     * The dimensionality `d` determines whether 2D or 3D factors are created.
     * Noise models are built from the measurement parameters (`sigma`, `kappa`, `tau`)
     * and the ambient dimension `p`.
     *
     * @param p Ambient dimension of the lifted representation.
     * @return  A `NonlinearFactorGraph` containing range and odometry factors.
     *
     * @note This version is for the non-lifted range factor formulation (no unit-sphere variables).
     */
    NonlinearFactorGraph buildGraphAtLevel_NonLifted(size_t p) {
        NonlinearFactorGraph inputGraph;
        size_t k = 0;
        for (const auto meas : measurements_.rangeMeasurements)
        {
            Vector sigmas = Vector::Zero(1);
            sigmas.head(1).setConstant(sqrt(1 * meas.sigma));
            auto noise = noiseModel::Diagonal::Sigmas(sigmas);
            if (d == 2) {
                inputGraph.emplace_shared<NormalRangeFactor2>( meas.i, Symbol('L', meas.j), meas.range, p , noise);
            }
            else if (d == 3)
            {
                inputGraph.emplace_shared<NormalRangeFactor3>( meas.i, Symbol('L', meas.j), meas.range, p , noise);
            }
            k++;
        } // range factor
        for (const auto meas : measurements_.poseMeasurements)
        {
            Vector sigmas = Vector::Zero(p * d + p);
            sigmas.head(p * d).setConstant(std::sqrt(1.0 / (1 * meas.kappa)));
            sigmas.tail(p).setConstant(std::sqrt(1.0 / (1 * meas.tau)));
            auto noise = noiseModel::Diagonal::Sigmas(sigmas);

            if (d == 2) {
                double theta = std::atan2(meas.R(1, 0), meas.R(0, 0)); // Extract angle
                auto odom = gtsam::Pose2(meas.t(0), meas.t(1), theta);
                inputGraph.emplace_shared<SEsyncFactor2>(meas.i,  meas.j, meas.R, meas.t, p, noise);
            }
            else if (d ==3 ) {
                gtsam::Rot3 rot(meas.R);
                gtsam::Point3 trans(meas.t(0), meas.t(1), meas.t(2));
                auto odom = gtsam::Pose3(rot, trans);
                inputGraph.emplace_shared<SEsyncFactor3>(meas.i, meas.j, meas.R, meas.t, p, noise);
            }
            else {
                std::cerr << "Un" << std::endl;
            }
        } // PGO Factor
        return inputGraph;
    }

    std::optional<CertificateResults> LocalSearchWithOdomInitials(size_t p)
    {

        Values Qstar;
        auto t6 = CFGStopwatch::tick();
        std::cout << "Starting Local search-based optimization "<< std::endl;
        auto lmParams = opts_.lmParams;
        lmParams.maxIterations    = opts_.maxIterations;
        lmParams.relativeErrorTol = opts_.relativeErrorTol;
        lmParams.absoluteErrorTol = opts_.absoluteErrorTol;
        lmParams.verbosityLM      = opts_.verbosityLM;
        auto lm = std::make_shared<LevenbergMarquardtOptimizer>(currentGraph_, currentValues_, lmParams);
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


// Convenience aliases
using CertifiableRangeSLAM2 = CertifiableRangeSLAM<2>;
using CertifiableRangeSLAM3 = CertifiableRangeSLAM<3>;


#endif //STIEFELMANIFOLDEXAMPLE_CERTIFIABLERANGESLAM_H