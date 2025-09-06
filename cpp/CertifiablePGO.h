//
// Created by jason on 4/30/25.
//

#ifndef STIEFELMANIFOLDEXAMPLE_CERTIFIABLEPGO_H
#define STIEFELMANIFOLDEXAMPLE_CERTIFIABLEPGO_H
#pragma once
#include "Certifiable_problem.h"

/**
 * @brief Certifiable Pose Graph Optimization (PGO)
 *
 * Implements a hierarchy of relaxations of PGO via rank‐d Stiefel manifolds,
 * certifies solutions using eigenvalue tests, and refines by lifting along
 * descent directions when necessary.
 *
 * @tparam d  Ambient pose dimension (must be 2 or 3).
 */
template <size_t d>
class CertifiablePGO : public CertifiableProblem {
    static_assert(d == 2 || d == 3, "CertifiablePGO only supports d = 2 or 3.");
public:
    /**
     * @brief Construct with initial rank and measurement data.
     * @param p             Initial relaxation rank.
     * @param measurements  Parsed measurement struct (num_poses, etc.).
     */
    CertifiablePGO(size_t p, const DataParser::Measurement& measurements)
        : CertifiableProblem(d, p, measurements)
    {
        certificateResults_.startingRank = p;
    }


    /**
     * @brief Initialize graph, data matrix, and random values; record init time.
     */
    void init() {
        auto t0 = CFGStopwatch::tick();
        if (opts_.initType == CertifiableProblemOpts::InitType::Random){
            currentValues_  = randomInitAtLevelP(currentRank_);
        }
        else if (opts_.initType == CertifiableProblemOpts::InitType::Odom){
            currentValues_  = InitializationFromOdom(measurements_, currentRank_);
        }
        else {
            throw std::runtime_error("Unknown initialization type in CertifiableProblemOpts::initType");
        }
        currentGraph_   = buildGraphAtLevel(currentRank_);
        M_              = recoverDataMatrixFromHessian();
        auto t1 = CFGStopwatch::tock(t0);
        certificateResults_.initialization_time.push_back(t1);
    }

    // Only needed for incremental, only implemented to make class non abstacrt
    NonlinearFactorGraph buildGraphAtLevelUsingGraph(size_t p) override {
        return NonlinearFactorGraph();
    }

    // Only needed for incremental, only implemented to make class non abstacrt
    Values buildInitialValuesAtLevel(size_t p) {
        return randomInitAtLevelP(p);
    }
    Values poseInitAtLevelP(Values poseInital,size_t Pmin)
    {
        Values lifted;
        using Pose = std::conditional_t<d==2, Pose2, Pose3>;
        for (auto& k : poseInital.keys()) {
            const Pose p = poseInital.template at<Pose>(k);
            lifted.insert(
              k,
              LiftedPoseDP(
                StiefelManifoldKP::Lift(Pmin, p.rotation().matrix()),
                p.translation().matrix()
              )
            );
        }
        return lifted;
    }

    Values ExportPose3(Matrix& R)
    {
        Values finalposes;
        // Insert SE3 poses
        for (auto key : currentValues_.keys()) {
            Symbol s(key);
            if (s.chr() != 'x') continue;
            size_t i  = s.index();
            finalposes.insert(i, Pose3( Rot3(R.block(0, num_pose_ + i * d, d, d)), R.block(0, i, d, 1)));
        }

        return finalposes;
    }

    Values ExportPose2(Matrix& R)
    {
        Values finalposes;
        for (auto key : currentValues_.keys()) {
            // Only handle X(i) keys (skip any landmark or other variables)
            gtsam::Symbol s(key);
            if (s.chr() != 'x') continue;

            size_t i = s.index();  // this is the “i” used when you built R

            // rotation block: rows [0..d), cols [num_pose_ + i*d .. num_pose_ + (i+1)*d)
            Rot2 rot = Rot2::fromCosSin(
                    R.block(0, num_pose_ + i * d, d, d)(0, 0),
                    R.block(0, num_pose_ + i * d, d, d)(1, 0)
            );
            finalposes.insert(key, Pose2(rot, R.block(0, i, d, 1)));
        }

        return finalposes;
    }

    /**
     * @brief Assemble the factor graph at relaxation level p.
     * @param p  Relaxation rank.
     * @return   Factor graph containing only pose‐to‐pose factors.
     */
    NonlinearFactorGraph buildGraphAtLevel(size_t p) override {
        NonlinearFactorGraph inputGraph;
        for (const auto& meas : measurements_.poseMeasurements) {
            // Construct diagonal noise sigmas from kappa and tau
            Vector sigmas = Vector::Zero(p * d + p);
            sigmas.head(p * d).setConstant(std::sqrt(1.0 / (2 * meas.kappa)));
            sigmas.tail(p).setConstant(std::sqrt(1.0 / (2 * meas.tau)));
            auto noise = noiseModel::Diagonal::Sigmas(sigmas);

            // Emplace the appropriate SE‐sync factor
            if constexpr (d == 2) {
                inputGraph.emplace_shared<SEsyncFactor2>(
                    meas.i, meas.j, meas.R, meas.t, p, noise
                );
            } else {
                inputGraph.emplace_shared<SEsyncFactor3>(
                    meas.i, meas.j, meas.R, meas.t, p, noise
                );
            }
        }
        return inputGraph;
    }

    /**
     * @brief Compute the dense block‐diagonal certificate matrix.
     * @param Y  Variable matrix.
     * @return   Dense block matrix of size (d × d·num_poses).
     */
    Matrix computeLambdaBlocks(const Matrix& Y) override {
        Matrix SY = M_ * Y;
        Matrix Yt = Y.transpose();
        Matrix LambdaBlocks(d, num_pose_ * d);
        size_t offset = num_pose_;

        for (size_t i = 0; i < num_pose_; ++i) {
            Matrix P = SY.block(offset + i * d, 0, d, Y.cols())
                     * Yt.block(0, offset + i * d, Y.cols(), d);
            LambdaBlocks.block(0, i * d, d, d) = 0.5 * (P + P.transpose());
        }
        return LambdaBlocks;
    }

    /**
     * @brief Convert dense Λ_blocks into a sparse Λ matrix for certification.
     * @param LambdaBlocks  Dense block matrix.
     * @return              Sparse certificate matrix Λ.
     */
    SparseMatrix computeLambdaFromLambdaBlocks(
        const Matrix& LambdaBlocks) override
    {
        std::vector<Eigen::Triplet<Scalar>> elements;
        elements.reserve(d * d * num_pose_);
        size_t offset = num_pose_;
        for (size_t i = 0; i < num_pose_; ++i) {
            for (size_t r = 0; r < d; ++r) {
                for (size_t c = 0; c < d; ++c) {
                    elements.emplace_back(
                        offset + i * d + r,
                        offset + i * d + c,
                        LambdaBlocks(r, i * d + c)
                    );
                }
            }
        }
        SparseMatrix Lambda(offset + d * num_pose_, offset + d * num_pose_);
        Lambda.setFromTriplets(elements.begin(), elements.end());
        return Lambda;
    }

    /**
     * @brief Build the element matrix S from current Values.
     * @param values  GTSAM Values containing LiftedPoseDP variables.
     * @return        Sparse matrix S of size.
     */
    SparseMatrix elementMatrix(const Values& values) override {
        const size_t N = num_pose_, p = currentRank_;
        std::vector<Eigen::Triplet<double>> triplets;
        triplets.reserve(N * p * (d + 1));

        for (const auto& kv : values.extract<LiftedPoseDP>()) {
//            size_t i = kv.first;
            gtsam::Key rawKey   = kv.first;        // full 64-bit key
            gtsam::Symbol sym(rawKey);            // decode the key
            size_t i = sym.index();

            const auto& pose = kv.second;
            const auto& t   = pose.get_TranslationVector(); // p × 1
            const auto& mat = pose.matrix();                // p × d

            // Translation entries
            for (size_t row = 0; row < p; ++row)
                triplets.emplace_back(i, row, t(row));

            // Rotation (Stiefel) entries
            size_t row0 = N + i * d;
            for (size_t row = 0; row < p; ++row)
                for (size_t col = 0; col < d; ++col)
                    triplets.emplace_back(row0 + col, row, mat(row, col));
        }

        SparseMatrix S(N + N * d, p);
        S.setFromTriplets(triplets.begin(), triplets.end());
        return S;
    }

    /**
     * @brief Randomly initialize all poses at level Pmin using uniform Stiefel and random translation.
     * @param Pmin  Target relaxation rank.
     * @return      GTSAM Values container with random LiftedPoseDP entries.
     */
    Values randomInitAtLevelP(const size_t Pmin) override {
        Values initial;
        for (size_t j = 0; j < num_pose_; ++j) {
            StiefelManifoldKP Y =
                StiefelManifoldKP::Random(std::default_random_engine::default_seed, d, Pmin);
            Vector trans = Vector::Random(Pmin);
            initial.insert(j, LiftedPoseDP(Y, trans));
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
        Matrix Ydot = Matrix::Zero(v.size(), p);
        Ydot.rightCols<1>() = v;

        for (const auto& kv : values.extract<LiftedPoseDP>()) {
            size_t idx = gtsam::symbolIndex(kv.first);
            const auto& Y = kv.second.get_Y();

            // Extract ambient gradient blocks
            Matrix tangM = Ydot.block(num_pose_ + idx * d, 0, d, p).transpose();
            Vector transV = Ydot.block(idx, 0, 1, p).transpose();

            Vector xi = StiefelManifoldKP::Vectorize(tangM);
            Vector tVec = Y.G_.transpose() * xi;

            Vector combined(tVec.size() + transV.size());
            combined << tVec, transV;
            delta.insert(idx, combined);
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
        Matrix result = Ydot;
        size_t offset_t  = num_pose_;
        size_t rot_mat_sz = d * num_pose_;
        result.block(offset_t, 0, rot_mat_sz, p) =
            StiefelManifoldKP::Proj(
                Y.block(offset_t, 0, rot_mat_sz, p).transpose(),
                result.block(offset_t, 0, rot_mat_sz, p).transpose()
            ).transpose();
        return result;
    }

    /**
     * @brief Recover the data matrix from the current factor graph.
     * @return Sparse data matrix L of size.
     */
    SparseMatrix recoverDataMatrixFromGraph() override {
        std::vector<Eigen::Triplet<Scalar>> triplets;
        triplets.reserve(currentGraph_.size() * (4 + 2*d*d + 2*d + 2));

        if (d == 2) {
            for (auto& f_ptr : currentGraph_) {
                if (auto factor = std::dynamic_pointer_cast<SEsyncFactor2>(f_ptr))
                    factor->appendBlocksFromFactor(num_pose_, triplets);
            }
        } else {
            for (auto& f_ptr : currentGraph_) {
                if (auto factor = std::dynamic_pointer_cast<SEsyncFactor3>(f_ptr))
                    factor->appendBlocksFromFactor(num_pose_, triplets);
            }
        }

        SparseMatrix L((d + 1) * num_pose_, (d + 1) * num_pose_);
        L.setFromTriplets(triplets.begin(), triplets.end());
        return L;
    }

    void recoverHessian(std::map<std::pair<Key,Key>,Matrix> &HMap)
    {
        if (d == 2) {
            for (auto& f_ptr : currentGraph_) {
                if (auto factor = std::dynamic_pointer_cast<SEsyncFactor2>(f_ptr))
                {
                    factor->computeHessian(HMap);
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
            }
        }
    }

    SparseMatrix recoverDataMatrixFromHessian() {

        std::map<std::pair<Key,Key>,Matrix> HMap;
        recoverHessian(HMap);

        Ordering order = Ordering::Natural(currentGraph_);
        const size_t N = order.size();
        std::vector<size_t> dims(N), offsets(N);

        std::map<Key,size_t> Key2index;
        for (size_t i = 0; i < N; ++i) {
            const Key k = order[i];
            // Use the same getter you use later for your Stiefel blocks:
            dims[i]    = currentValues_.at<LiftedPoseDP>(k).get_Cols();
            offsets[i] = (i == 0 ? 0u : offsets[i-1] + dims[i-1]);
            Key2index[k] = i;
        }
        const size_t totalDim = offsets.back() + dims.back();

        // 2) Gather triplets
        std::vector<Eigen::Triplet<Scalar>> triplets;
        triplets.reserve(HMap.size() * (4*d)*(4*d));
        appendPGOHessian(triplets,HMap, offsets, order, Key2index);

        SparseMatrix L(totalDim, totalDim);
        L.setFromTriplets(triplets.begin(), triplets.end());
        // 1) make sure the sparse structure is finalized
        L.makeCompressed();
        return L;
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
//            auto t0 = std::chrono::high_resolution_clock::now();
            Qstar = tryOptimizingAtLevel(p);
            setCurrentValues(Qstar);
//            auto t1 = std::chrono::high_resolution_clock::now();
//            certificateResults_.elapsed_optimization_times.push_back((std::chrono::duration<double, std::milli> (t1 - t0)).count());

//            auto t2 = std::chrono::high_resolution_clock::now();
//            NonlinearFactorGraph nonlinear_graph = buildGraphAtLevel(p);
//
//            auto linear_graph = nonlinear_graph.linearize(Qstar);
//            auto grad_norm = linear_graph->gradientAtZero();
//            std::cout << "Gradient norm at level p = " << p << " is : " << grad_norm.norm() << std::endl;
//            certificateResults_.gradnorm.push_back(grad_norm.norm());
//            auto t3 = std::chrono::high_resolution_clock::now();
//            certificateResults_.initialization_time.push_back((std::chrono::duration<double, std::milli> (t3 - t2)).count());

            auto t4 = CFGStopwatch::tick();
            SparseMatrix S = elementMatrix(Qstar);
            Matrix lambdaBlocks = computeLambdaBlocks(S);
            SparseMatrix Lambda = computeLambdaFromLambdaBlocks(lambdaBlocks);

            SparseMatrix M = getDataMatrix();
            Scalar obj = evaluateObjective(S, M);

            // For test
//            std::cout << "GTSAM's error : " << nonlinear_graph.error(Qstar) <<std::endl;
//            std::cout << "SDP's error : " << evaluateObjective(S, M)  <<std::endl;

            bool success = false;
            Scalar eta;
            if (opts_.useAbsoluteEta == true) {
                eta = opts_.eta;
            } else {
                eta = thresholdVal(obj * opts_.REL_CERT_ETA, opts_.MIN_CERT_ETA, opts_.MAX_CERT_ETA);
            }
//            Scalar eta = opts_.eta;
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
                increaseCurrentRank();
                currentValues_ = initializeWithDescentDirection(Qstar, M, v, theta, 1e-2);
            } else {
                std::cout << "Solution verified at level p = " << p << std::endl;
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
#pragma omp parallel for
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
    void ExportData(const string &path, const Eigen::MatrixXd &R, bool g2o) override {
        if (g2o) {
            Values finalposes;
            if (d == 3) {
                // Insert SE3 poses
                for (auto i = 0; i < num_pose_; ++i) {
                    finalposes.insert(
                        i,
                        Pose3(
                            Rot3(R.block(0, num_pose_ + i * d, d, d)),
                            R.block(0, i, d, 1)
                        )
                    );
                }
            } else {
                // Insert SE2 poses
                for (auto i = 0; i < num_pose_; ++i) {
                    Rot2 rot = Rot2::fromCosSin(
                        R.block(0, num_pose_ + i * d, d, d)(0, 0),
                        R.block(0, num_pose_ + i * d, d, d)(1, 0)
                    );
                    finalposes.insert(i, Pose2(rot, R.block(0, i, d, 1)));
                }
            }
            writeG2o(currentGraph_, finalposes, path + ".g2o");
        } else {
            std::ofstream file(path + ".txt");
            if (d == 2) {
                for (auto i = 0; i < num_pose_; ++i) {
                    Eigen::Matrix3d R3 = Eigen::Matrix3d::Identity();
                    R3.topLeftCorner<2,2>() =
                        R.block(0, num_pose_ + i*2, 2, 2);
                    auto q = Eigen::Quaterniond(R3);
                    Vector t = R.block(0, i, d, 1);
                    file << i << " " << t(0) << " " << t(1) << " " << 0
                         << " " << q.x() << " " << q.y()
                         << " " << q.z() << " " << q.w() << "\n";
                }
            } else {
                for (auto i = 0; i < num_pose_; ++i) {
                    Quaternion q(R.block<3,3>(0, num_pose_ + i*3));
                    Vector t = R.block(0, i, d, 1);
                    file << i << " " << t(0) << " " << t(1) << " " << 0
                         << " " << t(2) << " " << q.x() << " "
                         << q.y() << " " << q.z() << " "
                         << q.w() << "\n";
                }
            }
            file.close();
        }
    }



//    struct OdomEdge {
//        size_t i{0}, j{0};   // consecutive pose ids with j = i+1
//        Eigen::MatrixXd R;          // d x d
//        Eigen::VectorXd t;          // d
//    };
//
//    Values InitializationFromOdom(std::string g2oPath, size_t P) {
//        std::ifstream in(g2oPath);
//        if (!in.is_open()) {
//            throw std::runtime_error("InitializationFromOdom: cannot open file: " + g2oPath);
//        }
//
//        // Parse only odometry edges (j == i+1). Supports 2D (SE2) and 3D (SE3:QUAT).
//        std::vector<OdomEdge> edges;
////        int d = -1;  // rotation/translation dimension (2 or 3)
//
//        std::string line;
//        while (std::getline(in, line)) {
//            if (line.empty() || line[0] == '#') continue;
//            std::istringstream iss(line);
//            std::string tag;
//            iss >> tag;
//
//            if (tag == "EDGE_SE2") {
//                size_t i, j;
//                double dx, dy, dtheta;
//                if (!(iss >> i >> j >> dx >> dy >> dtheta)) continue;
//                if (j != i + 1) continue;  // skip loop-closures
//
////                d = (d == -1 ? 2 : d);
//                if (d != 2) continue;      // skip mixed-dim files
//
//                double c = std::cos(dtheta), s = std::sin(dtheta);
//                Eigen::MatrixXd Rij = Eigen::MatrixXd::Identity(2, 2);
//                Rij(0,0) =  c; Rij(0,1) = -s;
//                Rij(1,0) =  s; Rij(1,1) =  c;
//
//                Eigen::VectorXd tij(2);
//                tij << dx, dy;
//
//                edges.push_back({i, j, Rij, tij});
//            } else if (tag == "EDGE_SE3:QUAT") {
//                size_t i, j;
//                double dx, dy, dz, qx, qy, qz, qw;
//                if (!(iss >> i >> j >> dx >> dy >> dz >> qx >> qy >> qz >> qw)) continue;
//                if (j != i + 1) continue;  // skip loop-closures
//
////                d = (d == -1 ? 3 : d);
//                if (d != 3) continue;      // skip mixed-dim files
//
//                Eigen::Quaterniond q(qw, qx, qy, qz);
//                q.normalize();
//                Eigen::MatrixXd Rij = q.toRotationMatrix();
//
//                Eigen::VectorXd tij(3);
//                tij << dx, dy, dz;
//
//                edges.push_back({i, j, Rij, tij});
//            }
//            // ignore other tags (VERTEX_*, EDGE_BEARING*, etc.)
//        }
//        in.close();
//
//        if (edges.empty()) {
//            throw std::runtime_error("InitializationFromOdom: no odometry (consecutive) edges found.");
//        }
//        if (d <= 0) {
//            throw std::runtime_error("InitializationFromOdom: could not determine problem dimension.");
//        }
//        if (P < static_cast<size_t>(d)) {
//            throw std::invalid_argument("InitializationFromOdom: P must satisfy P >= d.");
//        }
//
//        // Sort edges by i to ensure forward chaining.
//        std::sort(edges.begin(), edges.end(), [](const OdomEdge& a, const OdomEdge& b) {
//            return a.i < b.i;
//        });
//
//        // Chain odometry: set the first pose to identity.
//        // World pose T_j = T_i * T_ij  => R_j = R_i * R_ij,  t_j = t_i + R_i * t_ij
//        struct Pose {
//            Eigen::MatrixXd R;  // d x d
//            Eigen::VectorXd t;  // d
//        };
//        std::map<size_t, Pose> poses;
//
//        const size_t startId = edges.front().i;
//        poses[startId] = Pose{Eigen::MatrixXd::Identity(d, d), Eigen::VectorXd::Zero(d)};
//
//        for (const auto& e : edges) {
//            auto it = poses.find(e.i);
//            if (it == poses.end()) {
//                // If a gap exists (e.g., file starts at 1->2 but we never saw 0),
//                // start a new chain with identity at e.i.
//                poses[e.i] = Pose{Eigen::MatrixXd::Identity(d, d), Eigen::VectorXd::Zero(d)};
//                it = poses.find(e.i);
//            }
//            const Pose& Pi = it->second;
//
//            Pose Pj;
//            Pj.R = Pi.R * e.R;
//            Pj.t = Pi.t + Pi.R * e.t;
//
//            poses[e.j] = std::move(Pj);
//        }
//
//        // Build lifted initialization.
//        Values initial;
//        for (const auto& kv : poses) {
//            const size_t id = kv.first;
//            const Pose& Pwd = kv.second;
//
//            // Lift rotation: embed R (d x d) into Y (P x d) as [R; 0] and (optionally) project.
//            Eigen::MatrixXd Ymat = Eigen::MatrixXd::Zero(static_cast<Eigen::Index>(P), d);
//            Ymat.topRows(d) = Pwd.R;  // since R is orthonormal, [R;0] is already on St(p,d)
//
//            // If your API prefers explicit projection, keep this; otherwise you can drop it.
//            auto Y = StiefelManifoldKP::projectToManifold(Ymat);
//
//            // Lift translation: put t in the first d coordinates, pad zeros to length P.
//            Eigen::VectorXd tP = Eigen::VectorXd::Zero(static_cast<Eigen::Index>(P));
//            tP.head(d) = Pwd.t;
//
//            initial.insert(id, LiftedPoseDP(Y, tP));
//        }
//
//        return initial;
//    }



/**
 * @brief Initialize lifted poses from odometry measurements.
 * @param M   Parsed measurement data (poses, edges, etc.).
 * @param P   Ambient dimension of the lifted representation.
 * @return    gtsam::Values containing initialized LiftedPoseDP variables.
 */
    gtsam::Values InitializationFromOdom(const gtsam::DataParser::Measurement& M, size_t P) {
        using namespace Eigen;

        // Minimal contract checks (no heavy logic on d)
        if (P < static_cast<size_t>(d))
            throw std::invalid_argument("InitializationFromOdom: P must satisfy P >= d.");

        struct OdomEdge {
            size_t i, j;   // expect j == i+1 for odometry
            MatrixXd R;    // d x d
            VectorXd t;    // d
        };

        std::vector<OdomEdge> edges;
        std::set<size_t> allPoseIds;

        // Collect consecutive edges and all pose keys
        for (const auto& m : M.poseMeasurements) {
            allPoseIds.insert(m.i);
            allPoseIds.insert(m.j);
            if (m.j == m.i + 1) {
                edges.push_back(OdomEdge{m.i, m.j, m.R, m.t});
            }
        }

        // Chain odometry: T_j = T_i * T_ij
        struct Pose { MatrixXd R; VectorXd t; };
        std::map<size_t, Pose> poses;

        if (!edges.empty()) {
            std::sort(edges.begin(), edges.end(),
                      [](const OdomEdge& a, const OdomEdge& b){ return a.i < b.i; });

            const size_t startId = edges.front().i;
            poses[startId] = Pose{ MatrixXd::Identity(d, d), VectorXd::Zero(d) };

            for (const auto& e : edges) {
                if (!poses.count(e.i)) {
                    poses[e.i] = Pose{ MatrixXd::Identity(d, d), VectorXd::Zero(d) };
                }
                const Pose& Pi = poses[e.i];
                Pose Pj;
                Pj.R = Pi.R * e.R;
                Pj.t = Pi.t + Pi.R * e.t;
                poses[e.j] = std::move(Pj);
            }
        }

        // Ensure every pose id seen in measurements has an entry (identity if not chained)
        for (size_t pid : allPoseIds) {
            if (!poses.count(pid)) {
                poses[pid] = Pose{ MatrixXd::Identity(d, d), VectorXd::Zero(d) };
            }
        }

        // Build lifted initialization
        gtsam::Values initial;
        for (const auto& kv : poses) {
            const size_t key = kv.first;
            const Pose& Pw   = kv.second;

            MatrixXd Ymat = MatrixXd::Zero(static_cast<Index>(P), d);
            Ymat.topRows(d) = Pw.R;                 // Y = [R; 0]
            auto Y = StiefelManifoldKP::projectToManifold(Ymat);

            VectorXd tP = VectorXd::Zero(static_cast<Index>(P));
            tP.head(d) = Pw.t;

            initial.insert(key, LiftedPoseDP(Y, tP));
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


// Convenience aliases
using CertifiablePGO2 = CertifiablePGO<2>;
using CertifiablePGO3 = CertifiablePGO<3>;

#endif // STIEFELMANIFOLDEXAMPLE_CERTIFIABLEPGO_H

