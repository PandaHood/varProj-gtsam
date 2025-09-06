//
// Created by nikolas on 6/30/25.
//
#include "Certifiable_problem.h"
#ifndef CERTIFIABLE_INCREMENTAL_PGO_H
#define CERTIFIABLE_INCREMENTAL_PGO_H
#include "gtsam/nonlinear/LinearContainerFactor.h"
template <size_t d>
class CertifiableIncrementalPGO : public CertifiableProblem {
    static_assert(d == 2 || d == 3, "CertifiablePGO only supports d = 2 or 3.");
    NonlinearFactorGraph betweengraph_;
public:
    /**
     * @brief Construct with initial rank and measurement data.
     * @param p             Initial relaxation rank.
     * @param measurements  Parsed measurement struct (num_poses, etc.).
     */

    CertifiableIncrementalPGO(size_t p, size_t num_poses)
    : CertifiableProblem(d, p, num_poses)
    {
        certificateResults_.startingRank = p;
    }

    CertifiableIncrementalPGO()
    : CertifiableProblem(d,d,0)
    {
        certificateResults_.startingRank = d;
    }

    /**
     * @brief Initialize graph, data matrix, and random values; record init time.
     */
    void init(LevenbergMarquardtParams params, NonlinearFactorGraph graph, Values initial) {
        // parameters
        opts_.lmParams = params;
        opts_.verbosityLM = params.verbosityLM;
        opts_.maxIterations = params.maxIterations;
        opts_.relativeErrorTol = params.relativeErrorTol;
        opts_.absoluteErrorTol = params.absoluteErrorTol;

        // initalize
        betweengraph_ = graph;
        incremental_ = true;
        currentValues_ = poseInitAtLevelP(initial, d_);
        currentGraph_ = buildLiftedGraphFromGraph(graph, d_);
        M_ = recoverDataMatrixFromGraph();
    }

    // generete values of lifted poses from pose2 or pose3 variables
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

    // Random inialization if we need
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




// build lifted factor graph from between factor graph
//    NonlinearFactorGraph buildLiftedGraphFromGraph(NonlinearFactorGraph betweenGraph, size_t p)
//    {
//        NonlinearFactorGraph out;
//        using Between = BetweenFactor<std::conditional_t<d==2, Pose2, Pose3>>;
//        using SEsync   = std::conditional_t<d==2, SEsyncFactor2, SEsyncFactor3>;
//        for (auto& factorPtr : betweenGraph) {
//            auto between = std::dynamic_pointer_cast<Between>(factorPtr);
//            if (!between) continue;
//
//            // pull out σ’s…
//            auto diag = std::dynamic_pointer_cast<noiseModel::Diagonal>(between->noiseModel());
//            auto sig6 = diag->sigmas();
//            double sigma_t = sig6[0], sigma_r = sig6[(d==2?2:3)];
//
//            Vector sigmas = Vector::Zero(p*d + p);
//            sigmas.head(p*d).setConstant(sigma_r);
//            sigmas.tail(p).setConstant(sigma_t);
//            auto newNoise = noiseModel::Diagonal::Sigmas(sigmas);
//
//            auto z = between->measured();
//            Matrix M = z.rotation().matrix();
//            Matrix V = z.translation().matrix();
//
//            out.emplace_shared<SEsync>(
//              between->key1(), between->key2(), M, V, p, newNoise
//            );
//        }
//        return out;
//    }


    gtsam::NonlinearFactorGraph buildLiftedGraphFromGraph(
            const gtsam::NonlinearFactorGraph& betweenGraph,
            size_t p)
    {
        using namespace gtsam;

        // Shorthand types
        using Between =
                BetweenFactor<std::conditional_t<d==2, Pose2, Pose3>>;
        using SEsync =
                std::conditional_t<d==2, SEsyncFactor2, SEsyncFactor3>;

        NonlinearFactorGraph out;
        for (const auto& factorPtr : betweenGraph) {

            // 1) Handle BetweenFactor<Pose>
            if (auto between = std::dynamic_pointer_cast<Between>(factorPtr)) {
                auto diag = std::dynamic_pointer_cast<noiseModel::Diagonal>(
                        between->noiseModel());
                auto sig6 = diag->sigmas();
                double sigma_t = sig6[0],
                        sigma_r = sig6[(d==2?2:3)];

                Eigen::VectorXd sigmas = Eigen::VectorXd::Zero(p*d + p);
                sigmas.head(p*d).setConstant(sigma_r);
                sigmas.tail(p).setConstant(sigma_t);
                auto newNoise = noiseModel::Diagonal::Sigmas(sigmas);

                auto z = between->measured();
                Eigen::MatrixXd M = z.rotation().matrix();
                Eigen::MatrixXd V = z.translation().matrix();

                out.emplace_shared<SEsync>(
                        between->key1(), between->key2(), M, V, p, newNoise);
                continue;
            }

//            // 2) Handle any LinearContainerFactor
            if (auto lcf = std::dynamic_pointer_cast<LinearContainerFactor>(factorPtr)) {
                // extract and print all keys
                for (auto key : lcf->keys()) {
                    std::cout << "LinearContainerFactor on key: " << key << std::endl;
                }
                // optionally, you can also copy it through unchanged:
                out.push_back(lcf);
                continue;
            }

            // 2) LinearContainerFactor → rebuild with new Values


            // 3) (optional) copy through other factors unchanged
            // out.push_back(factorPtr);
        }

        return out;
    }

    Values ExportPose3(Matrix& R)
    {
        Values finalposes;
        // Insert SE3 poses
        for (auto key : currentValues_.keys()) {
            Symbol s(key);
            if (s.chr() != 'x') continue;
            size_t i  = s.index();
            finalposes.insert(key, Pose3( Rot3(R.block(0, num_pose_ + i * d, d, d)), R.block(0, i, d, 1)));
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

    // Runs the batch optimization needed for CFS smoother
    Values runIncrementalBatch(NonlinearFactorGraph graph, Values initial, LevenbergMarquardtParams params)
    {

        init(params,graph,initial);
        // solve
        Solve(d_,d_+10);

        // move back to problem dimension
        auto R = RoundSolutionS();
        switch (d_) {
        case 2: return ExportPose2(R);
        case 3: return ExportPose3(R);
        default:
            throw std::runtime_error("Unsupported dimension");
        }
    }

    /**
     * @brief Assemble the factor graph at relaxation level p.
     * @param p  Relaxation rank.
     * @return   Factor graph containing only pose‐to‐pose factors.
     */
    NonlinearFactorGraph buildGraphAtLevel(size_t p) override {
        NonlinearFactorGraph inputGraph;
        return inputGraph;
    }

    // Using between factor to construct lifted factor graph
    NonlinearFactorGraph buildGraphAtLevelUsingGraph( size_t p) override
    {
        NonlinearFactorGraph inputGraph;
        using Pose = std::conditional_t<d==2, Pose2, Pose3>;
        using SEsync   = std::conditional_t<d==2, SEsyncFactor2, SEsyncFactor3>;
            for (auto& factorPtr : betweengraph_) {
                // 1) only handle Pose3–Between factors
                auto between =
                  std::dynamic_pointer_cast<gtsam::BetweenFactor<Pose>>(factorPtr);
                if (!between) continue;

                // 2) pull out the old diagonal noise
                auto oldNoise   = between->noiseModel();
                auto diagNoise  = std::dynamic_pointer_cast<gtsam::noiseModel::Diagonal>(oldNoise);
                if (!diagNoise)
                    throw std::runtime_error("Expected Diagonal noise in BetweenFactor");

                // 3) extract the two unique sigmas [σ_t,σ_t,σ_t, σ_r,σ_r,σ_r]
                auto diag = std::dynamic_pointer_cast<noiseModel::Diagonal>(between->noiseModel());
                auto sig6 = diag->sigmas();
                double sigma_t = sig6[0], sigma_r = sig6[(d==2?2:3)];

                // 4) broadcast into a (p*d + p)-vector
                const size_t dimY = p * d;        // lifted-rotation dims
                const size_t dimT = p;            // lifted-translation dims
                gtsam::Vector sigmas(dimY + dimT);
                sigmas.head(dimY).setConstant(sigma_r);
                sigmas.tail(dimT).setConstant(sigma_t);

                // 5) make the new noise model
                auto newNoise = gtsam::noiseModel::Diagonal::Sigmas(sigmas);

                // 6) grab measurement
                const auto& z = between->measured();
                Matrix M = z.rotation().matrix();
                Matrix V = z.translation().matrix();


                inputGraph.emplace_shared<SEsync>(between->key1(), between->key2(), M, V, p, newNoise);
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
        auto t6 = std::chrono::high_resolution_clock::now();
        for (size_t p = pMin; p <= pMax; ++p) {
            std::cout << "Starting optimization at rank = " << p << std::endl;
            auto t0 = std::chrono::high_resolution_clock::now();
            setCurrentRank(p);
            Qstar = tryOptimizingAtLevel(p);
            setCurrentValues(Qstar);
            auto t1 = std::chrono::high_resolution_clock::now();
            certificateResults_.elapsed_optimization_times.push_back((std::chrono::duration<double, std::milli> (t1 - t0)).count());

            auto t2 = std::chrono::high_resolution_clock::now();
            NonlinearFactorGraph nonlinear_graph = buildGraphAtLevelUsingGraph(p);

            auto linear_graph = nonlinear_graph.linearize(Qstar);
            auto grad_norm = linear_graph->gradientAtZero();
            std::cout << "Gradient norm at level p = " << p << " is : " << grad_norm.norm() << std::endl;
            certificateResults_.gradnorm.push_back(grad_norm.norm());
            auto t3 = std::chrono::high_resolution_clock::now();
            certificateResults_.initialization_time.push_back((std::chrono::duration<double, std::milli> (t3 - t2)).count());

            auto t4 = std::chrono::high_resolution_clock::now();
            SparseMatrix S = elementMatrix(Qstar);
            Matrix lambdaBlocks = computeLambdaBlocks(S);
            SparseMatrix Lambda = computeLambdaFromLambdaBlocks(lambdaBlocks);
            SparseMatrix M = getDataMatrix();
            Scalar obj = evaluateObjective(S, M);

            // For test
            std::cout << "GTSAM's error : " << nonlinear_graph.error(Qstar) <<std::endl;
            std::cout << "SDP's error : " << evaluateObjective(S, M)  <<std::endl;

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
            auto t5 = std::chrono::high_resolution_clock::now();
            certificateResults_.verification_times.push_back((std::chrono::duration<double, std::milli> (t5 - t4)).count());
            if (!success) {
                increaseCurrentRank();
                currentValues_ = initializeWithDescentDirection(Qstar, M, v, theta, 1e-2);
            } else {
                std::cout << "Solution verified at level p = " << p << std::endl;
                certificateResults_.Yopt = S;
                certificateResults_.Lambda = Lambda;
                certificateResults_.xhat = RoundSolutionS();
                auto t7 = std::chrono::high_resolution_clock::now();
                certificateResults_.total_computation_time = (std::chrono::duration<double, std::milli> (t7 - t6).count());
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
    inline gtsam::Pose2 addNoiseToPose2(
                const gtsam::Pose2& guess,
                double sigma_x     = 1,
                double sigma_y     = 1,
                double sigma_theta = M_PI / 180.0 * 1)
        {
            static std::mt19937 rng{std::random_device{}()};
            std::normal_distribution<double> N(0.0, 1.0);

            Eigen::Vector3d xi;
            xi << sigma_x * N(rng),
                    sigma_y * N(rng),
                    sigma_theta * N(rng);

            // Manifold-aware perturbation:  noisy = Exp(xi) ∘ guess
            return guess.retract(xi);     // or: gtsam::Pose2::Expmap(xi) * guess;
        }

        inline gtsam::Pose2 randomPose2(
                double min_x      = -1,  double max_x      =  1,
                double min_y      = -1,  double max_y      =  1,
                double min_theta  = -M_PI,  double max_theta  =   M_PI)
        {
            static std::mt19937 rng{std::random_device{}()};

            std::uniform_real_distribution<double> Ux(min_x,     max_x);
            std::uniform_real_distribution<double> Uy(min_y,     max_y);
            std::uniform_real_distribution<double> Ut(min_theta, max_theta);

            double x = Ux(rng);
            double y = Uy(rng);
            double t = Ut(rng);

            return gtsam::Pose2(x, y, t);
        }

        template <typename VALUE>
        bool read_incremental_g2o(string line, std::unordered_map<size_t, size_t> &poseIndex, size_t & nextPoseIdx,VALUE & Pose, std::shared_ptr<noiseModel::Diagonal> &noise , int &id1, int & id2)
        {
            string token;
            std::stringstream ss(line);
            ss >> token;
            if (token == "EDGE_SE2")
            {
                if constexpr (d == 2)
                {
                    // Parse g2o: id1 id2 dx dy dtheta sigmas... (we use information matrix entries)
                    double dx, dy, dtheta, I11, I12, I13, I22, I23, I33, tau, kappa;
                    ss >> id1 >> id2 >> dx >> dy >> dtheta
                       >> I11 >> I12 >> I13 >> I22 >> I23 >> I33;

                    if (!poseIndex.count(id1)) poseIndex[id1] = nextPoseIdx++;
                    if (!poseIndex.count(id2)) poseIndex[id2] = nextPoseIdx++;
                    // Build information matrix and noise model
                    Eigen::Matrix3d info;
                    info << I11, I12, I13,
                            I12, I22, I23,
                            I13, I23, I33;
                    //        auto noise = gtsam::noiseModel::Gaussian::Information(info);


                    Eigen::Matrix<Scalar, 2, 2> TranInfo;
                    TranInfo << I11, I12, I12, I22;
                    tau = 2 / TranInfo.inverse().trace();

                    kappa = I33;

                    Pose = VALUE(dx,dy,dtheta);

                    /* noise model goes into problem class */
                    Vector sigmas = Vector::Zero(3);
                    sigmas.head(2).setConstant(std::sqrt(1.0 / (1 * tau)));
                    sigmas.tail(2).setConstant(std::sqrt(1.0 / (1 * kappa)));
                    noise = noiseModel::Diagonal::Sigmas(sigmas);

                    return false;
                }
                return true;
            }
            if (token == "EDGE_SE3:QUAT")
            {
                if constexpr(d == 3)
                {
                    // Parse g2o: id1 id2 dx dy dtheta sigmas... (we use information matrix entries)
                    double dx, dy, dz, dtheta, dqx, dqy, dqz, dqw, I11, I12, I13, I14, I15, I16,
                        I22, I23, I24, I25, I26, I33, I34, I35, I36, I44, I45, I46, I55, I56, I66, tau, kappa;
                    ss >> id1 >> id2 >> dx >> dy >> dz >> dqx >> dqy >> dqz >> dqw >> I11 >>
                                I12 >> I13 >> I14 >> I15 >> I16 >> I22 >> I23 >> I24 >> I25 >> I26 >>
                                I33 >> I34 >> I35 >> I36 >> I44 >> I45 >> I46 >> I55 >> I56 >> I66;

                    if (!poseIndex.count(id1)) poseIndex[id1] = nextPoseIdx++;
                    if (!poseIndex.count(id2)) poseIndex[id2] = nextPoseIdx++;
                    // Build information matrix and noise model
                    // Raw measurements
                    Matrix t = Eigen::Matrix<Scalar, 3, 1>(dx, dy, dz);
                    Matrix R =
                            Eigen::Quaternion<Scalar>(dqw, dqx, dqy, dqz).toRotationMatrix();

                    // Compute precisions

                    // Compute and store the optimal (information-divergence-minimizing) value
                    // of the parameter tau
                    Eigen::Matrix<Scalar, 3, 3> TranInfo;
                    TranInfo << I11, I12, I13, I12, I22, I23, I13, I23, I33;
                    tau = 3 / TranInfo.inverse().trace();

                    Eigen::Matrix<Scalar, 3, 3> RotInfo;
                    RotInfo << I44, I45, I46, I45, I55, I56, I46, I56, I66;
                    kappa = 3 / (2 * RotInfo.inverse().trace());

                    Pose = VALUE(Rot3(R),t);

                    /* noise model goes into problem class */
                    Vector sigmas = Vector::Zero(6);
                    sigmas.tail(3).setConstant(std::sqrt(1.0 / (1 * tau)));
                    sigmas.head(3).setConstant(std::sqrt(1.0 / (1 * kappa)));
                    noise = noiseModel::Diagonal::Sigmas(sigmas);

                    return false;
                }
                return true;
            }
            return true;
        }


        template <size_t DCASE, typename V>
        struct PriorBuilder {
            static void apply(V&, SharedNoiseModel&) {
                static_assert(DCASE==2||DCASE==3, "Unsupported dimension");
            }
        };

        // 2D specialization
        template <typename V>
        struct PriorBuilder<2, V> {
            static void apply(V& pose, SharedNoiseModel& noise) {
                pose  = V(0.0, 0.0, 0.0);
                noise = noiseModel::Diagonal::Sigmas(Vector3(0.3, 0.3, 0.1));
            }
        };

        // 3D specialization
        template <typename V>
        struct PriorBuilder<3, V> {
            static void apply(V& pose, SharedNoiseModel& noise) {
                pose = V::Identity();
                Vector6 sigmas;
                sigmas << 0.3, 0.3, 0.3,   // σ_x, σ_y, σ_z
                          0.1, 0.1, 0.1;   // σ_roll, σ_pitch, σ_yaw
                noise = noiseModel::Diagonal::Sigmas(sigmas);
            }
        };

    template <typename VALUE>
    void setPrior(VALUE& pose, SharedNoiseModel& noise) {
        PriorBuilder<d, VALUE>::apply(pose, noise);
    }

        void ExportPose2(const Values& poses, const std::string& path, bool g2o) {
            if (g2o) {
                // Write g2o: keys must be integers starting at 0
                //        gtsam::writeG2o(poses, path + ".g2o");
            } else {
                std::ofstream file(path + ".tum");
                if (!file) {
                    std::cerr << "Error opening " << path << ".tum" << std::endl;
                    return;
                }
                // TUM format: timestamp tx ty tz qx qy qz qw
                // Use integer key as timestamp
                size_t id = 0;
                for (auto key : poses.keys()) {
                    Pose2 p = poses.at<Pose2>(key);
                    double tx = p.x();
                    double ty = p.y();
                    // embed into 3D: z = 0, quaternion from yaw
                    double theta = p.theta();
                    // Construct quaternion about Z-axis
                    Eigen::AngleAxisd aa(theta, Eigen::Vector3d::UnitZ());
                    Eigen::Quaterniond q(aa);
                    file << id << " "
                         << tx << " " << ty << " " << 0.0 << " "
                         << q.x() << " " << q.y() << " " << q.z() << " " << q.w() << "\n";
                    id++;
                }
                file.close();
            }
        }

        // Export Pose3 values: .g2o if g2o=true, else TUM-style (.tum)
        void ExportPose3(const Values& poses, const std::string& path, bool g2o) {
            if (g2o) {
                //        gtsam::writeG2o(poses, path + ".g2o");
            } else {
                std::ofstream file(path + ".tum");
                if (!file) {
                    std::cerr << "Error opening " << path << ".tum" << std::endl;
                    return;
                }
                // TUM format: timestamp tx ty tz qx qy qz qw
                size_t id = 0;
                for (auto key : poses.keys()) {
                    Pose3 p = poses.at<Pose3>(key);
                    auto t = p.translation();
                    Eigen::Quaterniond q(p.rotation().toQuaternion());
                    file << id << " "
                         << t.x() << " " << t.y() << " " << t.z() << " "
                         << q.x() << " " << q.y() << " " << q.z() << " " << q.w() << "\n";
                    id++;
                }
                file.close();
            }
        }

        // Convenience wrappers
        // Export 2D poses
        void ExportValues2D(const Values& poses, const std::string& path, bool g2o) {
            ExportPose2(poses, path, g2o);
        }

        // Export 3D poses
        void ExportValues3D(const Values& poses, const std::string& path, bool g2o) {
            ExportPose3(poses, path, g2o);
        }

        void ExportValues(const Values& poses, const std::string& path, bool g2o)
        {
            if constexpr (d == 2) {
                ExportValues2D(poses, path, g2o);
            } else if constexpr (d == 3)
            {
                ExportValues3D(poses, path, g2o);
            }
        }
};

// Convenience aliases
using CertifiableIncrementalPGO2 = CertifiableIncrementalPGO<2>;
using CertifiableIncrementalPGO3 = CertifiableIncrementalPGO<3>;

#endif //CERTIFIABLE_INCREMENTAL_PGO_H
