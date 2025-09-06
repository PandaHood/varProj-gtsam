//
// Created by nikolas on 5/5/25.
//

#ifndef CERTIFIABLEPROBLEMOPTS_H
#define CERTIFIABLEPROBLEMOPTS_H
#include "Certifiable_problem.h"

using namespace std;
using namespace gtsam;

namespace gtsam
{
    /**
     * @brief Configuration options for certifiable problem routines.
     *
     * Contains parameters for the Levenberg–Marquardt optimizer as well as
     * settings for the fast verification (eigenvalue) step.
     */
    struct CertifiableProblemOpts
    {
        /// Levenberg–Marquardt optimizer parameters (default: Ceres defaults).
        LevenbergMarquardtParams lmParams = LevenbergMarquardtParams::CeresDefaults();

        /// Maximum number of LM iterations before termination.
        int maxIterations = 200;

        /// Relative error tolerance for LM convergence.
        double relativeErrorTol = 1e-15;

        /// Absolute error tolerance for LM convergence.
        double absoluteErrorTol = 1e-15;

        /// Level of verbosity for LM output (SUMMARY, SILENT, etc.).
        LevenbergMarquardtParams::VerbosityLM verbosityLM = LevenbergMarquardtParams::SUMMARY;

        /// Use absolute regularization term eat, or relative one
        /// Usually use relative eta in range-related example
        bool useAbsoluteEta = false;

        /// Regularization parameter η for certificate matrix M = S + η·I.
        Scalar eta = 1e-3;

        /// Lower bound for the eta
        const Scalar MIN_CERT_ETA = 1e-4;

        /// Upper bound for the eta
        const Scalar MAX_CERT_ETA = 1e-1;

        /// Relative ratio for the eta, i.e., eta = objective_value * REL_CERT_ETA
        const Scalar REL_CERT_ETA =  4.04597e-8;

        /// Block size (number of eigenvectors) for LOBPCG verification.
        size_t nx = 4;

        /// Maximum number of iterations for the LOBPCG solver.
        size_t max_iters = 100;

        /// Maximum fill factor for the ILDL preconditioner.
        Scalar max_fill_factor = 3;

        /// Drop tolerance for the ILDL preconditioner.
        Scalar drop_tol = 1;

        // Initialization type
        enum class InitType { Odom, Random, LocalSearch };
        InitType initType = InitType::Random;
    };
}  // namespace gtsam



#endif //CERTIFIABLEPROBLEMOPTS_H
