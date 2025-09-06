//
// Created by Nikolas on 2/3/25.
//
/* ----------------------------------------------------------------------------

 * GTSAM Copyright 2010-2019, Georgia Tech Research Corporation,
 * Atlanta, Georgia 30332-0415
 * All Rights Reserved
 * Authors: Frank Dellaert, et al. (see THANKS for the full author list)

 * See LICENSE for the license information

 * -------------------------------------------------------------------------- */

/**
 * @file LandmarkFactor.cpp
 * @date
 * @author Nikolas Sanderson, Jason Xu
 * @brief
 */

#include <gtsam/base/timing.h>
#include <gtsam/nonlinear/NonlinearFactorGraph.h>
#include "LandmarkFactor.h"
#include <iostream>
#include <vector>

using namespace std;
namespace gtsam
{

    //******************************************************************************
    /**
     * @brief Constructor for a lifted landmark factor.
     * @tparam d  Ambient dimension (2 or 3).
     * @param j1    Key of the observing pose variable.
     * @param j2    Key of the landmark variable.
     * @param T12   Measured translation vector V_ from pose to landmark.
     * @param p     Relaxation rank.
     * @param model Shared noise model; must have dimension == p.
     * @throws std::invalid_argument if noise model dimension ≠ p.
     */
    template <size_t d>
    LiftedLandmarkFactor<d>::LiftedLandmarkFactor(Key j1, Key j2,
                                                  const Vector &T12,
                                                  size_t p,
                                                  const SharedNoiseModel &model)
      : NoiseModelFactorN<LiftedPoseDP, Vector>(model, j1, j2),
        V_(T12),                     // store measured offset
        p_(p),                       // store relaxation rank
        d_(d),                       // ambient dimension
        pd_(p * d)                   // product p*d for convenience
    {
        // ensure noise model matches expected dimension
        if (noiseModel()->dim() != p_)
            throw std::invalid_argument(
                "LandmarkFactor: model with incorrect dimension.");
    }

    //******************************************************************************
    /**
     * @brief Print this factor’s contents for debugging.
     * @param s             Prefix string.
     * @param keyFormatter  Function to convert Keys → printable strings.
     */
    template <size_t d>
    void LiftedLandmarkFactor<d>::print(const std::string &s,
                                        const KeyFormatter &keyFormatter) const
    {
        std::cout << s
                  << "LandmarkFactor<" << p_ << " * " << d
                  << ">(" << keyFormatter(key<1>()) << ","
                  << keyFormatter(key<2>()) << ")\n";
        traits<Matrix>::Print(V_, "  V: ");         // print measured vector
        noiseModel_->print("  noise model: ");      // print noise details
    }

    //******************************************************************************
    /**
     * @brief Check equality with another factor.
     * @param expected  Other factor to compare.
     * @param tol       Numeric tolerance.
     * @return          True if same type, keys, rank, and measurement.
     */
    template <size_t d>
    bool LiftedLandmarkFactor<d>::equals(const NonlinearFactor &expected,
                                         double tol) const
    {
        auto e = dynamic_cast<const LiftedLandmarkFactor *>(&expected);
        return e != nullptr
            && NoiseModelFactorN<LiftedPoseDP, Vector>::equals(*e, tol)
            && p_ == e->p_
            && V_ == e->V_;
    }

    //******************************************************************************
    /**
     * @brief Compute Jacobians of the error w.r.t. pose and landmark.
     * @param Q1  Current lifted pose variable.
     * @param L1  Current landmark variable vector.
     * @param H1  Optional output Jacobian w.r.t. Q1.
     * @param H2  Optional output Jacobian w.r.t. L1.
     */
    template <size_t d>
    void LiftedLandmarkFactor<d>::fillJacobians(const LiftedPoseDP &Q1,
                                                const Vector &L1,
                                                OptionalMatrixType H1,
                                                OptionalMatrixType H2) const
    {
        const StiefelManifoldKP Y1 = Q1.get_Y();    // rotation part
        const Vector t1 = Q1.get_t();              // translation from pose
        const Vector t2 = L1;                      // landmark translation
        const size_t St_OutDim = p_ * d_;          // flattened rotation block size
        const size_t St_dim    = StiefelManifoldKP::Dimension(d_, p_);
        const size_t t_OutDim  = p_;               // translation output dim
        const size_t t_dim     = p_;               // translation input dim

        if (H1) {
            // Jacobian w.r.t. the lifted pose Q1: [derror/dY | derror/dt1]
            H1->resize(p_, St_dim + t_dim);
            H1->setZero();

            // build derivative of translation error w.r.t. Y1
            Matrix dt_dY = Matrix::Zero(t_OutDim, St_OutDim);
            for (int j = 0; j < d_; ++j) {
                dt_dY.block(0, j * p_, p_, p_) -= V_(j)
                    * Matrix::Identity(p_, p_);
            }
            H1->block(0, 0, t_OutDim, St_dim) = dt_dY * Y1.G_;
            H1->block(0, St_dim, t_OutDim, t_dim) =
                -Matrix::Identity(t_OutDim, t_dim);
        }

        if (H2) {
            // Jacobian w.r.t. landmark L1: identity in translation
            H2->resize(p_, t_dim);
            H2->setZero();
            H2->block(0, 0, t_OutDim, t_dim) =
                Matrix::Identity(t_OutDim, t_dim);
        }
    }

    //******************************************************************************
    /**
     * @brief Evaluate the error vector for this factor.
     * @param Q1   Current lifted pose variable.
     * @param L1   Current landmark variable vector.
     * @param H1   Optional Jacobian w.r.t. Q1.
     * @param H2   Optional Jacobian w.r.t. L1.
     * @return     Error = t2 - t1 - Y1.matrix() * V_.
     */
    template <size_t d>
    Vector LiftedLandmarkFactor<d>::evaluateError(const LiftedPoseDP &Q1,
                                                  const Vector &L1,
                                                  OptionalMatrixType H1,
                                                  OptionalMatrixType H2) const
    {
        const StiefelManifoldKP Y1 = Q1.get_Y();
        // translation error: measured minus predicted
        Vector error_trans = L1 - Q1.get_t() - Y1.matrix() * V_;

        // fill H1/H2 if requested
        this->fillJacobians(Q1, L1, H1, H2);

        return error_trans;
    }

    template <size_t d>
    void LiftedLandmarkFactor<d>::computeHessian(std::map<std::pair<Key,Key>,Matrix> &HMap) const
    {
        Eigen::MatrixXd I = Eigen::MatrixXd::Identity(d, d);
        Eigen::MatrixXd H = Eigen::MatrixXd::Zero(3*d, 3*d);

        // Top row blocks: [ I , −I ,  t_ij·1ᵀ ]
        H.block(0,   0, d, d) =  I;
        H.block(0,   d, d, d) = -I;
        H.block(0, 2*d, d, 1) =  V_;

        // Middle row blocks: [ −I , I , −t_ij·1ᵀ ]
        H.block(d,   0, d, d) = -I;
        H.block(d,   d, d, d) =  I;
        H.block(d, 2*d, d, 1) = -V_;

        // Bottom row blocks: [ 1·t_ijᵀ , −1·t_ijᵀ , t_ij·t_ijᵀ ]
        H.block(2*d,   0, 1, d) =  V_.transpose();
        H.block(2*d,   d, 1, d) = -V_.transpose();
        H.block(2*d, 2*d, d, d) =  V_ * V_.transpose();
        double tau_sigmas = noiseModel_->sigmas()(0);
        double tau  = 1.0 / (1*tau_sigmas*tau_sigmas);

        // Store
        // finally store half of the reordered Hessian
        if (HMap.count({key1(), key2()}) == 0)
            HMap[{key1(),key2()}] = tau * H;
        else
            HMap[{key1(),key2()}] += tau * H;
    }
    // Explicit instantiations for d=2 and d=3
    template class LiftedLandmarkFactor<2>;
    template class LiftedLandmarkFactor<3>;

     //******************************************************************************
    /**
     * @brief Constructor for a lifted landmark factor.
     * @tparam d  Ambient dimension (2 or 3).
     * @param j1    Key of the observing pose variable.
     * @param j2    Key of the landmark variable.
     * @param T12   Measured translation vector V_ from pose to landmark.
     * @param p     Relaxation rank.
     * @param model Shared noise model; must have dimension == p.
     * @throws std::invalid_argument if noise model dimension ≠ p.
     */
    template <size_t d>
    BAFactor<d>::BAFactor(Key j1, Key j2 ,Key j3, const Vector &T12 , size_t p,
                     const SharedNoiseModel &model)
      : NoiseModelFactorN<ScaledStiefelKPD,Vector, Vector>(model, j1, j2,j3),
        V_(T12),                     // store measured offset
        p_(p),                       // store relaxation rank
        d_(d),                       // ambient dimension
        pd_(p * d)                   // product p*d for convenience
    {
        // ensure noise model matches expected dimension
        if (noiseModel()->dim() != p_)
            throw std::invalid_argument(
                "LandmarkFactor: model with incorrect dimension.");
    }

    //******************************************************************************
    /**
     * @brief Print this factor’s contents for debugging.
     * @param s             Prefix string.
     * @param keyFormatter  Function to convert Keys → printable strings.
     */
    template <size_t d>
    void BAFactor<d>::print(const std::string &s,
                                        const KeyFormatter &keyFormatter) const
    {
        std::cout << s
                  << "LandmarkFactor<" << p_ << " * " << d
                  << ">(" << keyFormatter(key<1>()) << ","
                  << keyFormatter(key<2>()) << ")\n";
        traits<Matrix>::Print(V_, "  V: ");         // print measured vector
        noiseModel_->print("  noise model: ");      // print noise details
    }

    //******************************************************************************
    /**
     * @brief Check equality with another factor.
     * @param expected  Other factor to compare.
     * @param tol       Numeric tolerance.
     * @return          True if same type, keys, rank, and measurement.
     */
    template <size_t d>
    bool BAFactor<d>::equals(const NonlinearFactor &expected,
                                         double tol) const
    {
        auto e = dynamic_cast<const BAFactor *>(&expected);
        return e != nullptr
            && NoiseModelFactorN<ScaledStiefelKPD,Vector, Vector>::equals(*e, tol)
            && p_ == e->p_
            && V_ == e->V_;
    }

    //******************************************************************************
    /**
     * @brief Compute Jacobians of the error w.r.t. pose and landmark.
     * @param Q1  Current lifted pose variable.
     * @param L1  Current landmark variable vector.
     * @param H1  Optional output Jacobian w.r.t. Q1.
     * @param H2  Optional output Jacobian w.r.t. L1.
     */
    template <size_t d>
    void BAFactor<d>::fillJacobians(const ScaledStiefelKPD& P1,const Vector& P2 ,const Vector &L1, OptionalMatrixType H1, OptionalMatrixType H2, OptionalMatrixType H3) const
    {
        // Shapes
          const Matrix Y = P1.toMatrix();                // p x k
          const size_t p = static_cast<size_t>(Y.rows());
          const size_t k = static_cast<size_t>(Y.cols());
          const double eps = 1e-12;

        // Use the variable’s exact Stiefel and cached basis (no re-SVD, no drift)
        const Matrix& Rhat = P1.R();                    // p x k (exact on Stiefel)
        const Matrix& G_R  = P1.stiefelBasis();         // (p*k) x dSt
        const size_t dSt   = static_cast<size_t>(G_R.cols());

        // Build the scale vector for w = s ⊙ V
        Vector svec = P1.perColumn()
                      ? P1.scales()                                        // size k
                      : Vector::Constant((int)Rhat.cols(), P1.scales()(0)); // broadcast scalar


          // w = s ⊙ V (works for scalar too, since all s_j == s)
          const Vector w = svec.cwiseProduct(V_.head((int)k));

          // M(w) : (p x pk) such that M(w)*vec(ΔR) = ΔR * w
          Matrix M = Matrix::Zero((int)p, (int)(p*k));
          for (size_t j = 0; j < k; ++j) {
            M.block(0, (int)(j*p), (int)p, (int)p).setIdentity();
            M.block(0, (int)(j*p), (int)p, (int)p) *= -w((int)j); // negative sign included
          }

          // H1: w.r.t. ScaledStiefel local coordinates = [scale part | Stiefel part]
          if (H1) {
            // ds = #scale coordinates: either 1 (scalar) or k (per-column).
            // We can deduce it from the manifold dimension: dim = ds + dSt
            const size_t ds = static_cast<size_t>(P1.dim()) - dSt;

            H1->resize((int)p, (int)(ds + dSt));
            H1->setZero();

            // ---- scale part ----
            if (ds == 1) {
              // scalar scale: de/ds = -(R V) = -(Y V)/s
              const double s_scalar = std::max(eps, svec(0));
              const Vector deds = -(Y * V_.head((int)k)) / s_scalar;  // p x 1
              H1->block(0, 0, (int)p, 1) = deds;
            } else {
              // per-column: de/ds_j = - V_j * r_j  (use Rhat.col(j))
              for (size_t j = 0; j < k; ++j) {
                H1->block(0, (int)j, (int)p, 1) = -V_((int)j) * Rhat.col((int)j);
              }
            }

            // ---- Stiefel part ----
            // Δe = - ΔR (s ⊙ V)  => vec(Δe) = - ( (w^T ⊗ I_p) vec(ΔR) )
            // Our M already maps vec(ΔR) -> ΔR w, so:
            H1->block(0, (int)ds, (int)p, (int)dSt) = M * G_R;
          }

          // H2: w.r.t. P2 (subtracts directly)
          if (H2) {
            H2->resize((int)p, (int)p);
            *H2 = -Matrix::Identity((int)p, (int)p);
          }

          // H3: w.r.t. L1 (adds directly)
          if (H3) {
            H3->resize((int)p, (int)p);
            *H3 = Matrix::Identity((int)p, (int)p);
          }
    }

    //******************************************************************************
    /**
     * @brief Evaluate the error vector for this factor.
     * @param Q1   Current lifted pose variable.
     * @param L1   Current landmark variable vector.
     * @param H1   Optional Jacobian w.r.t. Q1.
     * @param H2   Optional Jacobian w.r.t. L1.
     * @return     Error = t2 - t1 - Y1.matrix() * V_.
     */
    template <size_t d>
    Vector BAFactor<d>::evaluateError(const ScaledStiefelKPD& P1,const Vector& P2 ,const Vector &L1, OptionalMatrixType H1, OptionalMatrixType H2, OptionalMatrixType H3) const
    {// translation error: measured minus predicted
        Vector error_trans = L1 - P2 - P1.toMatrix() * V_;

        // fill H1/H2 if requested
        this->fillJacobians(P1,P2, L1, H1, H2,H3);

        return error_trans;
    }

    // Explicit instantiations for d=2 and d=3
    template class BAFactor<2>;
    template class BAFactor<3>;
} // namespace gtsam


