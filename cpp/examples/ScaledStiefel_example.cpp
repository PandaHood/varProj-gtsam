//
// Created by nikolas on 9/3/25.
//

// ScaledStiefel_SelfCheck.cpp
#include <iostream>
#include <iomanip>
#include <string>
#include <Eigen/Dense>

// ⚠️ Adjust this include to where your class lives:
#include "../StiefelManifold.h"

using std::cout;
using std::endl;
using std::string;

template<int K, int P, bool PerCol>
using SS = gtsam::ScaledStiefel<K, P, PerCol>;

using SS33Scalar = SS<3, 3, false>;
using SS33PerCol = SS<3, 3, true>;
using SSDynScalar = SS<Eigen::Dynamic, Eigen::Dynamic, false>;
using SSDynPerCol = SS<Eigen::Dynamic, Eigen::Dynamic, true>;

// ---------- helpers ----------
static double frob(const Eigen::MatrixXd& A) { return A.norm(); }
static double offDiagFrob(const Eigen::MatrixXd& A) {
  Eigen::MatrixXd B = A;
  B.diagonal().setZero();
  return B.norm();
}
static bool near(double a, double b, double tol) { return std::abs(a - b) <= tol; }
static bool matNear(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B, double tol) {
  return (A - B).norm() <= tol;
}

struct Suite {
  int passed{0}, total{0};
  void report(const string& name, bool ok, const string& msg = "") {
    ++total;
    cout << (ok ? "[PASS] " : "[FAIL] ") << name;
    if (!ok && !msg.empty()) cout << " :: " << msg;
    cout << endl;
    if (ok) ++passed;
  }
  void summary() const {
    cout << "\n===== Summary: " << passed << " / " << total << " passed =====\n";
  }
};

// ---------- checks ----------
void check_CompileTimeDimension(Suite& s) {
  bool ok1 = (SS33Scalar::dimension == 4);  // 1 + dSt(3,3) == 1 + 3
  bool ok2 = (SS33PerCol::dimension == 6);  // 3 + dSt(3,3) == 3 + 3
  bool ok3 = (SSDynScalar::dimension == Eigen::Dynamic);
  bool ok4 = (SSDynPerCol::dimension == Eigen::Dynamic);
  s.report("CompileTimeDimension(SS33Scalar)", ok1, "expected 4");
  s.report("CompileTimeDimension(SS33PerCol)", ok2, "expected 6");
  s.report("CompileTimeDimension(SSDynScalar)", ok3, "expected Dynamic");
  s.report("CompileTimeDimension(SSDynPerCol)", ok4, "expected Dynamic");
}

void check_RuntimeDim(Suite& s) {
  SSDynScalar a(3,3,false);
  SSDynPerCol b(3,3,true);
  bool oka = (a.dim() == (1 + 3)); // scalar: 1 + dSt(3,3)=1+3
  bool okb = (b.dim() == (3 + 3)); // per-col: k + dSt(3,3)=3+3
  s.report("RuntimeDim(SSDynScalar 3x3)", oka, "expected 4");
  s.report("RuntimeDim(SSDynPerCol 3x3)", okb, "expected 6");
}

void check_ProjectToManifold_Scalar(Suite& s) {
  const int p = 3, k = 3;
  Eigen::MatrixXd G = Eigen::MatrixXd::Random(p, k);
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(G, Eigen::ComputeThinU | Eigen::ComputeThinV);
  Eigen::MatrixXd R = svd.matrixU() * svd.matrixV().transpose();

  const double s_val = 2.5;
  const Eigen::MatrixXd Y = s_val * R;

  auto X = SSDynScalar::projectToManifold(Y, /*perColumn=*/false);
  const Eigen::MatrixXd Yhat = X.toMatrix();

  const double tol = 1e-10;
  bool okRec = near(frob(Y - Yhat), 0.0, tol);
  s.report("ProjectToManifold_Scalar: reconstruct", okRec, "||Y - Yhat|| too large");

  const Eigen::MatrixXd Gram = Yhat.transpose() * Yhat;
  const Eigen::MatrixXd I = Eigen::MatrixXd::Identity(k, k);
  bool okGram = near(frob(Gram - (s_val * s_val) * I), 0.0, tol);
  s.report("ProjectToManifold_Scalar: Gram=s^2 I", okGram, "Gram mismatch");
}

void check_ProjectToManifold_PerColumn(Suite& s) {
  const int p = 3, k = 3;
  Eigen::MatrixXd G = Eigen::MatrixXd::Random(p, k);
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(G, Eigen::ComputeThinU | Eigen::ComputeThinV);
  Eigen::MatrixXd R = svd.matrixU() * svd.matrixV().transpose();

  Eigen::Vector3d sc; sc << 1.2, 2.0, 3.5;
  const Eigen::MatrixXd Y = R * sc.asDiagonal();

  auto X = SSDynPerCol::projectToManifold(Y, /*perColumn=*/true);
  const Eigen::MatrixXd Yhat = X.toMatrix();

  const double tol = 1e-10;
  bool okRec = near(frob(Y - Yhat), 0.0, tol);
  s.report("ProjectToManifold_PerColumn: reconstruct", okRec, "||Y - Yhat|| too large");

  const Eigen::MatrixXd Gram = Yhat.transpose() * Yhat;
  bool okOff = near(offDiagFrob(Gram), 0.0, tol);
  bool okDiagPos = true;
  for (int i = 0; i < k; ++i) okDiagPos = okDiagPos && (Gram(i,i) > 0.0);
  s.report("ProjectToManifold_PerColumn: off-diag≈0", okOff, "off-diagonals not ~0");
  s.report("ProjectToManifold_PerColumn: diag>0", okDiagPos, "diagonal not > 0");
}

void check_RandomStructure(Suite& s) {
  const double tol = 1e-10;

  auto Xs = SSDynScalar::Random(/*seed=*/1234, /*k=*/3, /*p=*/3, /*perColumn=*/false);
  const Eigen::MatrixXd Ys = Xs.toMatrix();
  const Eigen::MatrixXd Grams = Ys.transpose() * Ys;
  bool okS_off = near(offDiagFrob(Grams), 0.0, tol);
  bool okS_diag = true; for (int i=0;i<Grams.rows();++i) okS_diag = okS_diag && (Grams(i,i) > 0.0);
  s.report("RandomStructure_Scalar: off-diag≈0", okS_off);
  s.report("RandomStructure_Scalar: diag>0", okS_diag);

  auto Xd = SSDynPerCol::Random(/*seed=*/42, /*k=*/3, /*p=*/3, /*perColumn=*/true);
  const Eigen::MatrixXd Yd = Xd.toMatrix();
  const Eigen::MatrixXd Gramd = Yd.transpose() * Yd;
  bool okD_off = near(offDiagFrob(Gramd), 0.0, tol);
  bool okD_diag = true; for (int i=0;i<Gramd.rows();++i) okD_diag = okD_diag && (Gramd(i,i) > 0.0);
  s.report("RandomStructure_PerColumn: off-diag≈0", okD_off);
  s.report("RandomStructure_PerColumn: diag>0", okD_diag);
}

void check_RetractLocal_RoundTrip(Suite& s) {
  const double tol = 5e-7;

  {
    SSDynScalar X0(3,3,false);
    const size_t d = X0.dim();
    Eigen::VectorXd xi = 1e-3 * Eigen::VectorXd::Random((int)d);
    auto X1 = X0.retract(xi);
    auto xi_hat = X0.localCoordinates(X1);
    bool okSize = (xi_hat.size() == xi.size());
    bool okFirstOrder = ((xi_hat - xi).norm() <= tol);
    s.report("RetractLocal_RoundTrip_Scalar: size", okSize);
    s.report("RetractLocal_RoundTrip_Scalar: first-order", okFirstOrder);
  }

  {
    SSDynPerCol X0(3,3,true);
    const size_t d = X0.dim();
    Eigen::VectorXd xi = 1e-3 * Eigen::VectorXd::Random((int)d);
    auto X1 = X0.retract(xi);
    auto xi_hat = X0.localCoordinates(X1);
    bool okSize = (xi_hat.size() == xi.size());
    bool okFirstOrder = ((xi_hat - xi).norm() <= tol);
    s.report("RetractLocal_RoundTrip_PerColumn: size", okSize);
    s.report("RetractLocal_RoundTrip_PerColumn: first-order", okFirstOrder);
  }
}

void check_ProjectorRespectsGram(Suite& s) {
  const int p = 3, k = 3;
  const double tol = 1e-10;

  Eigen::MatrixXd Y = Eigen::MatrixXd::Random(p, k) + 0.5 * Eigen::MatrixXd::Ones(p, k);

  auto Xs = SSDynScalar::projectToManifold(Y, false);
  auto Xd = SSDynPerCol::projectToManifold(Y, true);

  const Eigen::MatrixXd Ys = Xs.toMatrix();
  const Eigen::MatrixXd Yd = Xd.toMatrix();

  const Eigen::MatrixXd Grams = Ys.transpose() * Ys;
  const Eigen::MatrixXd Gramd = Yd.transpose() * Yd;

  bool okS = near(offDiagFrob(Grams), 0.0, tol);
  bool okD_off = near(offDiagFrob(Gramd), 0.0, tol);
  bool okD_diag = true; for (int i=0;i<k;++i) okD_diag = okD_diag && (Gramd(i,i) > 0.0);

  s.report("ProjectorRespectsGram: scalar off-diag≈0", okS);
  s.report("ProjectorRespectsGram: per-col off-diag≈0", okD_off);
  s.report("ProjectorRespectsGram: per-col diag>0", okD_diag);
}

void check_FixedSizeConstruction(Suite& s) {
  SS33Scalar Xs;
  Eigen::MatrixXd Ys = Xs.toMatrix();
  bool ok1 = (Ys.rows() == 3 && Ys.cols() == 3);

  SS33PerCol Xp;
  Eigen::MatrixXd Yp = Xp.toMatrix();
  bool ok2 = (Yp.rows() == 3 && Yp.cols() == 3);

  s.report("FixedSizeConstruction(SS33Scalar)", ok1, "expected 3x3");
  s.report("FixedSizeConstruction(SS33PerCol)", ok2, "expected 3x3");
}

void check_ProjectToManifold_Idempotent(Suite& s) {
  const int p=3,k=3; const double tol=1e-12;

  // scalar
  auto Xs0 = SSDynScalar::Random(1, k, p, false);
  auto Xs1 = SSDynScalar::projectToManifold(Xs0.toMatrix(), false);
  bool okS = (Xs0.toMatrix() - Xs1.toMatrix()).norm() <= tol;
  s.report("Idempotent project (scalar)", okS);

  // per-column
  auto Xp0 = SSDynPerCol::Random(2, k, p, true);
  auto Xp1 = SSDynPerCol::projectToManifold(Xp0.toMatrix(), true);
  bool okP = (Xp0.toMatrix() - Xp1.toMatrix()).norm() <= tol;
  s.report("Idempotent project (per-col)", okP);
}

void check_LiftedDimsAndGram(Suite& s) {
  const int k=3; const double tol=1e-12;
  for (int p : {4,5,6,7}) {
    SSDynScalar X(k,p,false);
    bool okDim = (X.dim() == size_t(1 + k*(p-k) + k*(k-1)/2)); // = 1 + (3(p-3)+3) = 3p-2
    s.report("Lifted dim p="+std::to_string(p), okDim);

    auto Y = X.toMatrix();
    // ---- in check_LiftedDimsAndGram ----
    Eigen::MatrixXd Gram = (Y.transpose() * Y).eval();              // force dense
    Eigen::MatrixXd D    = Gram.diagonal().asDiagonal();            // materializes to dense
    Eigen::MatrixXd off  = Gram - D;                                // OK now
    bool okOff = (off.norm() <= tol);
    bool okPos = (Gram.diagonal().minCoeff() > 0.0);

    s.report("Lifted Gram off≈0 p="+std::to_string(p), okOff);
    s.report("Lifted Gram diag>0 p="+std::to_string(p), okPos);
  }
}

void check_NegativeScalarInput(Suite& s) {
  const int p=3,k=3; const double tol=1e-10;
  Eigen::MatrixXd G = Eigen::MatrixXd::Random(p,k);
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(G, Eigen::ComputeThinU | Eigen::ComputeThinV);
  Eigen::MatrixXd R = svd.matrixU() * svd.matrixV().transpose();
  double sneg = -2.7;
  Eigen::MatrixXd Y = sneg * R;

  auto X = SSDynScalar::projectToManifold(Y, false);
  auto Yhat = X.toMatrix();
  bool ok = (Yhat - Y).norm() <= tol;  // sign should be absorbed by R, s>0
  s.report("Negative scalar handled", ok);
}
void check_PerColumn_PermutationSign(Suite& s) {
  const int p=3,k=3; const double tol=1e-10;

  // Build clean Y = R diag(sc)
  Eigen::MatrixXd G = Eigen::MatrixXd::Random(p,k);
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(G, Eigen::ComputeThinU|Eigen::ComputeThinV);
  Eigen::MatrixXd R = svd.matrixU()*svd.matrixV().transpose();
  Eigen::Vector3d sc; sc << 0.8, 1.7, 3.1;
  Eigen::MatrixXd Y = R * sc.asDiagonal();

  // Apply a signed permutation to columns
  Eigen::Matrix3d P = Eigen::Matrix3d::Zero();
  std::array<int,3> perm{0,1,2};
  std::mt19937 rng(123);
  std::shuffle(perm.begin(), perm.end(), rng);
  for (int j=0;j<3;++j) P(perm[j], j) = (j%2==0 ? 1.0 : -1.0);
  Eigen::MatrixXd Yperm = Y * P;

  // Project and compare invariants
  auto X = SSDynPerCol::projectToManifold(Yperm, true);
  Eigen::MatrixXd Gram_in  = (Yperm.transpose()*Yperm).eval();
  Eigen::MatrixXd Gram_out = (X.toMatrix().transpose()*X.toMatrix()).eval();

  // Off-diagonals ~ 0 for both
  bool okOff = (offDiagFrob(Gram_out) <= tol);

  // Sorted diag (squared column norms) must match
  Eigen::Vector3d din = Gram_in.diagonal();
  Eigen::Vector3d dout= Gram_out.diagonal();
  std::array<double,3> ain{din(0),din(1),din(2)}, aout{dout(0),dout(1),dout(2)};
  std::sort(ain.begin(), ain.end());
  std::sort(aout.begin(), aout.end());
  double diff = std::sqrt( (ain[0]-aout[0])*(ain[0]-aout[0]) +
                           (ain[1]-aout[1])*(ain[1]-aout[1]) +
                           (ain[2]-aout[2])*(ain[2]-aout[2]) );
  bool okDiagMatch = (diff <= 1e-12);

  s.report("Per-col perm+sign invariance (Gram diag set)", okOff && okDiagMatch);
}

void check_Retract_StaysOnManifold(Suite& s) {
  const double tol = 1e-10;
  // Scalar
  {
    SSDynScalar X0(3,3,false);
    Eigen::VectorXd xi = Eigen::VectorXd::Random((int)X0.dim()); // large step
    auto X1 = X0.retract(xi);
    auto Y = X1.toMatrix();
    Eigen::MatrixXd Gram = Y.transpose()*Y;
    Eigen::MatrixXd D    = Gram.diagonal().asDiagonal(); // dense diag
    // off-diagonal ~0; diag > 0
    bool okOff = ((Gram - D).norm() <= tol);
    bool okPos = (Gram.diagonal().minCoeff() > 0.0);
    s.report("Retract stays manifold (scalar)", okOff && okPos);
  }
  // Per-column
  {
    SSDynPerCol X0(3,3,true);
    Eigen::VectorXd xi = Eigen::VectorXd::Random((int)X0.dim());
    auto X1 = X0.retract(xi);
    auto Y = X1.toMatrix();
    // ---- in check_Retract_StaysOnManifold (both scalar & per-col blocks) ----
    Eigen::MatrixXd Gram = (Y.transpose() * Y).eval();              // force dense
    Eigen::MatrixXd D    = Gram.diagonal().asDiagonal();
    bool okOff = ((Gram - D).norm() <= tol);
    bool okPos = (Gram.diagonal().minCoeff() > 0.0);

    s.report("Retract stays manifold (per-col)", okOff && okPos);
  }
}
void check_TinyScales(Suite& s) {
  const int p=3,k=3; const double tol=1e-10;
  SSDynPerCol X0(k,p,true);
  // Set scales extremely small via retraction on scale coords
  Eigen::VectorXd xi = Eigen::VectorXd::Zero((int)X0.dim());
  // For perColumn, the first k entries are scale directions
  xi.head(k).setConstant(-50.0 * X0.toMatrix().col(0).norm()); // big negative to shrink
  auto X1 = X0.retract(xi);
  auto Y = X1.toMatrix();
  Eigen::MatrixXd Gram = Y.transpose()*Y;
  bool okOff = ((Gram - Gram.diagonal().asDiagonal().toDenseMatrix()).norm() <= tol);
  bool okFinite = Y.allFinite();
  s.report("Tiny scales remain finite & structured", okOff && okFinite);
}
void check_LocalZero_And_Additivity(Suite& s) {
  const double tol = 1e-6;
  SSDynScalar X0(3,3,false);
  auto z = X0.localCoordinates(X0);
  s.report("localCoordinates(X,X)==0", z.norm() <= 1e-15);

  Eigen::VectorXd xi1 = 1e-3 * Eigen::VectorXd::Random((int)X0.dim());
  Eigen::VectorXd xi2 = 1e-3 * Eigen::VectorXd::Random((int)X0.dim());
  auto X1 = X0.retract(xi1);
  auto X12 = X1.retract(xi2);
  auto xi12_from_X0 = X0.localCoordinates(X12);
  bool ok = ((xi12_from_X0 - (xi1 + xi2)).norm() <= tol);
  s.report("First-order additivity", ok);
}
void check_NoiseRobustProjection(Suite& s) {
  const int p=3,k=3; const double tol=1e-6;
  // clean per-col Y
  Eigen::MatrixXd G = Eigen::MatrixXd::Random(p,k);
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(G, Eigen::ComputeThinU|Eigen::ComputeThinV);
  Eigen::MatrixXd R = svd.matrixU()*svd.matrixV().transpose();
  Eigen::Vector3d sc; sc << 1.0, 2.0, 3.0;
  Eigen::MatrixXd Y = R * sc.asDiagonal();

  Eigen::MatrixXd N = 1e-3 * Eigen::MatrixXd::Random(p,k);
  auto X = SSDynPerCol::projectToManifold(Y + N, true);
  auto Gram = X.toMatrix().transpose()*X.toMatrix();
  double off = (Gram.matrix() - Gram.diagonal().asDiagonal().toDenseMatrix()).norm();
  s.report("Noisy projection keeps Gram ~ diagonal", off <= tol);
}
double f_of(const SSDynScalar& X, const Eigen::MatrixXd& Y0) {
  Eigen::MatrixXd Y = X.toMatrix();
  return 0.5 * (Y - Y0).squaredNorm();
}
void check_ManifoldDirectionalDerivative(Suite& s) {
  const int k=3,p=5;
  SSDynScalar X0(k,p,false);
  Eigen::MatrixXd Y0 = Eigen::MatrixXd::Random(p,k);
  Eigen::VectorXd xi = Eigen::VectorXd::Random((int)X0.dim());
  xi.normalize();

  auto phi = [&](double t){
    return f_of( X0.retract(t*xi), Y0 );
  };
  double h = 1e-6;
  double dnum = (phi(h) - phi(-h)) / (2*h);
  // Sanity: derivative should be finite and not absurdly large
  bool ok = std::isfinite(dnum) && std::abs(dnum) < 1e6;
  s.report("Directional derivative finite (scalar,p=5)", ok);
}



// ---------- main ----------
int main() {
  std::cout.setf(std::ios::fixed); std::cout << std::setprecision(6);
  Suite suite;

  check_CompileTimeDimension(suite);
  check_RuntimeDim(suite);
  check_ProjectToManifold_Scalar(suite);
  check_ProjectToManifold_PerColumn(suite);
  check_RandomStructure(suite);
  check_RetractLocal_RoundTrip(suite);
  check_ProjectorRespectsGram(suite);
  check_FixedSizeConstruction(suite);
  check_ProjectToManifold_Idempotent(suite);
  check_LiftedDimsAndGram(suite);
  check_NegativeScalarInput(suite);
  check_PerColumn_PermutationSign(suite);
  check_Retract_StaysOnManifold(suite);
  check_TinyScales(suite);
  check_LocalZero_And_Additivity(suite);
  check_NoiseRobustProjection(suite);
  check_ManifoldDirectionalDerivative(suite);

  suite.summary();
  return (suite.passed == suite.total) ? 0 : 1;
}
