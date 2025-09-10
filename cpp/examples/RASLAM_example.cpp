//
// Created by jason on 1/30/25.
//
#include <iomanip>
#include <gtsam/base/timing.h>
#include <gtsam/slam/InitializePose.h>
#include <gtsam/slam/dataset.h>
#include <unordered_set>
#include "../utils.h"
#include "../RaFactor.h"
#include "../LiftedPose.h"
#include "../SEsyncFactor.h"
#include "../LandmarkFactor.h"
#include "../LiftedRangeFactor.h"
#include <gtsam/nonlinear/LevenbergMarquardtOptimizer.h>
#include <iostream>
#include <fstream>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <ctime>
#include <limits>
#include <filesystem>
#include "../json.hpp"
#include <chrono>

using nlohmann::json;
using clock_type = std::chrono::steady_clock; // monotonic, not affected by NTP
// ISO-8601 UTC timestamp
static std::string isoNowUTC() {
  std::time_t t = std::time(nullptr);
  std::tm tm{};
#if defined(_WIN32)
  gmtime_s(&tm, &t);
#else
  gmtime_r(&t, &tm);
#endif
  char buf[32];
  std::strftime(buf, sizeof(buf), "%Y-%m-%dT%H:%M:%SZ", &tm);
  return std::string(buf);
}

static std::string datasetKeyFromPath(const std::string& input_file) {
  std::filesystem::path p(input_file);
  return p.filename().string(); // e.g., "tiers.pyfg"
}
std::string datasetNameFromPath(const std::string& p) {
  std::filesystem::path x(p);
  auto stem = x.stem().string();   // "single_drone"
  return stem;
}


void appendRunToResultsJsonFlat(const std::string& out_path,
                                const std::string& dataset_name,        // e.g., "single_drone"
                                const string formulation,                         // e.g., 2
                                const std::string& init_file,            // full path to init file
                                const std::vector<double>& costs,        // cost at each iteration (len >= 1)
                                const std::vector<double>& times)        // cumulative time stamps (seconds), len == costs.size()
{
  using nlohmann::json;

  // 1) Load existing file (must be an array). If missing/invalid, start a new array.
  json root = json::array();
  {
    std::ifstream ifs(out_path);
    if (ifs.good()) {
      try {
        ifs >> root;
        if (!root.is_array()) {
          // If previous format was object, discard to match the new flat-array requirement.
          root = json::array();
        }
      } catch (...) {
        root = json::array();
      }
    }
  }

  // 2) Compose this run’s record
  json entry = {
    {"costs",        costs},
    {"dataset_name", dataset_name},
    {"formulation",  formulation},
    {"init_file",    init_file},
    {"times",        times}
  };

  // 3) Append and save (pretty)
  std::filesystem::create_directories(std::filesystem::path(out_path).parent_path());
  root.push_back(std::move(entry));
  std::ofstream ofs(out_path);
  ofs << std::setw(2) << root << '\n';
}

static inline void printVec(const std::string& name, const gtsam::Vector& v) {
  std::cout << name << " [len=" << v.size() << "]: [";
  std::cout.setf(std::ios::fixed); std::cout << std::setprecision(6);
  for (int i = 0; i < v.size(); ++i) {
    std::cout << v(i); if (i + 1 < v.size()) std::cout << ", ";
  }
  std::cout << "]\n";
}

static inline void printMat(const std::string& name, const gtsam::Matrix& M) {
  std::cout << name << " [" << M.rows() << "x" << M.cols() << "]:\n";
  std::cout.setf(std::ios::fixed); std::cout << std::setprecision(6);
  for (int r = 0; r < M.rows(); ++r) {
    std::cout << "  ";
    for (int c = 0; c < M.cols(); ++c) {
      std::cout << M(r,c);
      if (c + 1 < M.cols()) std::cout << "  ";
    }
    std::cout << "\n";
  }
}



static bool isDuplicateRunRecorded(const std::string& out_path,
                                   const std::string& dataset_name,
                                   const std::string& formulation,
                                   const std::string& init_file,
                                   bool match_on_basename = true) {
  using nlohmann::json;

  std::ifstream ifs(out_path);
  if (!ifs.good()) return false;  // no file yet → no duplicates

  json root;
  try {
    ifs >> root;
  } catch (...) {
    return false; // bad/partial JSON → treat as no duplicates
  }
  if (!root.is_array()) return false; // old object format → ignore

  // Normalize init_file key the same way your writer does
  const std::string init_key =
      match_on_basename ? std::filesystem::path(init_file).filename().string()
                        : init_file;

  for (const auto& e : root) {
    if (!e.is_object()) continue;

    const std::string e_ds  = e.value("dataset_name", "");
    const std::string e_for = e.value("formulation",  "");
    const std::string e_ini = e.value("init_file",    "");

    // Compare both as-is and by basename to be robust
    const std::string e_ini_base = std::filesystem::path(e_ini).filename().string();

    if (e_ds == dataset_name && e_for == formulation &&
        (e_ini == init_key || e_ini_base == init_key)) {
      return true; // found a duplicate
    }
  }
  return false;
}

using namespace std;
using namespace gtsam;

#include <sstream>
#include <unordered_map>
// Check 2x2 orthonormality
static inline bool isOrthonormal2(const gtsam::Matrix& R, double tol=1e-3) {
  if (R.rows()!=2 || R.cols()!=2) return false;
  gtsam::Matrix RtR = R.transpose()*R;
  return (RtR - gtsam::Matrix::Identity(2,2)).norm() < tol;
}


static inline std::string idKey(char tag, size_t idx) {
  return std::string(1, tag) + std::to_string(idx);
}
static inline std::string pairKey(char ti, size_t i, char tj, size_t j) {
  return idKey(ti, i) + "|" + idKey(tj, j);
}
static inline void parseId(const std::string& token, char& tag, size_t& idx) {
  if (token.empty()) throw std::runtime_error("Empty ID token");
  tag = token[0];
  idx = std::stoul(token.substr(1));
}

struct InitParsed {
  // Poses from file
  std::unordered_map<gtsam::Key, LiftedPoseDP> poses;
  // Landmarks/vectors from file
  std::unordered_map<gtsam::Key, gtsam::Vector> points;
  // Unit bearings keyed by (poseID | targetID) -> p-dim vector
  std::unordered_map<std::string, gtsam::Vector> bearings_pdim;
};

// ---------- helpers (no orthogonalization) ----------
static inline gtsam::Vector embedToP(const gtsam::Vector& v_d, int d, int p) {
  gtsam::Vector v_p = gtsam::Vector::Zero(p);
  v_p.head(d) = v_d;
  return v_p;
}

// Extract planar yaw from quaternion (x,y,z,w) for d=2 case
static inline double yaw_from_quat_xyzw(double qx,double qy,double qz,double qw) {
  double siny_cosp = 2.0 * (qw*qz + qx*qy);
  double cosy_cosp = 1.0 - 2.0 * (qy*qy + qz*qz);
  return std::atan2(siny_cosp, cosy_cosp);
}

static inline double orthoErrRtR(const gtsam::Matrix& R) {
  gtsam::Matrix I = gtsam::Matrix::Identity(R.cols(), R.cols());
  return (R.transpose() * R - I).norm();
}

// ---------- parser (no orthogonalization anywhere) ----------
static InitParsed parseInitFile(const std::string& path, int d, int p, bool VERBOSE_PARSE = true) {
  std::ifstream in(path);
  if (!in.is_open()) throw std::runtime_error("Could not open init file: " + path);

  InitParsed out;
  std::string line;

  if (VERBOSE_PARSE) {
    std::cout << "[InitFile] Reading initials from: " << path
              << "  (d=" << d << ", p=" << p << ")\n";
  }

  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') continue;

    std::istringstream iss(line);
    std::string tag; iss >> tag;

    // ---------------------------
    // VERTEX_POSE
    // ---------------------------
     if (tag == "VERTEX_POSE") {
      std::string id; iss >> id;
      char T; size_t I; parseId(id, T, I);

      std::vector<double> nums; nums.reserve(64);
      double val; while (iss >> val) nums.push_back(val);

      if (VERBOSE_PARSE) {
        std::cout << "\n[InitFile] VERTEX_POSE " << id << "\n";
        std::cout << "  raw nums (" << nums.size() << "): [";
        std::cout.setf(std::ios::fixed); std::cout << std::setprecision(6);
        for (size_t k = 0; k < nums.size(); ++k) {
          std::cout << nums[k]; if (k + 1 < nums.size()) std::cout << ", ";
        }
        std::cout << "]\n";
      }

      const size_t needLiftedVec = static_cast<size_t>(p + d * p);   // t_p + vec(d×p)
      const size_t needProblem   = static_cast<size_t>(d*d + d);     // R + t_d (order unknown)
      const size_t n = nums.size();

      gtsam::Vector t_p_(p);
      gtsam::Matrix Y_final;   // p×d
      bool built = false;

      // -------- 1) Lifted vectorized: t_p | vec(d×p) --------
      if (n == needLiftedVec) {
        // t_p
        for (int i = 0; i < p; ++i) t_p_(i) = nums[i];
        if (VERBOSE_PARSE) printVec("  t_p (from file)", t_p_);

        // vec(d×p)
        std::vector<double> rem(nums.begin() + p, nums.end()); // length d*p

        // Build Y candidates from ROW/COL major assumptions
        gtsam::Matrix Ydp_row(d, p), Ydp_col(d, p);
        { size_t k = 0; for (int r=0; r<d; ++r) for (int c=0; c<p; ++c,++k) Ydp_row(r,c) = rem[k]; }
        { size_t k = 0; for (int c=0; c<p; ++c) for (int r=0; r<d; ++r,++k) Ydp_col(r,c) = rem[k]; }
        gtsam::Matrix Y_row = Ydp_row.transpose(); // p×d
        gtsam::Matrix Y_col = Ydp_col.transpose(); // p×d

        const double e_row = orthoErrRtR(Y_row);
        const double e_col = orthoErrRtR(Y_col);
        const bool useRow = (e_row <= e_col);
        const gtsam::Matrix& Y_raw = useRow ? Y_row : Y_col;

        if (VERBOSE_PARSE) {
          printMat(useRow ? "  Y_raw (ROW-major vec)" : "  Y_raw (COL-major vec)", Y_raw);
          std::cout << "  orthoErr(YᵀY−I): row=" << e_row << " col=" << e_col
                    << "  -> chose " << (useRow ? "ROW" : "COL") << " (diagnostic only)\n";
          printMat("  Y_raw^T Y_raw (diagnostic)", Y_raw.transpose() * Y_raw);
        }

        try {
          StiefelManifoldKP Yst(Y_raw);       // will validate; we do NOT change Y
          LiftedPoseDP P(Yst, t_p_);
          if (VERBOSE_PARSE) { printMat("  Y(final, as-is)", Yst.matrix()); printVec("  t(final)", t_p_); }
          out.poses.emplace(gtsam::Symbol(T, I), P);
          built = true;
        } catch (const std::exception& e) {
          std::cerr << "[InitFile][WARNING] " << id << ": provided Y is not Stiefel (" << e.what()
                    << "). Falling back to identity R via Lift.\n";
          gtsam::Matrix R_id = gtsam::Matrix::Identity(d, d);
          StiefelManifoldKP Yst = StiefelManifoldKP::Lift(p, R_id);
          LiftedPoseDP P(Yst, t_p_);
          out.poses.emplace(gtsam::Symbol(T, I), P);
          built = true;
        }
      }

      // -------- 2) Problem-dim: R (d×d) + t_d (d) (order unknown) --------
      if (!built && n == needProblem) {
        // Candidate A: R-first
        gtsam::Matrix RA(d, d); { size_t k=0; for (int r=0;r<d;++r) for (int c=0;c<d;++c,++k) RA(r,c)=nums[k]; }
        gtsam::Vector tA(d); for (int i=0;i<d;++i) tA(i)=nums[d*d+i];

        // Candidate B: t-first
        gtsam::Vector tB(d); for (int i=0;i<d;++i) tB(i)=nums[i];
        gtsam::Matrix RB(d, d); { size_t k=d; for (int r=0;r<d;++r) for (int c=0;c<d;++c,++k) RB(r,c)=nums[k]; }

        const double eA = orthoErrRtR(RA);
        const double eB = orthoErrRtR(RB);
        const bool useA = (eA <= eB);
        const gtsam::Matrix& R_d = useA ? RA : RB;
        const gtsam::Vector& t_d = useA ? tA : tB;

        if (VERBOSE_PARSE) {
          std::cout << "  problem-dim pose detected (" << needProblem << " numbers)\n";
          printMat(useA ? "  R_d (R-first chosen)" : "  R_d (R-last chosen)", R_d);
          std::cout << "  orthoErr(RᵀR−I): A=" << eA << "  B=" << eB << "  -> chose " << (useA ?"A":"B") << "\n";
          printVec("  t_d (chosen)", t_d);
        }

        // Build Y via Lift(p, R_d); no modification to R_d is performed
        try {
          StiefelManifoldKP Yst = StiefelManifoldKP::Lift(p, R_d);
          t_p_.setZero(); t_p_.head(d) = t_d;
          LiftedPoseDP P(Yst, t_p_);
          if (VERBOSE_PARSE) { printMat("  Y(final) from Lift(p,R_d)", Yst.matrix()); printVec("  t_p (embedded)", t_p_); }
          out.poses.emplace(gtsam::Symbol(T, I), P);
          built = true;
        } catch (const std::exception& e) {
          std::cerr << "[InitFile][WARNING] " << id << ": Lift(p,R) failed (" << e.what()
                    << "). Falling back to identity R via Lift.\n";
          gtsam::Matrix R_id = gtsam::Matrix::Identity(d, d);
          StiefelManifoldKP Yst = StiefelManifoldKP::Lift(p, R_id);
          t_p_.setZero(); t_p_.head(d) = t_d;
          LiftedPoseDP P(Yst, t_p_);
          out.poses.emplace(gtsam::Symbol(T, I), P);
          built = true;
        }
      }

      // -------- 3) Legacy lifted matrix: row-major p×(d+1) [Y | t] --------
      const size_t needLiftedMat = static_cast<size_t>(p * (d + 1));
      if (!built && n == needLiftedMat) {
        gtsam::Matrix P_raw(p, d+1);
        { size_t k=0; for (int r=0;r<p;++r) for (int c=0;c<d+1;++c,++k) P_raw(r,c)=nums[k]; }
        if (VERBOSE_PARSE) printMat("  P_raw (p x (d+1))", P_raw);

        gtsam::Matrix Y_raw = P_raw.leftCols(d);
        t_p_ = P_raw.col(d);

        if (VERBOSE_PARSE) {
          printMat("  Y_raw (as-is)", Y_raw);
          printVec("  t_p (as-is)", t_p_);
          printMat("  Y_raw^T Y_raw (diagnostic)", Y_raw.transpose()*Y_raw);
        }

        try {
          StiefelManifoldKP Yst(Y_raw);     // validate; no modification
          LiftedPoseDP P(Yst, t_p_);
          out.poses.emplace(gtsam::Symbol(T, I), P);
          built = true;
        } catch (const std::exception& e) {
          std::cerr << "[InitFile][WARNING] " << id << ": Y not Stiefel (" << e.what()
                    << "). Falling back to identity R via Lift.\n";
          gtsam::Matrix R_id = gtsam::Matrix::Identity(d, d);
          StiefelManifoldKP Yst = StiefelManifoldKP::Lift(p, R_id);
          LiftedPoseDP P(Yst, t_p_);
          out.poses.emplace(gtsam::Symbol(T, I), P);
          built = true;
        }
      }

      if (!built) {
        throw std::runtime_error(
          "VERTEX_POSE " + id + ": unsupported number count (" + std::to_string(n) + "). "
          "Expected p+d*p (t_p|vec(d×p)), or d*d+d (R+t), or p*(d+1) (row-major [Y|t])."
        );
      }
    } else if (tag == "VERTEX_POINT") {
      std::string id; iss >> id;
      char T; size_t I; parseId(id, T, I);

      std::vector<double> nums; nums.reserve(std::max(d, p));
      double val; while (iss >> val) nums.push_back(val);

      if (VERBOSE_PARSE) {
        std::cout << "\n[InitFile] VERTEX_POINT " << id << "\n";
        std::cout << "  raw nums (" << nums.size() << "): [";
        std::cout.setf(std::ios::fixed); std::cout << std::setprecision(6);
        for (size_t k = 0; k < nums.size(); ++k) {
          std::cout << nums[k]; if (k + 1 < nums.size()) std::cout << ", ";
        }
        std::cout << "]\n";
      }

      gtsam::Vector v;
      if (nums.size() == static_cast<size_t>(p)) {
        v = Eigen::Map<Eigen::VectorXd>(nums.data(), p);                       // as-is
        if (VERBOSE_PARSE) printVec("  v_p (from file, as-is)", v);
      } else if (nums.size() == static_cast<size_t>(d)) {
        gtsam::Vector vd = Eigen::Map<Eigen::VectorXd>(nums.data(), d);
        v = embedToP(vd, d, p);                                                // zero-pad
        if (VERBOSE_PARSE) {
          printVec("  v_d (from file)", vd);
          printVec("  v_p (embedded)", v);
        }
      } else {
        throw std::runtime_error("VERTEX_POINT " + id + ": expected " + std::to_string(d) +
                                 " or " + std::to_string(p) + " numbers, got " + std::to_string(nums.size()));
      }

      out.points.emplace(gtsam::Symbol(T, I), v);

    // ---------------------------
    // VERTEX_BEARING
    // ---------------------------
    } else if (tag == "VERTEX_BEARING") {
      std::string id1, id2; iss >> id1 >> id2;
      char T1, T2; size_t I1, I2; parseId(id1, T1, I1); parseId(id2, T2, I2);

      std::vector<double> nums; nums.reserve(std::max(d, p));
      double val; while (iss >> val) nums.push_back(val);

      if (VERBOSE_PARSE) {
        std::cout << "\n[InitFile] VERTEX_BEARING " << id1 << " -> " << id2 << "\n";
        std::cout << "  raw nums (" << nums.size() << "): [";
        std::cout.setf(std::ios::fixed); std::cout << std::setprecision(6);
        for (size_t k = 0; k < nums.size(); ++k) {
          std::cout << nums[k]; if (k + 1 < nums.size()) std::cout << ", ";
        }
        std::cout << "]\n";
      }

      gtsam::Vector b_p;
      if (nums.size() == static_cast<size_t>(p)) {
        b_p = Eigen::Map<Eigen::VectorXd>(nums.data(), p);                     // as-is
        if (VERBOSE_PARSE) printVec("  b_p (from file, as-is)", b_p);
      } else if (nums.size() == static_cast<size_t>(d)) {
        gtsam::Vector b_d = Eigen::Map<Eigen::VectorXd>(nums.data(), d);
        b_p = embedToP(b_d, d, p);                                             // zero-pad
        if (VERBOSE_PARSE) {
          printVec("  b_d (from file)", b_d);
          printVec("  b_p (embedded)", b_p);
        }
      } else {
        throw std::runtime_error("VERTEX_BEARING " + id1 + "->" + id2 + ": expected " +
                                 std::to_string(d) + " or " + std::to_string(p) +
                                 " numbers, got " + std::to_string(nums.size()));
      }

      // Store (normalization for UnitSphereD happens later when you create R_k in main)
      out.bearings_pdim.emplace(pairKey(T1, I1, T2, I2), b_p);

    } else {
      // Unknown tag -> ignore silently
    }
  }

  if (VERBOSE_PARSE) {
    std::cout << "\n[InitFile] Done. Parsed: poses=" << out.poses.size()
              << ", points=" << out.points.size()
              << ", bearings=" << out.bearings_pdim.size() << "\n";
  }

  return out;
}



LiftedPoseDP LiftedToP(const Pose3 &pose3_, const size_t p) {
    return LiftedPoseDP(StiefelManifoldKP::Lift(p, pose3_.rotation().matrix()), pose3_.translation());
}

int main(int argc, char* argv[]) {
  using namespace std;
  using namespace gtsam;

  if (argc < 4) {
    cerr << "Usage: " << argv[0] << " <d> <p> <dataset.pyfg> [init.txt]\n";
    return 1;
  }

  // Args
  const int d = std::stoi(argv[1]);
  const int p = std::stoi(argv[2]);
  const std::string inputFile = argv[3];
  const std::string initFile  = (argc > 4) ? std::string(argv[4]) : std::string();

  if (!(d == 2 || d == 3)) {
    cerr << "[ERROR] d must be 2 or 3\n";
    return 2;
  }
  if (p < d) {
    cerr << "[ERROR] p must satisfy p >= d\n";
    return 3;
  }


  // Keys must match what you later write:
  const std::string dataset_name = datasetNameFromPath(inputFile);          // e.g., "single_drone"
  const std::string formulation  = "gtsam";                                  // your tag
  const std::string init_key     = (argc > 4) ? datasetKeyFromPath(initFile) // filename only
                                              : std::string();               // empty init is valid key

  // Build the same results path you write to later
  const std::string out_path = "/home/alan/varProj-gtsam/data/raslam/" + dataset_name + "/results.json";

  if (isDuplicateRunRecorded(out_path, dataset_name, formulation, init_key, /*match_on_basename=*/true)) {
    std::cout << "[SKIP] Duplicate run exists in " << out_path
              << " for (dataset=" << dataset_name
              << ", formulation=" << formulation
              << ", init=" << init_key << ").\n";
    return 0; // skip optimization
  }
  // Load dataset
  auto measurements = DataParser::read_pycfg_file(inputFile);
  cout << "Loaded " << measurements.poseMeasurements.size()     << " pose measurements. "
       << "Loaded " << measurements.landmarkMeasurements.size() << " landmark measurements. "
       << "Loaded " << measurements.rangeMeasurements.size()    << " Range measuerments. "
       << "Total Number of measurements: "
       << (measurements.landmarkMeasurements.size() + measurements.poseMeasurements.size() + measurements.rangeMeasurements.size())
       << ". Loaded " << measurements.num_poses << " poses from file. Loaded "
       << measurements.num_landmarks << " Landmarks. Total number of variables: "
       << (measurements.num_poses + measurements.num_landmarks + measurements.num_ranges) << std::endl;

  // Parse init file (poses/points/bearings already embedded to p in parser)
  InitParsed parsed;
  bool haveInit = false;
  if (!initFile.empty()) {
    try {
      parsed = parseInitFile(initFile, d, p,false); // your verbose parser
      haveInit = true;
    } catch (const std::exception& e) {
      cerr << "[InitFile][WARNING] " << e.what() << " — continuing with fallbacks.\n";
    }
  }

  // Build factor graph
  NonlinearFactorGraph inputGraph;

  // Pose–pose factors
  for (const auto& meas : measurements.poseMeasurements) {
    Vector sigmas = Vector::Zero(p * d + p);
    sigmas.head(p * d).setConstant(std::sqrt(1.0 / (1.0 * meas.kappa)));
    sigmas.tail(p).setConstant(std::sqrt(1.0 / (1.0 * meas.tau)));
    auto noise = noiseModel::Diagonal::Sigmas(sigmas);

    const char tag_i = meas.ci ? meas.ci : 'X';
    const char tag_j = meas.cj ? meas.cj : 'X';
    Key Ki = Symbol(tag_i, meas.i);
    Key Kj = Symbol(tag_j, meas.j);

    if (d == 2) {
      inputGraph.emplace_shared<SEsyncFactor2>(Ki, Kj, meas.R, meas.t, p, noise);
    } else {
      inputGraph.emplace_shared<SEsyncFactor3>(Ki, Kj, meas.R, meas.t, p, noise);
    }
  }

  // Pose–landmark factors
  for (const auto& meas : measurements.landmarkMeasurements) {
    Vector sigmas = Vector::Constant(p, std::sqrt(1.0 / meas.nu));
    auto noise = noiseModel::Diagonal::Sigmas(sigmas);

    const char tag_i = meas.ci ? meas.ci : 'X';
    const char tag_j = meas.cj ? meas.cj : 'L';
    Key Xi = Symbol(tag_i, meas.i);
    Key Lj = Symbol(tag_j, meas.j);

    if (d == 2) {
      inputGraph.emplace_shared<LiftedLandmarkFactor2>(Xi, Lj, meas.l, p, noise);
    } else {
      inputGraph.emplace_shared<LiftedLandmarkFactor3>(Xi, Lj, meas.l, p, noise);
    }
  }

  // --- Determine which variable keys MUST exist in Values ---
  std::set<std::pair<char, size_t>> poseKeys;      // (tag, idx)
  std::set<std::pair<char, size_t>> landmarkKeys;  // (tag, idx)
  for (const auto& m : measurements.poseMeasurements) {
    poseKeys.emplace(m.ci ? m.ci : 'X', m.i);
    poseKeys.emplace(m.cj ? m.cj : 'X', m.j);
  }
  for (const auto& m : measurements.landmarkMeasurements) {
    poseKeys.emplace(m.ci ? m.ci : 'X', m.i);
    landmarkKeys.emplace(m.cj ? m.cj : 'L', m.j);
  }
  for (const auto& m : measurements.rangeMeasurements) {
    poseKeys.emplace(m.ci ? m.ci : 'X', m.i);
    const char tj = m.cj ? m.cj : 'L';
    if (tj == 'L') landmarkKeys.emplace('L', m.j);
    else           poseKeys.emplace(tj, m.j);
  }

  // --- Build initial Values from init file; fallback if missing ---
  Values initial;

  // Poses
  for (const auto& [T, I] : poseKeys) {
    const Key K = Symbol(T, I);
    if (haveInit) {
      auto it = parsed.poses.find(K);
      if (it != parsed.poses.end()) {
        initial.insert(K, it->second);
      }
    }
  }

  // Landmarks / vectors
  for (const auto& [T, I] : landmarkKeys) {
    const Key K = Symbol(T, I);
    if (haveInit) {
      auto it = parsed.points.find(K);
      if (it != parsed.points.end()) {
        initial.insert(K, it->second);
      }
    }
  }

  // --- Range factors and *aligned* R_k initialization ---
  size_t k = 0;
  for (const auto& meas : measurements.rangeMeasurements) {
    Vector sigmas = Vector::Constant(p, std::sqrt(meas.sigma));
    auto noise = noiseModel::Diagonal::Sigmas(sigmas);

    const char tag_i = meas.ci ? meas.ci : 'X';
    const char tag_j = meas.cj ? meas.cj : 'L';

    const Key Pi = Symbol(tag_i, meas.i);
    const Key Rk = Symbol('R', k);

    // Initialize R_k *now*, using bearing from init file if available
    if (!initial.exists(Rk)) {
      // Look up bearing by (pose|target) key
      const std::string bkey = pairKey(tag_i, meas.i, tag_j, meas.j);
      Eigen::VectorXd unit(p);
      bool found = false;

      if (haveInit) {
        auto itB = parsed.bearings_pdim.find(bkey);
        if (itB != parsed.bearings_pdim.end()) {
          unit = itB->second; // already p-dim
        }
      }
      Matrix X0(p,1);
      X0.col(0) = unit;      // UnitSphereD expects a p×1 with unit-norm column
      UnitSphereD Y(X0);
      initial.insert(Rk, Y);
    }

    // Add the range factor that uses this exact R_k
    if (tag_j == 'L') {
      const Key Lj = Symbol('L', meas.j);
      if (d == 2) {
        inputGraph.emplace_shared<LiftedRangeFactor2>(Pi, Lj, Rk, meas.range, p, noise);
      } else {
        inputGraph.emplace_shared<LiftedRangeFactor3>(Pi, Lj, Rk, meas.range, p, noise);
      }
    } else {
      const Key Pj = Symbol(tag_j, meas.j);
      if (d == 2) {
        inputGraph.emplace_shared<LiftedRangeFactorPosePose<2>>(Pi, Pj, Rk, meas.range, p, noise);
      } else {
        inputGraph.emplace_shared<LiftedRangeFactorPosePose<3>>(Pi, Pj, Rk, meas.range, p, noise);
      }
    }

    ++k; // IMPORTANT: increments in lock-step with factor creation
  }

  cout << "Total Measurments from graph " << inputGraph.size() << ".\n";

  // --- Diagnostics ---
  auto printKey = [](gtsam::Key K) {
    gtsam::Symbol s(K);
    std::cout << s.chr() << s.index();
  };

  const size_t numRangeAux = k;

  cout << "\n[Init] ================= INITIALIZATION SUMMARY =================\n";
  cout << "[Init] d=" << d << "  p=" << p << "\n";

  cout << "[Init] Pose keys (" << poseKeys.size() << "): ";
  for (const auto& [T, I] : poseKeys) {
    const Key K = Symbol(T, I);
    cout << (initial.exists(K) ? "" : "[MISSING!] ");
    printKey(K); cout << " ";
  }
  cout << "\n";

  cout << "[Init] Landmark keys (" << landmarkKeys.size() << "): ";
  for (const auto& [T, I] : landmarkKeys) {
    const Key K = Symbol(T, I);
    cout << (initial.exists(K) ? "" : "[MISSING!] ");
    printKey(K); cout << " ";
  }
  cout << "\n";

  cout << "[Init] Range auxiliaries R_k (" << numRangeAux << "): ";
  for (size_t kk = 0; kk < numRangeAux; ++kk) {
    const Key K = Symbol('R', kk);
    cout << (initial.exists(K) ? "" : "[MISSING!] ");
    printKey(K); cout << " ";
  }
  cout << "\n";

  const size_t expected = poseKeys.size() + landmarkKeys.size() + numRangeAux;
  cout << "[Init] Values.size() = " << initial.size()
       << "  (expected ~ " << expected << ")\n";
  if (initial.size() != expected) {
    cerr << "[Init][WARNING] Values size mismatch; check duplicate keys or insertion logic.\n";
  }
  cout << "[Init] ===========================================================\n\n";

  // --- Optimize ---
  auto lmParams = LevenbergMarquardtParams::CeresDefaults();
  lmParams.maxIterations    = 1000;
  lmParams.relativeErrorTol = 1e-20;
  lmParams.verbosityLM      = LevenbergMarquardtParams::SUMMARY;

  auto lm = std::make_shared<LevenbergMarquardtOptimizer>(inputGraph, initial, lmParams);
  auto t2 = clock_type::now();

  std::vector<double> costs;
  std::vector<double> times;  // cumulative seconds
  costs.reserve(lmParams.maxIterations + 1);
  times.reserve(lmParams.maxIterations + 1);
  costs.push_back(lm->error());
  times.push_back(0.0);
  double prevErr = std::numeric_limits<double>::infinity();
  for (size_t it = 0; it < lmParams.maxIterations; ++it) {
    auto t0 = clock_type::now();

    // one step of LM
    // (In most GTSAM builds this is public; if your version returns bool, use it to detect convergence.)
    lm->iterate();

    auto t1 = clock_type::now();
    times.push_back( std::chrono::duration<double>(t1 - t0).count());

    // your own convergence check (mirrors NonlinearOptimizer)
    double err = lm->error();
    costs.push_back(err);
    auto t3 = clock_type::now();
    if ((std::abs(prevErr - err) < lmParams.relativeErrorTol * err) || std::chrono::duration<double>(t3 - t2).count() > 20) break;
    prevErr = err;
  }



  // Example metadata for the JSON writer:
  appendRunToResultsJsonFlat(
    "/home/alan/varProj-gtsam/data/raslam/" + dataset_name + "/results.json",
    dataset_name, formulation, init_key, costs, times);

  cout << "[Done] iterations=" << lm->iterations()
       << "  final_error=" << lm->error()<<"\n";
  return 0;
}
