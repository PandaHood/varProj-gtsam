//
// Created by nikolas on 9/3/25.
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
#include <../json.hpp>
#include <chrono>
#include <random>

using namespace std;
using namespace gtsam;

// Monotonic clock alias
using clock_type = std::chrono::steady_clock;

static void printKeyCompact(gtsam::Key K) {
  gtsam::Symbol s(K);
  std::cout << s.chr() << s.index();
}

// -----------------------------
// Helpers for init file parsing
// -----------------------------
static inline void parseId(const std::string& token, char& tag, size_t& idx) {
  if (token.empty()) throw std::runtime_error("Empty ID token");
  tag = token[0];
  idx = std::stoul(token.substr(1));
}
static inline std::string idKey(char tag, size_t idx) {
  return std::string(1, tag) + std::to_string(idx);
}
static inline std::string pairKey(char ti, size_t i, char tj, size_t j) {
  return idKey(ti, i) + "|" + idKey(tj, j);
}
static inline gtsam::Vector embedToP(const gtsam::Vector& v_d, int d, int p) {
  gtsam::Vector v_p = gtsam::Vector::Zero(p);
  v_p.head(d) = v_d;
  return v_p;
}
static std::string datasetKeyFromPath(const std::string& input_file) {
  std::filesystem::path p(input_file);
  return p.filename().string(); // e.g., "tiers.pyfg"
}
static inline gtsam::Vector takeLastPAsVector(const std::vector<double>& nums, int p) {
  if ((int)nums.size() < p) throw std::runtime_error("Not enough numbers to take last p elements");
  gtsam::Vector v(p);
  const int start = static_cast<int>(nums.size()) - p;
  for (int i = 0; i < p; ++i) v(i) = nums[start + i];
  return v;
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

// Add these two includes near your other includes:
#include <sstream>
#include <unordered_map>
// If you want the optional Stiefel projection helper:
#include <Eigen/Dense>
#include <Eigen/SVD>

// -----------------------------
// BA init structure (numeric only)
// -----------------------------
struct BAInit {
  // Pose A<i>/B<i>/...: translation t_p (R^p) and Y (p×d)
  std::unordered_map<gtsam::Key, gtsam::Vector> t_p;
  std::unordered_map<gtsam::Key, gtsam::Matrix> Y_pd;
  // Landmarks L<j>: p-dim vectors
  std::unordered_map<gtsam::Key, gtsam::Vector> points;
};

// -----------------------------
// Small helpers (mirrors SNL parser style)
// -----------------------------
static inline void baParseId(const std::string& token, char& tag, size_t& idx) {
  if (token.empty()) throw std::runtime_error("Empty ID token");
  tag = token[0];
  idx = std::stoul(token.substr(1));
}
static inline gtsam::Vector baEmbedToP(const gtsam::Vector& v_d, int d, int p) {
  gtsam::Vector v_p = gtsam::Vector::Zero(p);
  v_p.head(d) = v_d;
  return v_p;
}
static inline gtsam::Vector baTakeLastPAsVector(const std::vector<double>& nums, int p) {
  if ((int)nums.size() < p) throw std::runtime_error("Not enough numbers to take last p elements");
  gtsam::Vector v(p);
  const int start = static_cast<int>(nums.size()) - p;
  for (int i = 0; i < p; ++i) v(i) = nums[start + i];
  return v;
}

// Optional: nearest Stiefel (p×d) via SVD: M ≈ Q = U(:,1:d) * Vᵀ
static inline gtsam::Matrix baProjectToStiefel(const gtsam::Matrix& M) {
  const int p = (int)M.rows(), d = (int)M.cols();
  Eigen::JacobiSVD<gtsam::Matrix> svd(M, Eigen::ComputeFullU | Eigen::ComputeFullV);
  gtsam::Matrix U = svd.matrixU();              // p×p
  gtsam::Matrix V = svd.matrixV();              // d×d
  return U.leftCols(d) * V.transpose();         // p×d, with orthonormal columns
}

// -----------------------------
// parseBAInitFile
// -----------------------------
// Accepts lines like:
//   VERTEX_POINT L5  <p or d or >p numbers>
//   VERTEX_POSE  A81 <p (or d) numbers for t> <p*d numbers for Y>
// Robustness:
//  - If more than needed numbers are present, takes the last (p + p*d) for poses,
//    and the last p for points.
//  - If translation t is provided with length d, it's zero-padded to p.
//  - Y is read row-major: the stream order is [r0c0 r0c1 ... r0c(d-1) r1c0 ... r(p-1)c(d-1)].
static BAInit parseBAInitFile(const std::string& path, int d, int p,
                              bool orthonormalizeY = false, bool VERBOSE = false) {
  std::ifstream in(path);
  if (!in.is_open()) throw std::runtime_error("Could not open BA init file: " + path);

  if (VERBOSE) {
    std::cout << "[BA Init] Reading: " << path << "  (d=" << d << ", p=" << p << ")\n";
  }

  BAInit out;
  std::string line;

  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') continue;

    std::istringstream iss(line);
    std::string tag; iss >> tag;

    if (tag == "VERTEX_POINT") {
      std::string id; iss >> id;
      char T; size_t I; baParseId(id, T, I);
      if (T != 'L') continue;  // only landmarks here

      std::vector<double> nums;
      double v; while (iss >> v) nums.push_back(v);

      gtsam::Vector vp;
      if (nums.size() == (size_t)p) {
        vp = Eigen::Map<Eigen::VectorXd>(nums.data(), p);
      } else if (nums.size() == (size_t)d) {
        gtsam::Vector vd = Eigen::Map<Eigen::VectorXd>(nums.data(), d);
        vp = baEmbedToP(vd, d, p);
      } else if (nums.size() > (size_t)p) {
        vp = baTakeLastPAsVector(nums, p);
      } else {
        throw std::runtime_error("VERTEX_POINT " + id + ": expected d or p numbers (or trailing format), got " + std::to_string(nums.size()));
      }
      out.points.emplace(gtsam::Symbol(T, I), std::move(vp));
      continue;
    }

    if (tag == "VERTEX_POSE") {
      std::string id; iss >> id;
      char T; size_t I; baParseId(id, T, I);

      std::vector<double> nums;
      double x; while (iss >> x) nums.push_back(x);

      const size_t need = (size_t)p + (size_t)(p * d);
      const size_t needAlt = (size_t)d + (size_t)(p * d);

      if (nums.size() < need && nums.size() != needAlt) {
        throw std::runtime_error(
          "VERTEX_POSE " + id + ": need p + p*d or d + p*d numbers, got " + std::to_string(nums.size()));
      }

      // Decide the slice we use:
      size_t start = 0;           // where the (t, Y) block starts in nums
      size_t tCount = p;          // how many entries to take for t
      bool tIsD = false;

      if (nums.size() == need) {      // exactly p + p*d
        start = 0; tCount = p;
      } else if (nums.size() == needAlt) { // exactly d + p*d
        start = 0; tCount = d; tIsD = true;
      } else if (nums.size() > need) {
        // Be liberal: take the last (p + p*d) numbers
        start = nums.size() - need; tCount = p;
      } else {
        // (nums.size() < need) but equals needAlt was handled above
        start = 0; tCount = d; tIsD = true;
      }

      // Extract translation t
      gtsam::Vector t_vec;
      if (!tIsD) {
        t_vec = gtsam::Vector::Map(nums.data() + start, p);
      } else {
        gtsam::Vector t_d = gtsam::Vector::Map(nums.data() + start, d);
        t_vec = baEmbedToP(t_d, d, p);
      }

      // Extract Y (row-major) immediately after t
      const size_t yOffset = start + tCount;
      if (yOffset + (size_t)(p * d) > nums.size()) {
        throw std::runtime_error("VERTEX_POSE " + id + ": not enough numbers for Y (p*d)");
      }

      gtsam::Matrix Y(p, d);
      // Row-major fill:
      // nums[yOffset + r*d + c] -> Y(r,c)
      for (int r = 0; r < p; ++r) {
        for (int c = 0; c < d; ++c) {
          Y(r, c) = nums[yOffset + (size_t)r * (size_t)d + (size_t)c];
        }
      }

      if (orthonormalizeY) {
        Y = baProjectToStiefel(Y);  // optional: nearest Stiefel
      }

      const gtsam::Key K = gtsam::Symbol(T, I);
      out.t_p.emplace(K, std::move(t_vec));
      out.Y_pd.emplace(K, std::move(Y));
      continue;
    }

    // Unknown line type: ignore
  }

  if (VERBOSE) {
    std::cout << "[BA Init] Parsed: poses=" << out.t_p.size()
              << ", landmarks=" << out.points.size() << "\n";
  }
  return out;
}

int main(int argc, char* argv[])
{
  if (argc < 4) {
    std::cerr << "Usage: " << argv[0] << " <d> <p> <dataset.pyfg>\n";
    return 1;
  }

  // Args
  const int d = std::stoi(argv[1]);      // spatial dimension (2 or 3)
  const int p = std::stoi(argv[2]);      // lifted dimension (>= d)
  const std::string inputFile = argv[3];

  std::string initFile = (argc > 4) ? std::string(argv[4]) : std::string();


  if (!(d == 2 || d == 3)) {
    std::cerr << "[ERROR] d must be 2 or 3\n";
    return 2;
  }
  if (p < d) {
    std::cerr << "[ERROR] p must satisfy p >= d (got p=" << p << ", d=" << d << ")\n";
    return 3;
  }

  // Load dataset (expects DataParser::read_pycfg_file in your project)
  auto measurements = DataParser::read_pycfg_file(inputFile);
  std::cout << "Loaded " << measurements.poseMeasurements.size()     << " pose measurements. "
            << "Loaded " << measurements.landmarkMeasurements.size() << " landmark measurements. "
            << "Loaded " << measurements.rangeMeasurements.size()    << " Range measuerments. "
            << "Total Number of measurements: "
            << (measurements.landmarkMeasurements.size() + measurements.poseMeasurements.size() + measurements.rangeMeasurements.size())
            << ". Loaded " << measurements.num_poses << " poses from file. Loaded "
            << measurements.num_landmarks << " Landmarks. Total number of variables: "
            << (measurements.num_poses + measurements.num_landmarks + measurements.num_ranges) << std::endl;

  BAInit baInit;
  bool haveBAInit = false;
  if (!initFile.empty()) {
    try {
      baInit = parseBAInitFile(initFile, d, p, /*orthonormalizeY=*/false, /*VERBOSE=*/false);
      haveBAInit = true;
    } catch (const std::exception& e) {
      std::cerr << "[BA Init][WARNING] " << e.what() << " — continuing with random fallbacks.\n";
    }
  }
  // -----------------------------
  // Build a BA graph: ONLY LiftedLandmark factors (pose-landmark)
  // -----------------------------
  NonlinearFactorGraph graph;

  size_t landmarkFactorCount = 0;
  for (const auto& meas : measurements.landmarkMeasurements) {
    // Residual ~ p (translation-like in lifted space)
    // If nu is a precision, sigma = sqrt(1/nu); adjust to your convention
    Vector sigmas = Vector::Constant(p, std::sqrt(1.0 / std::max(1e-12, meas.nu)));
    auto noise = noiseModel::Diagonal::Sigmas(sigmas);

    // Respect tags if provided; fall back to 'X' for poses and 'L' for landmarks
    const char tag_i = meas.ci ? meas.ci : 'X';  // pose
    const char tag_j = meas.cj ? meas.cj : 'L';  // landmark
    Key Xi = Symbol(tag_i, meas.i);
    Key Lj = Symbol(tag_j, meas.j);

    if (d == 2) {
      graph.emplace_shared<LiftedLandmarkFactor2>(Xi, Lj, meas.l, p, noise);
    } else {
      graph.emplace_shared<LiftedLandmarkFactor3>(Xi, Lj, meas.l, p, noise);
    }
    ++landmarkFactorCount;
  }

  std::cout << "[Graph] Built with " << landmarkFactorCount
            << " LiftedLandmark factors (BA-only).\n";

  // -----------------------------
  // Random initialization (poses + landmarks used in the BA graph)
  // -----------------------------
  // Collect exactly the keys that appear in the landmark factors we added
  std::set<std::pair<char, size_t>> poseKeys;      // (tag, idx)
  std::set<std::pair<char, size_t>> landmarkKeys;  // (tag, idx)

  for (const auto& m : measurements.landmarkMeasurements) {
    poseKeys.emplace(m.ci ? m.ci : 'X', m.i);
    landmarkKeys.emplace(m.cj ? m.cj : 'L', m.j);
  }

  // 2) Initialization
  Values initial;

  // Poses
  for (const auto& [tag, idx] : poseKeys) {
    const gtsam::Key K = Symbol(tag, idx);
    Eigen::VectorXd t_p(p); t_p.setRandom();
    if (haveBAInit) {
      if (auto it = baInit.t_p.find(K); it != baInit.t_p.end()) t_p = it->second;
    }

    StiefelManifoldKP Yst = StiefelManifoldKP::Random(std::default_random_engine::default_seed, d, p);
    if (haveBAInit) {
      if (auto itM = baInit.Y_pd.find(K); itM != baInit.Y_pd.end()) {
        gtsam::Matrix Q = baProjectToStiefel(itM->second);  // ensures on Stiefel
        // ⬇️ Replace with your actual constructor-from-matrix, if available:
        Yst = StiefelManifoldKP(Q);
        // e.g., Yst = StiefelManifoldKP::FromMatrix(Q);   // if you have such an API
      }
    }

    initial.insert(K, LiftedPoseDP(Yst, t_p));
  }

  // Landmarks
  for (const auto& [tag, idx] : landmarkKeys) {
    const gtsam::Key K = Symbol(tag, idx);
    Eigen::VectorXd v(p); v.setRandom();
    if (haveBAInit) {
      if (auto it = baInit.points.find(K); it != baInit.points.end()) v = it->second;
    }
    initial.insert(K, v);
  }

  // Diagnostics on initialization
  std::cout << "\n[Init] ================== RANDOM INITIALIZATION SUMMARY ==================\n";
  std::cout << "[Init] d=" << d << "  p=" << p << "\n";
  std::cout << "[Init] Pose keys (" << poseKeys.size() << "): ";
  for (const auto& [tag, idx] : poseKeys) {
    const Key K = Symbol(tag, idx);
    std::cout << (initial.exists(K) ? "" : "[MISSING!] ");
    printKeyCompact(K); std::cout << " ";
  }
  std::cout << "\n[Init] Landmark keys (" << landmarkKeys.size() << "): ";
  for (const auto& [tag, idx] : landmarkKeys) {
    const Key K = Symbol(tag, idx);
    std::cout << (initial.exists(K) ? "" : "[MISSING!] ");
    printKeyCompact(K); std::cout << " ";
  }
  std::cout << "\n[Init] Values.size() = " << initial.size()
            << "  (expected ~ " << (poseKeys.size() + landmarkKeys.size()) << ")\n";
  std::cout << "[Init] ====================================================================\n\n";

  if (graph.empty()) {
    std::cerr << "[ERROR] No LiftedLandmark factors found in dataset; nothing to optimize.\n";
    return 4;
  }

  // --- Optimize ---
  auto lmParams = LevenbergMarquardtParams::CeresDefaults();
  lmParams.maxIterations    = 1000;
  lmParams.relativeErrorTol = 1e-20;
  lmParams.verbosityLM      = LevenbergMarquardtParams::SUMMARY;

  auto lm = std::make_shared<LevenbergMarquardtOptimizer>(graph, initial, lmParams);
  auto t2 = clock_type::now();

  std::vector<double> costs;
  std::vector<double> times;  // cumulative seconds
  costs.reserve(lmParams.maxIterations + 1);
  times.reserve(lmParams.maxIterations + 1);
  costs.push_back(lm->error());
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
  std::string dataset_name = datasetNameFromPath(std::string(argv[3]));
  string formulation = "gtsam"; // your code’s tag for “gtsam” / formulation type
  std::string init_file   = datasetKeyFromPath(std::string(argv[4]));

  appendRunToResultsJsonFlat(
    "/home/nikolas/varProj-gtsam/data/sfm/" + dataset_name + "/results.json",
    dataset_name, formulation, init_file, costs, times);

  cout << "[Done] iterations=" << lm->iterations()
       << "  final_error=" << lm->error()<<"\n";
  return 0;
}



