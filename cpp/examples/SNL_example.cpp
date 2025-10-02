
//
// Created by nikolas on 9/3/25.
//

#include <iomanip>
#include <gtsam/base/timing.h>
#include <gtsam/slam/InitializePose.h>
#include <gtsam/slam/dataset.h>
#include <unordered_set>
#include "../utils.h"
#include "../LiftedPose.h"
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
#include <random>
#include <sstream>
#include <unordered_map>

using namespace std;
using namespace gtsam;

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

struct SNLInit {
  // Landmarks L<i> initialized to p-dim vectors
  std::unordered_map<gtsam::Key, gtsam::Vector> points;
  // Unit bearings keyed "Li|Lj" -> p-dim unit vector
  std::unordered_map<std::string, gtsam::Vector> bearings_pdim;
};

// Parse a lightweight init file for SNL:
//   VERTEX_POINT L<i> [p or d numbers]
//   VERTEX_BEARING L<i> L<j> [p or d numbers]   (if more than p numbers, takes the last p)
// Notes:
//  - Any d-dim entries are zero-padded to p.
//  - Bearings are normalized to unit length.
static SNLInit parseSNLInitFile(const std::string& path, int d, int p, bool VERBOSE) {
  std::ifstream in(path);
  if (!in.is_open()) throw std::runtime_error("Could not open init file: " + path);

  if (VERBOSE) {
    std::cout << "[InitFile] Reading SNL initials from: " << path
              << "  (d=" << d << ", p=" << p << ")\n";
  }

  SNLInit out;
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream iss(line);
    std::string tag; iss >> tag;

    if (tag == "VERTEX_POINT") {
      std::string id; iss >> id;
      char T; size_t I; parseId(id, T, I);
      if (T != 'L') continue; // SNL landmarks are 'L'

      std::vector<double> nums;
      double v; while (iss >> v) nums.push_back(v);

      gtsam::Vector vp;
      if (nums.size() == static_cast<size_t>(p)) {
        vp = Eigen::Map<Eigen::VectorXd>(nums.data(), p);
      } else if (nums.size() == static_cast<size_t>(d)) {
        gtsam::Vector vd = Eigen::Map<Eigen::VectorXd>(nums.data(), d);
        vp = embedToP(vd, d, p);
      } else if (nums.size() > static_cast<size_t>(p)) {
        // Allow extra trailing numbers; take the last p as the point (robust to weird formats)
        vp = takeLastPAsVector(nums, p);
      } else {
        throw std::runtime_error("VERTEX_POINT " + id + ": expected d or p numbers (or trailing format), got " + std::to_string(nums.size()));
      }
      out.points.emplace(gtsam::Symbol(T, I), vp);

    } else if (tag == "VERTEX_BEARING") {
      std::string id1, id2; iss >> id1 >> id2;
      char T1, T2; size_t I1, I2; parseId(id1, T1, I1); parseId(id2, T2, I2);
      if (T1 != 'L' || T2 != 'L') continue;

      std::vector<double> nums;
      double v; while (iss >> v) nums.push_back(v);

      gtsam::Vector bp;
      if (nums.size() == static_cast<size_t>(p)) {
        bp = Eigen::Map<Eigen::VectorXd>(nums.data(), p);
      } else if (nums.size() == static_cast<size_t>(d)) {
        gtsam::Vector bd = Eigen::Map<Eigen::VectorXd>(nums.data(), d);
        bp = embedToP(bd, d, p);
      } else if (nums.size() > static_cast<size_t>(p)) {
        // Be liberal: take the last p numbers as the bearing vector.
        bp = takeLastPAsVector(nums, p);
      } else {
        throw std::runtime_error("VERTEX_BEARING " + id1 + " " + id2 + ": expected d or p numbers (or trailing format), got " + std::to_string(nums.size()));
      }

      const double nrm = bp.norm();
      if (nrm > 1e-12) bp /= nrm;
      out.bearings_pdim.emplace(pairKey(T1, I1, T2, I2), bp);
    }
  }

  if (VERBOSE) {
    std::cout << "[InitFile] Parsed: points=" << out.points.size()
              << ", bearings=" << out.bearings_pdim.size() << "\n";
  }
  return out;
}

// Monotonic clock alias
using clock_type = std::chrono::steady_clock;

static void printKeyCompact(gtsam::Key K) {
  gtsam::Symbol s(K);
  std::cout << s.chr() << s.index();
}

int main(int argc, char* argv[])
{
  if (argc < 4) {
    std::cerr << "Usage: " << argv[0] << " <d> <p> <dataset.pyfg> [init.txt]\n";
    return 1;
  }

  // Args
  const int d = std::stoi(argv[1]);      // spatial dimension (2 or 3)
  const int p = std::stoi(argv[2]);      // lifted dimension (>= d)
  const std::string inputFile = argv[3];
  const std::string initFile  = (argc > 4) ? std::string(argv[4]) : std::string();

  if (!(d == 2 || d == 3)) {
    std::cerr << "[ERROR] d must be 2 or 3\n";
    return 2;
  }
  if (p < d || p > 7) { // per your spec: up to rank 7
    std::cerr << "[ERROR] p must satisfy d <= p <= 7 (got p=" << p << ", d=" << d << ")\n";
    return 3;
  }

  // Keys must match what you later write:
  const std::string dataset_name = datasetNameFromPath(inputFile);          // e.g., "single_drone"
  const std::string formulation  = "gtsam";                                  // your tag
  const std::string init_key     = (argc > 4) ? datasetKeyFromPath(initFile) // filename only
                                              : std::string();               // empty init is valid key

  // Build the same results path you write to later
  const std::string out_path = "/home/alan/varProj-gtsam/data/snl/" + dataset_name + "/results.json";

  if (isDuplicateRunRecorded(out_path, dataset_name, formulation, init_key, /*match_on_basename=*/true)) {
    std::cout << "[SKIP] Duplicate run exists in " << out_path
              << " for (dataset=" << dataset_name
              << ", formulation=" << formulation
              << ", init=" << init_key << ").\n";
    return 0; // skip optimization
  }

  // Load dataset
  auto measurements = DataParser::read_pycfg_file(inputFile);
  std::cout << "Loaded " << measurements.poseMeasurements.size()     << " pose measurements. "
            << "Loaded " << measurements.landmarkMeasurements.size() << " landmark measurements. "
            << "Loaded " << measurements.rangeMeasurements.size()    << " Range measuerments. "
            << "Total Number of measurements: "
            << (measurements.landmarkMeasurements.size() + measurements.poseMeasurements.size() + measurements.rangeMeasurements.size())
            << ". Loaded " << measurements.num_poses << " poses from file. Loaded "
            << measurements.num_landmarks << " Landmarks. Total number of variables: "
            << (measurements.num_poses + measurements.num_landmarks + measurements.num_ranges) << std::endl;

  // Optional init file
  SNLInit parsed;
  bool haveInit = false;
  if (!initFile.empty()) {
    try {
      parsed = parseSNLInitFile(initFile, d, p, /*VERBOSE=*/false);
      haveInit = true;
    } catch (const std::exception& e) {
      std::cerr << "[InitFile][WARNING] " << e.what() << " — continuing with random fallbacks.\n";
    }
  }

  // -----------------------------
  // Build a pure SNL graph
  //   - variables:   L<i> (p-dim vectors), R<k> (UnitSphereD in R^p)
  //   - factors:     LiftedRangeFactorPointPoint<d>(L_i, L_j, R_k, range, p, noise)
  //   - init:        from file when present; otherwise random
  // -----------------------------
  NonlinearFactorGraph inputGraph;
  Values initial;

  // 1) Collect exactly the landmark ids that appear in the ranges
  std::set<size_t> landmarkIds;
  for (const auto& m : measurements.rangeMeasurements) {
    landmarkIds.insert(m.i);
    landmarkIds.insert(m.j);
  }

  // 2) Initialize those landmarks (use init file when possible, else random)
  for (size_t lid : landmarkIds) {
    const Key K = Symbol('L', lid);
    bool set = false;
    if (haveInit) {
      auto it = parsed.points.find(K);
      if (it != parsed.points.end()) {
        initial.insert(K, it->second); // p-dim already (or embedded)
        set = true;
      }
    }
    if (!set) {
      gtsam::Vector x = gtsam::Vector::Random(p);  // uniform in [-1,1]
      initial.insert(K, x);
    }
  }

 // 3) Add one UnitSphereD R_k per edge and the LiftedRangePointPoint factor
  size_t k = 0;
  for (const auto& m : measurements.rangeMeasurements) {
    // noise: if sigma is variance, use sqrt; if already std-dev, use as-is
    gtsam::Vector sigmas = gtsam::Vector::Constant(p, std::sqrt(std::max(1e-12, m.sigma)));
    auto noise = gtsam::noiseModel::Diagonal::Sigmas(sigmas);

    Key Li_meas = Symbol('L', m.i);
    Key Lj_meas = Symbol('L', m.j);
    Key Rk      = Symbol('R', k++);

    // Decide orientation based on available bearing lines (do NOT reverse vectors)
    bool haveBearing = false;
    bool swapped = false;          // if true, we'll use (Lj, Li) order
    gtsam::Vector dir;

    if (haveInit) {
      const std::string fwd = pairKey('L', m.i, 'L', m.j); // Li -> Lj
      const std::string rev = pairKey('L', m.j, 'L', m.i); // Lj -> Li
      auto itF = parsed.bearings_pdim.find(fwd);
      if (itF != parsed.bearings_pdim.end()) {
        dir = itF->second;         // already unit
        haveBearing = (dir.norm() > 1e-12);
        swapped = false;           // keep (Li, Lj)
      } else {
        auto itR = parsed.bearings_pdim.find(rev);
        if (itR != parsed.bearings_pdim.end()) {
          dir = itR->second;       // already unit; keep as-is (no flipping)
          haveBearing = (dir.norm() > 1e-12);
          swapped = true;          // use (Lj, Li)
        }
      }
    }

    // Choose final variable order consistent with the selected bearing (if any)
    Key A = Li_meas, B = Lj_meas;
    if (swapped) std::swap(A, B);  // respect Lj->Li if that's what we found

    if (!haveBearing) {
      // Try using point difference in the chosen (A,B) order
      if (initial.exists(A) && initial.exists(B)) {
        const auto &vA = initial.at<gtsam::Vector>(A);
        const auto &vB = initial.at<gtsam::Vector>(B);
        dir = vB - vA;             // direction from A to B
        if (dir.norm() > 1e-12) {
          dir.normalize();
          haveBearing = true;
        }
      }
    }

    if (!haveBearing) {
      // Random unit vector in R^p
      dir = gtsam::Vector::Random(p);
      if (dir.norm() < 1e-12) dir.setConstant(1.0);
      dir.normalize();
    }

    // Insert R_k (UnitSphereD over R^p). Adjust ctor if your UnitSphereD differs.
    UnitSphereD Y(dir);
    initial.insert(Rk, Y);

    // Add the factor in the chosen (A, B) order
    if (d == 2) {
      inputGraph.emplace_shared<LiftedRangeFactorPointPoint<2>>(A, B, Rk, m.range, p, noise);
    } else { // d == 3
      inputGraph.emplace_shared<LiftedRangeFactorPointPoint<3>>(A, B, Rk, m.range, p, noise);
    }
  }

  if (inputGraph.empty()) {
    std::cerr << "[ERROR] No range factors found; nothing to optimize.\n";
    return 4;
  }

  std::cout << "Total range factors: " << inputGraph.size()
            << " | variables init: " << initial.size() << "\n";

  // Diagnostics (optional)
  // {
  //   std::cout << "[Init] Landmarks initialized: ";
  //   for (auto lid : landmarkIds) { printKeyCompact(Symbol('L', lid)); std::cout << " "; }
  //   std::cout << "\n[Init] Range auxiliaries R_k: ";
  //   for (size_t kk = 0; kk < k; ++kk) { printKeyCompact(Symbol('R', kk)); std::cout << " "; }
  //   std::cout << "\n";
  // }

  // --- Optimize ---
  auto lmParams = LevenbergMarquardtParams::CeresDefaults();
  lmParams.maxIterations    = 1000;
  lmParams.relativeErrorTol = 1e-20;
  lmParams.verbosityLM      = LevenbergMarquardtParams::SILENT;

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
    "/home/alan/varProj-gtsam/data/snl/" + dataset_name + "/results.json",
    dataset_name, formulation, init_key, costs, times);

  cout << "[Done] iterations=" << lm->iterations()
       << "  final_error=" << lm->error()<<"\n";
  return 0;
}