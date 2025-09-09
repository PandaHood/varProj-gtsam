//
// Created by nikolas on 9/5/25.
//
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
#include "../json.hpp"
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
    Key Ti = Symbol('T', meas.i);
    Key Lj = Symbol(tag_j, meas.j);

    if (d == 2) {
      graph.emplace_shared<BAFactor<2>>(Xi, Ti,Lj, meas.l, p, noise);
    } else {
      graph.emplace_shared<BAFactor<3>>(Xi, Ti,Lj, meas.l, p, noise);
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

  Values initial;

  // Poses: random Stiefel (p x d) + random t in R^p
  for (const auto& [tag, idx] : poseKeys) {
    // Stiefel random (constructor provided by your class)
    ScaledStiefelKPD Y = ScaledStiefelKPD::Random(idx, d,p,false);

    // Important: avoid inserting temporary expressions directly (clangd/traits issue)
    Eigen::VectorXd t_p(p);
    t_p.setRandom(); // uniform in [-1, 1]

    initial.insert(Symbol(tag, idx), Y);
    initial.insert(Symbol('T', idx), t_p);
  }

  // Landmarks: random vector in R^p
  for (const auto& [tag, idx] : landmarkKeys) {
    Eigen::VectorXd v(p);
    v.setRandom();
    initial.insert(Symbol(tag, idx), v);
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

  // -----------------------------
  // Optimize with LM
  // -----------------------------
  auto lmParams = LevenbergMarquardtParams::CeresDefaults();
  lmParams.maxIterations    = 200;      // tweak as you wish
  lmParams.relativeErrorTol = 1e-10;    // tighter tolerances if needed
  lmParams.verbosityLM      = LevenbergMarquardtParams::SUMMARY;

  std::cout << "[Optimize] Starting LM with " << graph.size()
            << " factors and " << initial.size() << " variables...\n";

  auto t0 = clock_type::now();
  LevenbergMarquardtOptimizer lm(graph, initial, lmParams);
  Values results = lm.optimize();
  auto t1 = clock_type::now();

  const double solve_time_sec = std::chrono::duration<double>(t1 - t0).count();

  std::cout << std::fixed << std::setprecision(6);
  std::cout << "[Optimize] Done. iterations=" << lm.iterations()
            << "  final_error=" << lm.error()
            << "  time=" << solve_time_sec << "s\n";

  return 0;
}