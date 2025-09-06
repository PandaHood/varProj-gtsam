//
// Created by nikolas on 9/3/25.
//

#include <iomanip>
#include <gtsam/base/timing.h>
//#include <gtsam/sfm/ShonanAveraging.h>
#include <gtsam/slam/InitializePose.h>
#include <gtsam/slam/dataset.h>
#include "../utils.h"
#include "../RaFactor.h"
#include "../LiftedPose.h"
#include "../SEsyncFactor.h"
#include <gtsam/nonlinear/LevenbergMarquardtOptimizer.h>
#include <iostream>
#include <fstream>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <limits>
#include <filesystem>
#include <nlohmann/json.hpp>
#include <chrono>


using clock_type = std::chrono::steady_clock; // monotonic, not affected by NTP
using nlohmann::json;
using namespace std;
using namespace gtsam;

static std::string datasetKeyFromPath(const std::string& input_file) {
  std::filesystem::path p(input_file);
  return p.filename().string(); // e.g., "goats_16.pyfg"
}

// Writes:
// results[<dataset>]["rank<rank>"]["rank<rank>_initN.txt"] = {cost, iterations, time, formulation}
void appendSESyncPGORunToJson(const std::string& out_path,
                              const std::string& input_file,
                              int rank,
                              int iterations,
                              double cost,
                              double time_sec,
                              const std::string& formulation) // e.g. "SESync" or "PGO"
{
  // 1) Load or start fresh
  json root;
  {
    std::ifstream ifs(out_path);
    if (ifs.good()) {
      try { ifs >> root; } catch (...) { root = json::object(); }
    } else {
      root = json::object();
    }
  }

  // 2) Ensure containers
  if (!root.contains("results") || !root["results"].is_object())
    root["results"] = json::object();

  const std::string dsKey   = datasetKeyFromPath(input_file);           // "goats_16.pyfg"
  const std::string rankKey = "rank" + std::to_string(rank);            // "rank3"

  json& dsObj = root["results"][dsKey];
  if (!dsObj.is_object()) dsObj = json::object();

  // If a legacy "summary" exists, remove it to keep the requested shape clean
  if (dsObj.contains("summary")) dsObj.erase("summary");

  json& rankObj = dsObj[rankKey];
  if (!rankObj.is_object()) rankObj = json::object();

  // 3) Find next available init index: rank<r>_init1.txt, _init2.txt, ...
  int next_init = 1;
  for (;; ++next_init) {
    std::string candidate = rankKey + "_init" + std::to_string(next_init) + ".txt";
    if (!rankObj.contains(candidate)) break;
  }
  const std::string initKey = rankKey + "_init" + std::to_string(next_init) + ".txt";

  // 4) Write this run
  rankObj[initKey] = {
    {"cost",        cost},
    {"iterations",  iterations},
    {"time",        time_sec},
    {"formulation", formulation}
  };

  // 5) Save (pretty)
  std::filesystem::create_directories(std::filesystem::path(out_path).parent_path());
  std::ofstream ofs(out_path);
  ofs << std::setw(2) << root << '\n';
}

LiftedPoseDP LiftedToP(const Pose3 &pose3_, const size_t p) {
    return LiftedPoseDP(StiefelManifoldKP::Lift(p, pose3_.rotation().matrix()), pose3_.translation());
}

int main(int argc, char* argv[]) {

    /*
    * argv[1] = d
    * argv[2] = p
    * argv[3] = dataset
    * argv[4] = output file
    *
    * */

    std::string inputFile;
    std::string outputFile;

    // Parse spatial dimension (d) and lifted dimension (p) from arguments
    int d = std::stoi(argv[1]);
    int p = std::stoi(argv[2]);

    inputFile  = std::string(argv[3]);  // Path to the input dataset

    // Range of lifted dimensions to test
    size_t pMin = p;
    size_t pMax = d + 10;

    // Number of poses in the dataset
    size_t num_poses;

    auto measurements = DataParser::read_g2o_file(inputFile, num_poses);
    std::cout << "Loaded " << measurements.poseMeasurements.size() << " measurements between "
              << num_poses << " poses from file "  << inputFile << std::endl;

    // Random generator from Dave's code
    std::default_random_engine generator(std::default_random_engine::default_seed);
    std::normal_distribution<double> g;

    NonlinearFactorGraph inputGraph;

    for (const auto & meas : measurements.poseMeasurements) {
        double Kappa = meas.kappa;
        double tau = meas.tau;
        Vector sigmas = Vector::Zero(p * d + p);
        sigmas.head(p * d).setConstant(sqrt(1/ ( 2 * Kappa)));
        sigmas.tail(p).setConstant( sqrt(1/  (2 * tau)));
        noiseModel::Diagonal::shared_ptr noise = noiseModel::Diagonal::Sigmas(sigmas);


        if (d == 2) {
            inputGraph.emplace_shared<SEsyncFactor2>(meas.i, meas.j, meas.R, meas.t, p, noise);

        }
        else if (d ==3 ) {
            inputGraph.emplace_shared<SEsyncFactor3>(meas.i, meas.j, meas.R, meas.t, p, noise);

        }
        else {
            std::cerr << "Un" << std::endl;
        }

    }

    NonlinearFactorGraph::shared_ptr graph;
    Values::shared_ptr initials;
    bool is3D = true;
    std::tie(graph, initials) = readG2o(inputFile, is3D);

    Values initials_g2o;
    for (size_t k = 0; k < initials->size(); k++) {
        Pose3 retrieved_pose;
        if (initials->exists(k))
            retrieved_pose = initials->at<gtsam::Pose3>(k);
        initials_g2o.insert(k, LiftedToP(retrieved_pose, p));
    }

    Values initial;
    for (size_t j = 0; j < num_poses; j++) {
        StiefelManifoldKP Y = StiefelManifoldKP::Random(std::default_random_engine::default_seed, d, p);

        Vector trans = Vector::Zero(p);
        for (int i = 0; i < p; i++) {
            trans(i) = g(generator);
        }
        initial.insert(j, LiftedPoseDP(Y, trans));
    }


    Values::shared_ptr posesInFile;
    Values poses;
    auto lmParams = LevenbergMarquardtParams::CeresDefaults();
//    LevenbergMarquardtParams lmParams;
    lmParams.maxIterations = 1000;
    lmParams.relativeErrorTol = 1e-5;
    lmParams.verbosityLM = LevenbergMarquardtParams::SUMMARY;

    auto lm = std::make_shared<LevenbergMarquardtOptimizer>(inputGraph, initial, lmParams);
    auto t0 = clock_type::now();
    auto results = lm->optimize();
    auto t1 = clock_type::now();

    double solve_time_sec =
        std::chrono::duration<double>(t1 - t0).count(); // seconds as double
    // after optimize():
    const int iterations = lm->iterations();
    const double final_cost = lm->error();             // or graph.error(results) if you prefer
    appendSESyncPGORunToJson("/home/nikolas/StiefelManifold/results/PGO/results.json", inputFile, p, iterations, final_cost, solve_time_sec, "gtsam"); // or "PGO"



    return 0;
}