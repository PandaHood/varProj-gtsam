//
// Created by jason on 6/23/25.
//


#include "ConcurrentCertifiableBatchSmoother.h"
#include <gtsam/nonlinear/LinearContainerFactor.h>
#include <gtsam/linear/GaussianJunctionTree.h>
#include <gtsam/base/timing.h>
#include <gtsam/base/debug.h>
#include <gtsam/geometry/Pose2.h>
#include <gtsam/geometry/Pose3.h>

namespace gtsam {

   /* ************************************************************************* */
   void ConcurrentCertifiableBatchSmoother::print(const std::string& s, const KeyFormatter& keyFormatter) const {
     std::cout << s;
     std::cout << "  Factors:" << std::endl;
     for(const NonlinearFactor::shared_ptr& factor: factors_) {
       PrintNonlinearFactor(factor, "    ", keyFormatter);
     }
     theta_.print("Values:\n");
   }

   /* ************************************************************************* */
   bool ConcurrentCertifiableBatchSmoother::equals(const ConcurrentSmoother& rhs, double tol) const {
     const ConcurrentCertifiableBatchSmoother* smoother = dynamic_cast<const ConcurrentCertifiableBatchSmoother*>(&rhs);
     return smoother
            && factors_.equals(smoother->factors_)
            && theta_.equals(smoother->theta_)
            && ordering_.equals(smoother->ordering_)
            && delta_.equals(smoother->delta_)
            && separatorValues_.equals(smoother->separatorValues_);
   }

   /* ************************************************************************* */
   template<typename VALUE, size_t d, class Problem>
   ConcurrentCertifiableBatchSmoother::Result
   ConcurrentCertifiableBatchSmoother::update(
       const NonlinearFactorGraph& newFactors,
       const Values& newTheta,
       const std::optional<std::vector<size_t>>& removeFactorIndices)
   {
     gttic(update);

     // Create the return result meta-data
     Result result;

     // Update all of the internal variables with the new information
     gttic(augment_system);
     {
       // Add the new variables to theta
       theta_.insert(newTheta);

       // Add new variables to the end of the ordering
       for (const auto key : newTheta.keys()) {
         ordering_.push_back(key);
       }

       // Augment Delta
       delta_.insert(newTheta.zeroVectors());

       // Add the new factors to the graph, updating the variable index
       insertFactors(newFactors);

       if (removeFactorIndices)
         removeFactors(*removeFactorIndices);
     }
     gttoc(augment_system);

     if (factors_.size() > 0) {
       // Reorder the system to ensure efficient optimization (and marginalization) performance
       gttic(reorder);
       reorder();
       gttoc(reorder);

       // Optimize the factors using the templated LM
       gttic(optimize);
       result = optimize<VALUE,d, Problem>();
       gttoc(optimize);
     }

     // TODO: The following code does considerable work, much of which could be redundant given the previous optimization step
     // Refactor this code to reduce computational burden

     // Calculate the marginal on the separator from the smoother factors
     if (separatorValues_.size() > 0) {
       gttic(presync);
       updateSmootherSummarization();
       gttoc(presync);
     }

     gttoc(update);
     return result;
   }
   /* ************************************************************************* */
   void ConcurrentCertifiableBatchSmoother::presync() {

     gttic(presync);

     gttoc(presync);
   }

   /* ************************************************************************* */
   void ConcurrentCertifiableBatchSmoother::getSummarizedFactors(NonlinearFactorGraph& summarizedFactors, Values& separatorValues) {

     gttic(get_summarized_factors);

     // Copy the previous calculated smoother summarization factors into the output
     summarizedFactors.push_back(smootherSummarization_);

     // Copy the separator values into the output
     separatorValues.insert(separatorValues_);

     gttoc(get_summarized_factors);
   }

   /* ************************************************************************* */
   void ConcurrentCertifiableBatchSmoother::synchronize(const NonlinearFactorGraph& smootherFactors, const Values& smootherValues,
                                             const NonlinearFactorGraph& summarizedFactors, const Values& separatorValues) {

     gttic(synchronize);

     // Remove the previous filter summarization from the graph
     removeFactors(filterSummarizationSlots_);

     // Insert new linpoints into the values, augment the ordering, and store new dims to augment delta
     for(const auto key: smootherValues.keys()) {
       if(!theta_.exists(key)) {
         // If this a new key for theta_, also add to ordering and delta.
         const auto& value = smootherValues.at(key);
         delta_.insert(key, Vector::Zero(value.dim()));
         theta_.insert(key, value);
         ordering_.push_back(key);
       } else {
         // If the key already existed in theta_, just update.
         const auto& value = smootherValues.at(key);
         theta_.update(key, value);
       }
     }
     for(const auto key: separatorValues.keys()) {
       if(!theta_.exists(key)) {
         // If this a new key for theta_, also add to ordering and delta.
         const auto& value = separatorValues.at(key);
         delta_.insert(key, Vector::Zero(value.dim()));
         theta_.insert(key, value);
         ordering_.push_back(key);
       } else {
         // If the key already existed in theta_, just update.
         const auto& value = separatorValues.at(key);
         theta_.update(key, value);
       }
     }

     // Insert the new smoother factors
     insertFactors(smootherFactors);

     // Insert the new filter summarized factors
     filterSummarizationSlots_ = insertFactors(summarizedFactors);

     // Update the list of root keys
     separatorValues_ = separatorValues;

     gttoc(synchronize);
   }

   /* ************************************************************************* */
   void ConcurrentCertifiableBatchSmoother::postsync() {

     gttic(postsync);

     gttoc(postsync);
   }

   /* ************************************************************************* */
   std::vector<size_t> ConcurrentCertifiableBatchSmoother::insertFactors(const NonlinearFactorGraph& factors) {

     gttic(insert_factors);

     // create the output vector
     std::vector<size_t> slots;
     slots.reserve(factors.size());

     // Insert the factor into an existing hole in the factor graph, if possible
     for(const NonlinearFactor::shared_ptr& factor: factors) {
       size_t slot;
       if(availableSlots_.size() > 0) {
         slot = availableSlots_.front();
         availableSlots_.pop();
         factors_.replace(slot, factor);
       } else {
         slot = factors_.size();
         factors_.push_back(factor);
       }
       slots.push_back(slot);
     }

     gttoc(insert_factors);

     return slots;
   }

   /* ************************************************************************* */
   void ConcurrentCertifiableBatchSmoother::removeFactors(const std::vector<size_t>& slots) {

     gttic(remove_factors);

     // For each factor slot to delete...
     for(size_t slot: slots) {

       // Remove the factor from the graph
       factors_.remove(slot);

       // Mark the factor slot as available
       availableSlots_.push(slot);
     }

     gttoc(remove_factors);
   }

   /* ************************************************************************* */
   void ConcurrentCertifiableBatchSmoother::reorder() {

     // Recalculate the variable index
     variableIndex_ = VariableIndex(factors_);

     KeyVector separatorKeys = separatorValues_.keys();
     ordering_ = Ordering::ColamdConstrainedLast(variableIndex_, KeyVector(separatorKeys.begin(), separatorKeys.end()));

   }

   /* ************************************************************************* */
//   template<typename VALUE, size_t d,class Problem>
//   ConcurrentCertifiableBatchSmoother::Result
//   ConcurrentCertifiableBatchSmoother::optimize() {
//     using namespace gtsam;
//
//     // 1) Prepare Result & initial error check
//     Result result;
//     result.nonlinearVariables = theta_.size() - separatorValues_.size();
//     result.linearVariables    = separatorValues_.size();
//
//     // Create the evaluation point and its error
//     Values evalPoint = theta_.retract(delta_);
//     result.error = factors_.error(evalPoint);
//     if (result.error < parameters_.errorTol)
//       return result;
//
//     // 2) Build an augmented graph that pins separator keys of type VALUE
//     NonlinearFactorGraph graph = factors_;
//     // For debug
////       for (const auto& factor : graph) {
////           // Try the common Pose2 factor types:
////           if (   !std::dynamic_pointer_cast<gtsam::BetweenFactor<gtsam::Pose2>>(factor)
////                  && !std::dynamic_pointer_cast<gtsam::PriorFactor<gtsam::Pose2>>(factor)
////               /* add other Pose2‐based factor casts here if needed */)
////           {
////               // If none of the casts succeeded, it’s “something else”
////
////           }
////       }
//     for (Key k : separatorValues_.keys()) {
//       const VALUE& x = theta_.at<VALUE>(k);
//       graph.emplace_shared<NonlinearEquality<VALUE>>(k, x);
//     }
//
//     // 3) Delegate the entire LM loop to GTSAM’s built‐in optimizer
////     LevenbergMarquardtOptimizer optimizer(graph, theta_, parameters_);
////     Values optimizedValues = optimizer.optimize();
//
//
//     auto problem = std::make_unique<Problem>(d, theta_.size());
//
//
//     Values optimizedValues = problem->runIncrementalBatch(graph, theta_, parameters_);
//
//     // auto Problem = CertifiablePGO3(3,theta_.size());
//     // Values optimizedValues = Problem.runIncrementalBatch(graph,theta_,parameters_);
//
//     // ——— SANITY CHECK: separator keys must not have moved ———
////     const double eps = 1e-9;
////     for (Key k : separatorValues_.keys()) {
////       const VALUE& orig = theta_.at<VALUE>(k);
////       const VALUE& opt  = optimizedValues.at<VALUE>(k);
////       Vector d = orig.localCoordinates(opt);
////       if (d.norm() > eps) {
////         std::ostringstream msg;
////         msg << "Error: separator key " << k
////             << " moved by norm " << d.norm();
////         throw std::runtime_error(msg.str());
////       }
////     }
//
//     // 4) Update internal state exactly as before
//     theta_ = optimizedValues;  // new linearization point for next round
//     delta_.setZero();          // reset all deltas
//
//     // 5) Fill in the rest of Result
////     result.iterations = optimizer.iterations();
//// For debug
////       printValuesKeys(optimizedValues);
////       printGraphKeys(factors_);
//
//     result.lambdas    = 0;                         // LMParams doesn’t expose lambda count
//     result.error      = factors_.error(optimizedValues);
//     return result;
//   }


 //***************************************NEW***************************************
 template<typename VALUE, size_t d, class Problem>
 ConcurrentCertifiableBatchSmoother::Result
 ConcurrentCertifiableBatchSmoother::optimize() {
     using namespace gtsam;

     Result result;
     result.nonlinearVariables = theta_.size() - separatorValues_.size();
     result.linearVariables    = separatorValues_.size();

     // early exit on already-good solution
     Values evalPoint = theta_.retract(delta_);
     result.error = factors_.error(evalPoint);
     if (result.error < parameters_.errorTol)
         return result;

     // Build the graph of nonlinear factors (no separator‐fix factors here)
     NonlinearFactorGraph graph = factors_;

     // ----> call your external batch solver <----
     auto problem = std::make_unique<Problem>(d, theta_.size());
     Values optimizedValues =
             problem->runIncrementalBatch(graph, theta_, parameters_);

     // ----> re-inject (overwrite) all separator keys with their fixed values
     optimizedValues.update(separatorValues_);

     // Commit the new linearization point and reset delta
     theta_    = optimizedValues;
     delta_.setZero();

     // Fill out the rest of the result
     result.lambdas = 0;  // (or copy from problem if it exposes it)
     result.error   = factors_.error(theta_);
     return result;
 }



   /* ************************************************************************* */
   void ConcurrentCertifiableBatchSmoother::updateSmootherSummarization() {

     // The smoother summarization factors are the resulting marginal factors on the separator
     // variables that result from marginalizing out all of the other variables
     // These marginal factors will be cached for later transmission to the filter using
     // linear container factors

     // Create a nonlinear factor graph without the filter summarization factors
     NonlinearFactorGraph graph(factors_);
     for(size_t slot: filterSummarizationSlots_) {
       graph.remove(slot);
     }

     // Get the set of separator keys
     const KeySet separatorKeys = separatorValues_.keySet();

     // Calculate the marginal factors on the separator
     smootherSummarization_ = internal::calculateMarginalFactors(graph, theta_, separatorKeys, parameters_.getEliminationFunction());
   }

   /* ************************************************************************* */
   void ConcurrentCertifiableBatchSmoother::PrintNonlinearFactor(const NonlinearFactor::shared_ptr& factor, const std::string& indent, const KeyFormatter& keyFormatter) {
     std::cout << indent;
     if(factor) {
       if(std::dynamic_pointer_cast<LinearContainerFactor>(factor)) {
         std::cout << "l( ";
       } else {
         std::cout << "f( ";
       }
       for(Key key: *factor) {
         std::cout << keyFormatter(key) << " ";
       }
       std::cout << ")" << std::endl;
     } else {
       std::cout << "{ nullptr }" << std::endl;
     }
   }

   /* ************************************************************************* */
   void ConcurrentCertifiableBatchSmoother::PrintLinearFactor(const GaussianFactor::shared_ptr& factor, const std::string& indent, const KeyFormatter& keyFormatter) {
     std::cout << indent;
     if(factor) {
       std::cout << "g( ";
       for(Key key: *factor) {
         std::cout << keyFormatter(key) << " ";
       }
       std::cout << ")" << std::endl;
     } else {
       std::cout << "{ nullptr }" << std::endl;
     }
   }

   template ConcurrentCertifiableBatchSmoother::Result
   ConcurrentCertifiableBatchSmoother::update<gtsam::Pose2,2,CertifiableIncrementalPGO2>(
       const gtsam::NonlinearFactorGraph&,
       const gtsam::Values&,
       const std::optional<std::vector<size_t>>&);

   template ConcurrentCertifiableBatchSmoother::Result
   ConcurrentCertifiableBatchSmoother::update<gtsam::Pose3,3,CertifiableIncrementalPGO3>(
       const gtsam::NonlinearFactorGraph&,
       const gtsam::Values&,
       const std::optional<std::vector<size_t>>&);

   /* ************************************************************************* */
 }/// namespace gtsam
