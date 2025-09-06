//
// Created by jason on 1/27/25.
//

/**
 * @file utils.cpp
 * @date
 * @author Jason Xu, Nikolas Sanderson
 * @brief
 */

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <Eigen/Geometry>
#include "utils.h"
#include <unordered_set>
#include <set>

namespace gtsam {
    namespace DataParser {
        Matrix fromAngle(double angle_rad) {
            Matrix rotation_matrix_2d(2, 2);
            rotation_matrix_2d << cos(angle_rad), -sin(angle_rad), sin(angle_rad),
                    cos(angle_rad);
            return rotation_matrix_2d;
        }

        Matrix fromQuat(double qx, double qy, double qz, double qw) {
            Eigen::Quaterniond q(qw, qx, qy, qz);
            auto rot_mat = q.toRotationMatrix();
            // Not sure why we can't cast it directly?
            Matrix result(3, 3);
            result << rot_mat(0, 0), rot_mat(0, 1), rot_mat(0, 2), rot_mat(1, 0),
                    rot_mat(1, 1), rot_mat(1, 2), rot_mat(2, 0), rot_mat(2, 1), rot_mat(2, 2);
            return result;
        }

        Measurement read_g2o_file(const std::string &filename, size_t &num_poses_) {
            std::unordered_set<size_t> pose_ids, landmark_ids;

            // Preallocate output vector
            Measurement measurements;
            RelativePoseMeasurement posemeasurement;
            RelativeLandmarkMeasurement landmarkmeasurement;


            // A string used to contain the contents of a single line
            std::string line;

            // A string used to extract tokens from each line one-by-one
            std::string token;

            // Preallocate various useful quantities
            Scalar dx, dy, dz, dtheta, dqx, dqy, dqz, dqw, I11, I12, I13, I14, I15, I16,
                    I22, I23, I24, I25, I26, I33, I34, I35, I36, I44, I45, I46, I55, I56, I66;

            size_t i, j;

            // Open the file for reading
            std::ifstream infile(filename);
            size_t &num_poses = measurements.num_poses;
            size_t &num_landmarks = measurements.num_landmarks;

            std::unordered_map<size_t, size_t> poses;
            std::unordered_map<size_t, size_t> landmarks;

            while (std::getline(infile, line)) {
                // Construct a stream from the string
                std::stringstream strstrm(line);

                // Extract the first token from the string
                strstrm >> token;

                if (token == "EDGE_SE2") {
                    // This is a 2D pose measurement

                    /** The g2o format specifies a 2D relative pose measurement in the
                     * following form:
                     *
                     * EDGE_SE2 id1 id2 dx dy dtheta, I11, I12, I13, I22, I23, I33
                     *
                     */

                    // Extract formatted output
                    strstrm >> i >> j >> dx >> dy >> dtheta >> I11 >> I12 >> I13 >> I22 >>
                            I23 >> I33;
                    if (poses.insert({i, num_poses}).second) num_poses++;

                    if (poses.insert({j, num_poses}).second) num_poses++;

                    // Pose ids
                    posemeasurement.i = poses[i];
                    posemeasurement.j = poses[j];
                    pose_ids.insert(i);
                    pose_ids.insert(j);

                    // Raw measurements
                    posemeasurement.t = Eigen::Matrix<Scalar, 2, 1>(dx, dy);
                    posemeasurement.R = Eigen::Rotation2D<Scalar>(dtheta).toRotationMatrix();

                    Eigen::Matrix<Scalar, 2, 2> TranInfo;
                    TranInfo << I11, I12, I12, I22;
                    posemeasurement.tau = 2 / TranInfo.inverse().trace();

                    posemeasurement.kappa = I33;



                    measurements.poseMeasurements.push_back(posemeasurement);

                } else if (token == "EDGE_SE3:QUAT") {

                    // This is a 3D pose measurement

                    /** The g2o format specifies a 3D relative pose measurement in the
                     * following form:
                     *
                     * EDGE_SE3:QUAT id1, id2, dx, dy, dz, dqx, dqy, dqz, dqw
                     *
                     * I11 I12 I13 I14 I15 I16
                     *     I22 I23 I24 I25 I26
                     *         I33 I34 I35 I36
                     *             I44 I45 I46
                     *                 I55 I56
                     *                     I66
                     */

                    // Extract formatted output
                    strstrm >> i >> j >> dx >> dy >> dz >> dqx >> dqy >> dqz >> dqw >> I11 >>
                            I12 >> I13 >> I14 >> I15 >> I16 >> I22 >> I23 >> I24 >> I25 >> I26 >>
                            I33 >> I34 >> I35 >> I36 >> I44 >> I45 >> I46 >> I55 >> I56 >> I66;

                    // Fill in elements of the measurement

                    // Pose ids
                    if (poses.insert({i, num_poses}).second) num_poses++;

                    if (poses.insert({j, num_poses}).second) num_poses++;
                    posemeasurement.i = poses[i];
                    posemeasurement.j = poses[j];
                    pose_ids.insert(i);
                    pose_ids.insert(j);

                    // Raw measurements
                    posemeasurement.t = Eigen::Matrix<Scalar, 3, 1>(dx, dy, dz);
                    posemeasurement.R =
                            Eigen::Quaternion<Scalar>(dqw, dqx, dqy, dqz).toRotationMatrix();

                    // Compute precisions

                    // Compute and store the optimal (information-divergence-minimizing) value
                    // of the parameter tau
                    Eigen::Matrix<Scalar, 3, 3> TranInfo;
                    TranInfo << I11, I12, I13, I12, I22, I23, I13, I23, I33;
                    posemeasurement.tau = 3 / TranInfo.inverse().trace();
                    posemeasurement.trans_precision = 3 / TranInfo.inverse().trace();
//                    measurement.trans_precision_true << TranInfo.inverse()(0, 0), TranInfo.inverse()(1, 1), TranInfo.inverse()(2, 2);
                    // Compute and store the optimal (information-divergence-minimizing value
                    // of the parameter kappa

                    Eigen::Matrix<Scalar, 3, 3> RotInfo;
                    RotInfo << I44, I45, I46, I45, I55, I56, I46, I56, I66;
                    posemeasurement.kappa = 3 / (2 * RotInfo.inverse().trace());
//                    measurement.kappa = 3 / (RotInfo.inverse().trace());
                    posemeasurement.rot_precision = 3 / (2* RotInfo.inverse().trace());
//                    measurement.rot_precision_true << RotInfo.inverse()(0, 0), RotInfo.inverse()(1, 1), RotInfo.inverse()(2, 2);



                    measurements.poseMeasurements.push_back(posemeasurement);

                } else if (token == "LANDMARK2")
                {
                    strstrm >> i >> j >> dx >> dy >> I11 >> I12 >> I22;

                    if (poses.insert({i, num_poses}).second) num_poses++;

                    if (landmarks.insert({j, num_landmarks}).second) num_landmarks++;

                    landmarkmeasurement.i = poses[i];     // pose index (if you need a mapping there too, do the same)
                    landmarkmeasurement.j = landmarks[j];
                    pose_ids.insert(i);
                    landmark_ids.insert(j);    // the “j” is a landmark

                    landmarkmeasurement.l = Eigen::Matrix<Scalar,2,1>(dx,dy);

                    Eigen::Matrix<Scalar,2,2> TranCov;
                    TranCov << I11, I12,
                               I12, I22;
                    landmarkmeasurement.nu = 2.0 / TranCov.inverse().trace();

                    measurements.landmarkMeasurements.push_back(landmarkmeasurement);
                }else if (token == "LANDMARK3")
                {
                    strstrm >> i >> j >> dx >> dy >> dz >> dqx >> dqy >> dqz >> dqw >> I11 >>
                    I12 >> I13 >> I14 >> I15 >> I16 >> I22 >> I23 >> I24 >> I25 >> I26 >>
                    I33 >> I34 >> I35 >> I36 >> I44 >> I45 >> I46 >> I55 >> I56 >> I66;

                    if (poses.insert({i, num_poses}).second) num_poses++;

                    if (landmarks.insert({j, num_landmarks}).second) num_landmarks++;
                    landmarkmeasurement.i = poses[i];
                    landmarkmeasurement.j = landmarks[j];
                    pose_ids.insert(i);        // the “i” is still a pose
                    landmark_ids.insert(j);    // the “j” is a landmark
                    landmarkmeasurement.l = Eigen::Matrix<Scalar, 3, 1>(dx, dy, dz);

                    Eigen::Matrix<Scalar, 3, 3> TranCov;
                    TranCov << I11, I12, I13, I12, I22, I23, I13, I23, I33;

                    landmarkmeasurement.nu = 3 / (2 * TranCov.inverse().trace());



                    measurements.landmarkMeasurements.push_back(landmarkmeasurement);
                }else if ((token == "VERTEX_SE2") || (token == "VERTEX_SE3:QUAT") || (token == "VERTEX_XY")) {
                    // This is just initialization information, so do nothing
                    continue;
                }
                else {
                    std::cout << "Error: unrecognized type: " << token << "!" << std::endl;
                    assert(false);
                }


            } // while

            infile.close();
            measurements.num_poses = pose_ids.size();
            measurements.num_landmarks = landmark_ids.size();
            return measurements;
        }

         Measurement read_pycfg_file(const std::string &filename) {
            std::ifstream infile(filename);
            if (!infile) {
                throw std::runtime_error("Could not open " + filename);
            }

            // Scalars read from file
            Scalar timestamp, dx, dy, dz, dtheta, dqx, dqy, dqz, dqw;
            Scalar I11, I12, I13, I14, I15, I16, I22, I23, I24, I25, I26, I33, I34, I35, I36, I44, I45, I46, I55, I56, I66;
            Scalar range, sigma;

            Measurement M;

            // Track seen indices to build global reindexing
            std::map<char, std::set<size_t>> robot_pose_indices; // per-robot local pose indices
            std::set<char> robot_identity;                       // which robot letters exist
            std::set<size_t> landmark_ids;                       // landmark local IDs (from L<number>)

            std::string ida, idb, idl; // endpoint tokens as read from file
            std::string line;

            auto parse_id = [](const std::string& s, char& c, size_t& idx) {
                if (s.empty()) throw std::runtime_error("Empty ID token");
                c   = s[0];
                idx = std::stoull(s.substr(1));
            };

            while (std::getline(infile, line)) {
                // Skip blanks / comments
                if (line.empty()) continue;
                // Allow comments starting with '#' or '%'
                char first = 0;
                for (char ch : line) { if (!isspace(static_cast<unsigned char>(ch))) { first = ch; break; } }
                if (first == '#' || first == '%') continue;

                std::stringstream ss(line);
                std::string token;
                ss >> token;
                if (token.empty()) continue;

                if (token == "EDGE_SE2") {
                    ss >> timestamp >> ida >> idb >> dx >> dy >> dtheta
                       >> I11 >> I12 >> I13 >> I22 >> I23 >> I33;

                    char ca, cb; size_t ia, ib;
                    parse_id(ida, ca, ia);
                    parse_id(idb, cb, ib);

                    RelativePoseMeasurement m;
                    m.ci = ca; m.cj = cb;
                    m.i  = ia; m.j  = ib;

                    // track robots/poses
                    robot_identity.insert(ca);
                    robot_identity.insert(cb);
                    robot_pose_indices[ca].insert(ia);
                    robot_pose_indices[cb].insert(ib);

                    // translation (2D)
                    m.t.resize(2);
                    m.t << dx, dy;

                    // rotation from angle
                    m.R.resize(2,2);
                    m.R = Eigen::Rotation2D<Scalar>(dtheta).toRotationMatrix();

                    // translational precision (tau) from 2x2 info trace
                    // information → tau, kappa
                    Eigen::Matrix2d Tinfo;
                    Tinfo << I11, I12,
                             I12, I22;
                    m.tau   = 2.0 / Tinfo.trace();
                    m.kappa = 1.0 / I33;
                    M.poseMeasurements.push_back(m);

                } else if (token == "EDGE_SE3:QUAT") {
                    ss >> timestamp >> ida >> idb
                       >> dx >> dy >> dz >> dqx >> dqy >> dqz >> dqw
                       >> I11 >> I12 >> I13 >> I14 >> I15 >> I16
                       >> I22 >> I23 >> I24 >> I25 >> I26
                       >> I33 >> I34 >> I35 >> I36
                       >> I44 >> I45 >> I46 >> I55 >> I56 >> I66;

                    char ca, cb; size_t ia, ib;
                    parse_id(ida, ca, ia);
                    parse_id(idb, cb, ib);

                    RelativePoseMeasurement m;
                    m.ci = ca; m.cj = cb;
                    m.i  = ia; m.j  = ib;

                    robot_identity.insert(ca);
                    robot_identity.insert(cb);
                    robot_pose_indices[ca].insert(ia);
                    robot_pose_indices[cb].insert(ib);

                    // translation (3D)
                    m.t.resize(3);
                    m.t << dx, dy, dz;

                    // rotation from quaternion (w,x,y,z)
                    m.R.resize(3,3);
                    m.R = Eigen::Quaternion<Scalar>(dqw, dqx, dqy, dqz).toRotationMatrix();

                    // translational precision from 3x3 info block trace

                    Eigen::Matrix<Scalar, 3, 3> Tinfo;
                    Tinfo<< I11, I12, I13, I12, I22, I23, I13, I23, I33;
                    m.tau = 3 / Tinfo.trace();
                    m.trans_precision = 3 / Tinfo.trace();


                    Eigen::Matrix<Scalar, 3, 3> RotInfo;
                    RotInfo << I44, I45, I46, I45, I55, I56, I46, I56, I66;
                    m.kappa = 3 / (2 * RotInfo.trace());

                    m.rot_precision = 3 / ( 2* RotInfo.trace());

                    M.poseMeasurements.push_back(m);

                } else if (token == "EDGE_SE2_XY") {
                    // Robot-to-landmark relative landmark measurement in 2D
                    // Format (per example):
                    // EDGE_SE2_XY <t> <RobotPoseId> <LandmarkId> lx ly I11 I12 I22
                    Scalar lx, ly;
                    ss >> timestamp >> ida >> idl >> lx >> ly >> I11 >> I12 >> I22;

                    char cr, cl; size_t ir, il;
                    parse_id(ida, cr, ir);  // robot pose endpoint
                    parse_id(idl, cl, il);  // landmark endpoint (must be 'L')

                    RelativeLandmarkMeasurement m;
                    m.ci = cr; m.cj = cl;   // keep prefixes ('A','B',...,'L')
                    m.i  = ir; m.j  = il;

                    // track robot pose & landmark usage
                    robot_identity.insert(cr);
                    robot_pose_indices[cr].insert(ir);
                    if (cl == 'L') landmark_ids.insert(il);

                    // measured landmark position relative to robot pose
                    m.l.resize(2);
                    m.l << lx, ly;

                    // convert 2x2 info to scalar precision (nu) via trace
                    Eigen::Matrix<Scalar,2,2> Linfo;
                    Linfo << I11, I12,
                             I12, I22;
                    m.nu = Scalar(2.0) / Linfo.trace();

                    M.landmarkMeasurements.push_back(m);

                } else if (token == "EDGE_SE3_XYZ") {
                    // Robot-to-landmark relative landmark measurement in 3D
                    // Format:
                    // EDGE_SE3_XYZ <t> <RobotPoseId> <LandmarkId> lx ly lz I11 I12 I13 I22 I23 I33
                    Scalar timestamp;
                    std::string ida, idl;
                    Scalar lx, ly, lz;
                    Scalar I11, I12, I13, I22, I23, I33;

                    ss >> timestamp >> ida >> idl >> lx >> ly >> lz
                       >> I11 >> I12 >> I13 >> I22 >> I23 >> I33;

                    char cr, cl; size_t ir, il;
                    parse_id(ida, cr, ir);  // robot pose endpoint
                    parse_id(idl, cl, il);  // landmark endpoint (must be 'L')

                    RelativeLandmarkMeasurement m;
                    m.ci = cr; m.cj = cl;   // keep prefixes ('A','B',...,'L')
                    m.i  = ir; m.j  = il;

                    // track robot pose & landmark usage
                    robot_identity.insert(cr);
                    robot_pose_indices[cr].insert(ir);
                    if (cl == 'L') landmark_ids.insert(il);

                    // measured landmark position relative to robot pose
                    m.l.resize(3);
                    m.l << lx, ly, lz;

                    // convert 3x3 info to scalar precision (nu) via trace (consistent with 2D)
                    Eigen::Matrix<Scalar,3,3> Linfo;
                    Linfo << I11, I12, I13,
                             I12, I22, I23,
                             I13, I23, I33;
                    m.nu = Scalar(3.0) / Linfo.trace();

                    M.landmarkMeasurements.push_back(m);
                } else if (token == "EDGE_RANGE") {
                    // Range between robot↔landmark OR robot↔robot
                    // EDGE_RANGE <t> <id_a> <id_b_or_L> range sigma
                    ss >> timestamp >> ida >> idl >> range >> sigma;

                    char ca, cb; size_t ia, ib;
                    parse_id(ida, ca, ia);
                    parse_id(idl, cb, ib);

                    RangeMeasurement m;
                    m.ci    = ca; m.cj    = cb;
                    m.i     = ia; m.j     = ib;
                    m.range = range;
                    m.sigma = sigma;

                    // Track endpoints
                    // Pose side (ci) is always a robot pose in your convention
                    robot_identity.insert(ca);
                    robot_pose_indices[ca].insert(ia);

                    if (cb == 'L') {
                        landmark_ids.insert(ib);
                    } else {
                        // robot-robot range
                        robot_identity.insert(cb);
                        robot_pose_indices[cb].insert(ib);
                    }

                    M.rangeMeasurements.push_back(m);
                }
                // Unknown tokens are ignored (silently). Add 'else { ... }' if you want warnings.
            }
            infile.close();
            // --- Final bookkeeping WITHOUT reindexing ---
            M.robot_pose_indices = robot_pose_indices;
            M.robot_identity     = robot_identity;

            // Count total poses (sum of per-robot unique local indices)
            // NOTE: you won't use M.num_poses for keying anymore, but keep it if needed.
            size_t total_poses = 0;
            for (const auto& kv : robot_pose_indices) total_poses += kv.second.size();
            M.num_poses     = total_poses;

            // Landmarks are the distinct L IDs you saw in file
            M.num_landmarks = landmark_ids.size();

            // Ranges = exactly how many EDGE_RANGE lines we pushed
            M.num_ranges    = M.rangeMeasurements.size();

            return M;

            return M;
        }



    }// namespace DataParser
} // namespace gtsam