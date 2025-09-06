#!/bin/bash

# Declare array with arguments formatted correctly
args_list=(
    "2 2 /home/nikolas/StiefelManifold/data/pgo/ais2klinik.g2o /home/nikolas/StiefelManifold/results/CPGO/ais2kliniktest_results"
    "2 2 /home/nikolas/StiefelManifold/data/pgo/city10000.g2o /home/nikolas/StiefelManifold/results/CPGO/city10000test_results"
    "2 2 /home/nikolas/StiefelManifold/data/pgo/CSAIL.g2o /home/nikolas/StiefelManifold/results/CPGO/CSAIL_results"
    "2 2 /home/nikolas/StiefelManifold/data/pgo/intel.g2o /home/nikolas/StiefelManifold/results/CPGO/intel_results"
    "2 2 /home/nikolas/StiefelManifold/data/pgo/kitti_00.g2o /home/nikolas/StiefelManifold/results/CPGO/kitti_00_results"
    "2 2 /home/nikolas/StiefelManifold/data/pgo/kitti_02.g2o /home/nikolas/StiefelManifold/results/CPGO/kitti_02_results"
    "2 2 /home/nikolas/StiefelManifold/data/pgo/kitti_05.g2o /home/nikolas/StiefelManifold/results/CPGO/kitti_05_results"
    "2 2 /home/nikolas/StiefelManifold/data/pgo/manhattan.g2o /home/nikolas/StiefelManifold/results/CPGO/manhattan_results"
    "2 2 /home/nikolas/StiefelManifold/data/pgo/MIT.g2o /home/nikolas/StiefelManifold/results/CPGO/MIT_results"
    "3 3 /home/nikolas/StiefelManifold/data/pgo/cubicle.g2o /home/nikolas/StiefelManifold/results/CPGO/cubicle_results"
    "3 3 /home/nikolas/StiefelManifold/data/pgo/grid3D.g2o /home/nikolas/StiefelManifold/results/CPGO/grid3D_results"
    "3 3 /home/nikolas/StiefelManifold/data/pgo/parking-garage.g2o /home/nikolas/StiefelManifold/results/CPGO/parkinggarage_results"
    "3 3 /home/nikolas/StiefelManifold/data/pgo/rim.g2o /home/nikolas/StiefelManifold/results/CPGO/rim_results"
    "3 3 /home/nikolas/StiefelManifold/data/pgo/sphere2500.g2o /home/nikolas/StiefelManifold/results/CPGO/sphere2500_results"
    "3 3 /home/nikolas/StiefelManifold/data/pgo/smallGrid3D.g2o /home/nikolas/StiefelManifold/results/CPGO/smallGrid3D_results"
    "3 3 /home/nikolas/StiefelManifold/data/pgo/torus3D.g2o /home/nikolas/StiefelManifold/results/CPGO/torus3D_results"
)

# Correct variable assignment without spaces
SESYNC_PATH="/home/nikolas/StiefelManifold/cmake-build-release/Preconditioners/bin/Certifiable_PGO_example"

# Loop through the array and run the command for each set of arguments
for args in "${args_list[@]}"; do
    # Split each string into its respective parts
    read -r d p dataset output <<< "$args"

    # Echo the arguments for debugging purposes
    echo "Running with arguments: d=$d, p=$p, dataset=$dataset, output=$output"

    # Call the program with the correct arguments
    $SESYNC_PATH "$d" "$p" "$dataset" "$output"
done
