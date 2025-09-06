#!/bin/bash

# Declare array with arguments formatted correctly
args_list=(
    "2 2 data/range/goats_14 data/range/goats_14_test"
    "2 2 data/range/goats_15 data/range/goats_15_test"
    "2 2 data/range/goats_16 data/range/goats_16_test"
    "2 2 data/range/marine data/range/marine_test"
    "2 2 data/range/plaza1 data/range/plaza1_test"
    "2 2 data/range/plaza2 data/range/plaza2_test"
    "3 3 data/range/single_drone data/range/single_drone_test"
)

# Correct variable assignment without spaces
SESYNC_PATH="./cmake-build-release/Preconditioners/bin/Certifiable_RangeSLAM_example"

# Loop through the array and run the command for each set of arguments
for args in "${args_list[@]}"; do
    # Split each string into its respective parts
    read -r d p dataset output <<< "$args"

    # Echo the arguments for debugging purposes
    echo "Running with arguments: d=$d, p=$p, dataset=$dataset, output=$output"

    # Call the program with the correct arguments
    $SESYNC_PATH "$d" "$p" "$dataset" "$output"
done