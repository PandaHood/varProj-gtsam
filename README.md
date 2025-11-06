# GTSAM Experiments for the paper Sparse Variable Projection in Robotic Perception: Exploiting Separable Structure for Efficient Nonlinear Optimization

This is the code for the gtsam experiments in the paper Sparse Variable Projection in Robotic Perception: Exploiting SeparableStructure for Efficient Nonlinear Optimization. For the rest of the experiments please see ___

## Install

The C++ implementation can be built and exported as a CMake project.

#### C++ quick installation guide

The following installation instructions have been verified on Ubuntu 22.04:

*Step 1:*  Install dependencies

```
sudo apt-get install liblapack-dev libblas-dev libsuitesparse-dev
```
[Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page)  (Make sure using the same version o external as GTSAM's to avoid conflicts, i.e., you can enable `-DGTSAM_USE_SYSTEM_EIGEN=ON` when compile GTSAM)

[GTSAM](https://github.com/borglab/gtsam) (tested on release/[4.3a0](https://github.com/borglab/gtsam/releases/tag/4.3a0), should also compatible with a series of versions that after **Boost Removal**, we use OptionalMatrixType but not boost optional matrix for customer factors)

(If you have multiple GTSAM versions in your machine, and do not want to destroy the environment. The easiest thing to do may be just specify the path to your GTSAM compatible with our project by :

```
set(GTSAM_DIR "path to your /gtsam/build")
set(GTSAM_UNSTABLE_DIR "path to your /gtsam/build")
## Above your find_package(GTSAM REQUIRED) and :
set(GTSAMCMakeTools_DIR "path to your /gtsam/cmake")
## Above your find_package(GTSAMCMakeTools REQUIRED)
```

)
To run the experiments:
```bash
git clone git@github.com:PandaHood/varProj-gtsam.git
cd varProj-gtsam 
mkdir cmake-build-default
cd cmake-build-default
cmake ..
make -j
```

You can then run the bash scripts inside the data folders (Provided you have added the data).

```
./data/raslam/run_raslam.sh
./data/sfm/run_sfm.sh
./data/snl/run_snl.sh
```

