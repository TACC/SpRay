#!/bin/bash

# set the following four variables
embree_library_path=
intel_tbb_library_path=
spray_home_path=
spray_bin_path=

export LD_LIBRARY_PATH=$embree_library_path:$intel_tbb_library_path:$LD_LIBRARY_PATH
export SPRAY_HOME_PATH=$spray_home_path
export SPRAY_BIN_PATH=$spray_bin_path
