#!/bin/bash

# set the following three variables
# embree_library_path=not required on macOS
intel_tbb_library_path=
spray_home_path=
spray_bin_path=

export DYLD_LIBRARY_PATH=$intel_tbb_library_path:$DYLD_LIBRARY_PATH
export SPRAY_HOME_PATH=$spray_home_path
export SPRAY_BIN_PATH=$spray_bin_path
