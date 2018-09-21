
# SpRay: a distributed-memory speculative ray tracer for out-of-core and in situ rendering

## Copyright

Copyright (c) 2017-2018 The University of Texas at Austin. All rights reserved.

SpRay is licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License. A copy of the License is included with this software in the file `LICENSE`. If your copy does not contain the License, you may obtain a copy of the License at: [Apache License Version 2.0][1]. 
Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.  

## Contributors

  * Hyungman Park, ECE and TACC, UT Austin (Developer)
  * Paul Navratil, TACC, UT Austin (Advisor)
  * Donald Fussell, CS, UT Austin (Advisor)

## Funding
 * National Science Foundation grant ACI-1339863
 * An Intel Visualization Center of Excellence award through the IPCC program


## Building Spray

Build and install Embree. We tested Spray with Embree v2.17.1.
```bash
git clone https://github.com/embree/embree.git
```
Check out Spray with all submodules.
```bash
git clone https://github.com/TACC/SpRay.git
git submodule init --update --recursive
```
Or simply,
```bash
git clone --recursive https://github.com/TACC/SpRay.git
```

Create a build directory and run cmake.
```bash
cd SpRay
mkdir build
cd build
cmake -DEMBREE_INSTALL_DIR=<path_to_embree_install> -DCMAKE_INSTALL_PREFIX=<path_to_spray_install> ..
make
```

If you wish to install Spray,
```bash
make install
```

## Running Spray

### Rendering isosurfaces of 64 domains

```bash
export SPRAY_HOME_PATH=<path_to_spray_home>
```
For film mode,
```bash
source $SPRAY_HOME_PATH/examples/wavelets64/wavelets64_film.sh
display spray.ppm
```

Notice that this bash script launches 4 MPI tasks using `mpirun -n 4`. You may have to modify the MPI command.

For glfw mode,
```bash
source $SPRAY_HOME_PATH/examples/wavelets64/wavelets64_glfw.sh
```
Type the `q-key` to close the window.

If you wish to use installed binaries, you'll have to modify the variable `SPRAY_BIN` in the bash scripts and set runtime search paths for the shared libraries used.

You should see the following as a result:

![wavelets.jpg](images/wavelets64.jpg)



[1]: https://www.apache.org/licenses/LICENSE-2.0

