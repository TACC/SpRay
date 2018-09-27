
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


## Building Spray (Linux and Mac)

If you are using Apple Clang on a Mac, install `libomp`.
```bash
brew install libomp
```
Build and install [Embree][2]. We tested Spray with Embree v2.17.1.

Check out Spray with all submodules.
```bash
git clone https://github.com/TACC/SpRay.git
cd SpRay
git submodule init
git submodule update
```
Or simply,
```bash
git clone --recurse-submodules https://github.com/TACC/SpRay.git
cd SpRay
```
Create a build directory.
```bash
mkdir build
cd build
```
If you are using Apple Clang on a Mac,
```bash
cmake -DOpenMP_INSTALL_DIR=$(brew --prefix libomp) -DEMBREE_INSTALL_DIR=<path_to_embree_install> ..
```
Otherwise,
```bash
cmake -DEMBREE_INSTALL_DIR=<path_to_embree_install> ..
```
```bash
make
```
If you wish to install Spray,
```bash
make install
```

## Running Spray (Linux and Mac)

### Rendering isosurfaces of 64 domains

```bash
export SPRAY_HOME_PATH=<path_to_spray_home>
```
For film mode,
```bash
source $SPRAY_HOME_PATH/examples/wavelets64/wavelets64_film.sh
display spray.ppm
```

Notice that this bash script launches 2 MPI tasks using `mpirun -n 2`. You may have to modify the MPI command as needed.

For glfw mode,
```bash
source $SPRAY_HOME_PATH/examples/wavelets64/wavelets64_glfw.sh
```
Type the `q-key` to close the window.

If you wish to use installed binaries, you'll have to modify the variable `SPRAY_BIN` in the bash scripts and set runtime search paths, LD_LIBRARY_PATH or DYLD_LIBRARY_PATH, as needed.

You should see the following as a result:

![wavelets.jpg](images/wavelets64.jpg)



[1]: https://www.apache.org/licenses/LICENSE-2.0
[2]: https://github.com/embree/embree

