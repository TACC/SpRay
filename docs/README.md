# SpRay: a distributed-memory speculative ray tracer for out-of-core and in situ rendering

## Introduction

This repository contains the source code of the speculative ray scheduling technique described in the following paper:
```
SpRay: Speculative Ray Scheduling for Large Data Visualization
Hyungman Park, Donald Fussell, Paul Navratil
IEEE Symposium on Large Data Analysis and Visualization 2018
```

It also includes our implementations of the baseline algorithm described in the following paper:
```
Exploring the Spectrum of Dynamic Scheduling Algorithms for Scalable Distributed-Memory Ray Tracing
Paul A. Navratil, Hank Childs, Donald S. Fussell, Calvin Lin
IEEE Transactions on Visualization and Computer Graphics 2013
```

You can find our paper and slides on <a href="https://hyungman.bitbucket.io/projects/spray/" target="_blank">[the project page]</a>.

## Documentation

Please read our preliminary documentation to get started with SpRay. You can refer to either the `docs` directory or <a href="https://TACC.github.io/SpRay/" target="_blank">the documentation page</a>.

* [Building SpRay on Linux](build_linux.md)
* [Building SpRay on macOS](build_mac.md)
* [Example 1: rendering a scene of 64 wavelet domains](example1.md)
* [Example 2: generating a scene file and rendering the scene](example2.md)

## Contributors
* Hyungman Park, ECE and TACC, UT Austin (Developer)
* Paul Navratil, TACC, UT Austin (Advisor)
* Donald Fussell, CS, UT Austin (Advisor)

## License
Copyright (c) 2017-2018 The University of Texas at Austin. All rights reserved.

SpRay is licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License. A copy of the License is included with this software in the file `LICENSE`. If your copy does not contain the License, you may obtain a copy of the License at: [Apache License Version 2.0][1].
Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.  

## Acknowledgments
* National Science Foundation grant ACI-1339863
* An Intel Visualization Center of Excellence award through the IPCC program

[1]: https://www.apache.org/licenses/LICENSE-2.0
[2]: https://github.com/embree/embree
[3]: https://www.cs.utexas.edu/~lin/papers/tvcg13.pdf
[4]: https://hyungman.bitbucket.io/projects/spray/
