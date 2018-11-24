# SpRay: a distributed-memory speculative ray tracer for out-of-core and in situ rendering

## Introduction

This repository contains the source code of the speculative ray scheduling technique described in the following [paper][ldav-paper]:
```
SpRay: Speculative Ray Scheduling for Large Data Visualization
Hyungman Park, Donald Fussell, Paul Navratil
IEEE Symposium on Large Data Analysis and Visualization 2018
```

It also includes our implementations of the baseline algorithm described in the following [paper][tvcg-paper]:
```
Exploring the Spectrum of Dynamic Scheduling Algorithms for Scalable Distributed-Memory Ray Tracing
Paul A. Navratil, Hank Childs, Donald S. Fussell, Calvin Lin
IEEE Transactions on Visualization and Computer Graphics 2013
```

You can find our [LDAV 2018 paper][ldav-paper] and [talk slides][ldav-slides] on the [project page][spray-project-page].

## Documentation

Please read our preliminary documentation to get started with SpRay.

* [Building SpRay on Linux](build_linux.md)
* [Building SpRay on macOS](build_mac.md)
* [Example 1: rendering a scene of 64 wavelet domains](example1.md)
* [Example 2: generating a scene file and rendering the scene](example2.md)

## Contributors
* [Hyungman Park][park], [ECE][utece] and [TACC][tacc], [UT Austin][utexas] (Developer)
* [Paul Navratil][navratil], [TACC][tacc], [UT Austin][utexas] (Advisor)
* [Donald Fussell][fussell], [CS][utcs], [UT Austin][utexas] (Advisor)

## License
Copyright (c) 2017-2018 The University of Texas at Austin. All rights reserved.

SpRay is licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License. A copy of the License is included with this software in the file `LICENSE`. If your copy does not contain the License, you may obtain a copy of the License at: [Apache License Version 2.0][apache].
Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.  

## Acknowledgments
* National Science Foundation grant ACI-1339863
* An Intel Visualization Center of Excellence award through the IPCC program


[ldav-paper]: https://hyungman.bitbucket.io/projects/spray/spray_preprint.pdf
[tvcg-paper]: https://www.cs.utexas.edu/~lin/papers/tvcg13.pdf
[ldav-slides]: https://hyungman.bitbucket.io/projects/spray/ldav18_spray_slides.pdf
[spray-project-page]: https://hyungman.bitbucket.io/projects/spray/
[spray-doc]: https://tacc.github.io/SpRay/
[park]: https://hyungman.bitbucket.io/
[navratil]: http://pages.tacc.utexas.edu/~pnav/
[fussell]: https://www.cs.utexas.edu/users/fussell/
[utece]: http://www.ece.utexas.edu/
[tacc]: https://www.tacc.utexas.edu/
[utcs]: https://www.cs.utexas.edu/
[utexas]: https://www.utexas.edu/
[apache]: https://www.apache.org/licenses/LICENSE-2.0

