\brief Simulation of transport through the Nuclear Pore Complex

[![Build Status](https://github.com/salilab/npctransport/workflows/build/badge.svg?branch=develop)](https://github.com/salilab/npctransport/actions?query=workflow%3Abuild)
[![codecov](https://codecov.io/gh/salilab/npctransport/branch/develop/graph/badge.svg)](https://codecov.io/gh/salilab/npctransport)

The npctransport module is a module for simulating transport through the NPC.

_Authors_: Barak Raveh, Daniel Russel

_License_: [LGPL](http://www.gnu.org/licenses/old-licenses/lgpl-2.1.html)
This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2 of the License, or (at your option) any later version.

_Publications_:
Timney*, Raveh* et al., JCB 2016

Building from source in a nutshell (see https://integrativemodeling.org/latest/doc/manual/installation.html for general IMP installation instruction):
- Make sure IMP prerequisites are installed
- Install the CGAL package
- Install the protobuf package
- Download IMP source code as explained in the installation manual
- Build IMP according to online instructions
- The npctransport installation can be tested using ctest, eg - "ctest -R npctransport"

_Versions_:
4.5
- optimization of ball size used to define close pair range results in faster runs for larger particles

fg_simulation: run NPC transport simulations {#fg_simulation_bin}
============================================

The `fg_simulation` command line tool can be used to run simulations of
NPC FG repeat domains using this module. For an example of its use, see
[our 2018 study in Nature](https://salilab.org/npc_fg_2018).
