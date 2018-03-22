The npctransport module is a module for simulating transport through the NPC.

Authors: Barak Raveh, Daniel Russel

License:

Publications:
Timney*, Raveh* et al., JCB 2016

Building from source in a nutshell (see https://integrativemodeling.org/latest/doc/manual/installation.html for general IMP installation instruction):
- Make sure IMP prerequisites are installed
- Install the CGAL package
- Install the protobuf package
- Download IMP source code as explained in the installation manual
- Build IMP according to online instructions
- The npctransport installation can be tested using ctest, eg - "ctest -R npctransport"

Versions:
4.5
- optimization of ball size used to define close pair range results in faster runs for larger particles
