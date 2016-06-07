The npctransport module is a module for simulating transport through the NPC.

Authors: Barak Raveh, Daniel Russel

License:

Publications: Timney*, Raveh* et al., in review

Building from source in a nutshell (see https://integrativemodeling.org/latest/doc/manual/installation.html for general IMP installation instruction): 
- Make sure IMP prerequisites are installed
- Install the CGAL package 
- Install the protobuf package
- Download IMP source code as explained in the installation manual
- add npctransport as a subfolder in the 'modules' subfolder (e.g. "cd <IMP-REPOSITORY>/modules && git clone https://github.com/salilab/npctransport/"
- Build IMP according to online instructions as usual
- The npctransport installation can be tested using ctest, eg - "ctest -R npctransport"
