# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /usr/local/google/home/drussel/extern_src/avro-cpp-1.7.2

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /usr/local/google/home/drussel/extern_src/avro-cpp-1.7.2

# Include any dependencies generated for this target.
include CMakeFiles/testgentest.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/testgentest.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/testgentest.dir/flags.make

CMakeFiles/testgentest.dir/test/testgentest.cc.o: CMakeFiles/testgentest.dir/flags.make
CMakeFiles/testgentest.dir/test/testgentest.cc.o: test/testgentest.cc
	$(CMAKE_COMMAND) -E cmake_progress_report /usr/local/google/home/drussel/extern_src/avro-cpp-1.7.2/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/testgentest.dir/test/testgentest.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/testgentest.dir/test/testgentest.cc.o -c /usr/local/google/home/drussel/extern_src/avro-cpp-1.7.2/test/testgentest.cc

CMakeFiles/testgentest.dir/test/testgentest.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/testgentest.dir/test/testgentest.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /usr/local/google/home/drussel/extern_src/avro-cpp-1.7.2/test/testgentest.cc > CMakeFiles/testgentest.dir/test/testgentest.cc.i

CMakeFiles/testgentest.dir/test/testgentest.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/testgentest.dir/test/testgentest.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /usr/local/google/home/drussel/extern_src/avro-cpp-1.7.2/test/testgentest.cc -o CMakeFiles/testgentest.dir/test/testgentest.cc.s

CMakeFiles/testgentest.dir/test/testgentest.cc.o.requires:
.PHONY : CMakeFiles/testgentest.dir/test/testgentest.cc.o.requires

CMakeFiles/testgentest.dir/test/testgentest.cc.o.provides: CMakeFiles/testgentest.dir/test/testgentest.cc.o.requires
	$(MAKE) -f CMakeFiles/testgentest.dir/build.make CMakeFiles/testgentest.dir/test/testgentest.cc.o.provides.build
.PHONY : CMakeFiles/testgentest.dir/test/testgentest.cc.o.provides

CMakeFiles/testgentest.dir/test/testgentest.cc.o.provides.build: CMakeFiles/testgentest.dir/test/testgentest.cc.o

# Object files for target testgentest
testgentest_OBJECTS = \
"CMakeFiles/testgentest.dir/test/testgentest.cc.o"

# External object files for target testgentest
testgentest_EXTERNAL_OBJECTS =

testgentest: CMakeFiles/testgentest.dir/test/testgentest.cc.o
testgentest: libavrocpp.so.1.7.2.0
testgentest: /usr/lib/libboost_filesystem-mt.so
testgentest: /usr/lib/libboost_system-mt.so
testgentest: /usr/lib/libboost_program_options-mt.so
testgentest: CMakeFiles/testgentest.dir/build.make
testgentest: CMakeFiles/testgentest.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable testgentest"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/testgentest.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/testgentest.dir/build: testgentest
.PHONY : CMakeFiles/testgentest.dir/build

CMakeFiles/testgentest.dir/requires: CMakeFiles/testgentest.dir/test/testgentest.cc.o.requires
.PHONY : CMakeFiles/testgentest.dir/requires

CMakeFiles/testgentest.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/testgentest.dir/cmake_clean.cmake
.PHONY : CMakeFiles/testgentest.dir/clean

CMakeFiles/testgentest.dir/depend:
	cd /usr/local/google/home/drussel/extern_src/avro-cpp-1.7.2 && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /usr/local/google/home/drussel/extern_src/avro-cpp-1.7.2 /usr/local/google/home/drussel/extern_src/avro-cpp-1.7.2 /usr/local/google/home/drussel/extern_src/avro-cpp-1.7.2 /usr/local/google/home/drussel/extern_src/avro-cpp-1.7.2 /usr/local/google/home/drussel/extern_src/avro-cpp-1.7.2/CMakeFiles/testgentest.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/testgentest.dir/depend
