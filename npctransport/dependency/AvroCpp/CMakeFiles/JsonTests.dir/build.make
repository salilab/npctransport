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
include CMakeFiles/JsonTests.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/JsonTests.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/JsonTests.dir/flags.make

CMakeFiles/JsonTests.dir/test/JsonTests.cc.o: CMakeFiles/JsonTests.dir/flags.make
CMakeFiles/JsonTests.dir/test/JsonTests.cc.o: test/JsonTests.cc
	$(CMAKE_COMMAND) -E cmake_progress_report /usr/local/google/home/drussel/extern_src/avro-cpp-1.7.2/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/JsonTests.dir/test/JsonTests.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/JsonTests.dir/test/JsonTests.cc.o -c /usr/local/google/home/drussel/extern_src/avro-cpp-1.7.2/test/JsonTests.cc

CMakeFiles/JsonTests.dir/test/JsonTests.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/JsonTests.dir/test/JsonTests.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /usr/local/google/home/drussel/extern_src/avro-cpp-1.7.2/test/JsonTests.cc > CMakeFiles/JsonTests.dir/test/JsonTests.cc.i

CMakeFiles/JsonTests.dir/test/JsonTests.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/JsonTests.dir/test/JsonTests.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /usr/local/google/home/drussel/extern_src/avro-cpp-1.7.2/test/JsonTests.cc -o CMakeFiles/JsonTests.dir/test/JsonTests.cc.s

CMakeFiles/JsonTests.dir/test/JsonTests.cc.o.requires:
.PHONY : CMakeFiles/JsonTests.dir/test/JsonTests.cc.o.requires

CMakeFiles/JsonTests.dir/test/JsonTests.cc.o.provides: CMakeFiles/JsonTests.dir/test/JsonTests.cc.o.requires
	$(MAKE) -f CMakeFiles/JsonTests.dir/build.make CMakeFiles/JsonTests.dir/test/JsonTests.cc.o.provides.build
.PHONY : CMakeFiles/JsonTests.dir/test/JsonTests.cc.o.provides

CMakeFiles/JsonTests.dir/test/JsonTests.cc.o.provides.build: CMakeFiles/JsonTests.dir/test/JsonTests.cc.o

# Object files for target JsonTests
JsonTests_OBJECTS = \
"CMakeFiles/JsonTests.dir/test/JsonTests.cc.o"

# External object files for target JsonTests
JsonTests_EXTERNAL_OBJECTS =

JsonTests: CMakeFiles/JsonTests.dir/test/JsonTests.cc.o
JsonTests: libavrocpp.so.1.7.2.0
JsonTests: /usr/lib/libboost_filesystem-mt.so
JsonTests: /usr/lib/libboost_system-mt.so
JsonTests: /usr/lib/libboost_program_options-mt.so
JsonTests: CMakeFiles/JsonTests.dir/build.make
JsonTests: CMakeFiles/JsonTests.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable JsonTests"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/JsonTests.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/JsonTests.dir/build: JsonTests
.PHONY : CMakeFiles/JsonTests.dir/build

CMakeFiles/JsonTests.dir/requires: CMakeFiles/JsonTests.dir/test/JsonTests.cc.o.requires
.PHONY : CMakeFiles/JsonTests.dir/requires

CMakeFiles/JsonTests.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/JsonTests.dir/cmake_clean.cmake
.PHONY : CMakeFiles/JsonTests.dir/clean

CMakeFiles/JsonTests.dir/depend:
	cd /usr/local/google/home/drussel/extern_src/avro-cpp-1.7.2 && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /usr/local/google/home/drussel/extern_src/avro-cpp-1.7.2 /usr/local/google/home/drussel/extern_src/avro-cpp-1.7.2 /usr/local/google/home/drussel/extern_src/avro-cpp-1.7.2 /usr/local/google/home/drussel/extern_src/avro-cpp-1.7.2 /usr/local/google/home/drussel/extern_src/avro-cpp-1.7.2/CMakeFiles/JsonTests.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/JsonTests.dir/depend
