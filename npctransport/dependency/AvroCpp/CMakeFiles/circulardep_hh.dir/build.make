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

# Utility rule file for circulardep_hh.

# Include the progress variables for this target.
include CMakeFiles/circulardep_hh.dir/progress.make

CMakeFiles/circulardep_hh: circulardep.hh

circulardep.hh: avrogencpp
circulardep.hh: jsonschemas/circulardep
	$(CMAKE_COMMAND) -E cmake_progress_report /usr/local/google/home/drussel/extern_src/avro-cpp-1.7.2/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold "Generating circulardep.hh"
	./avrogencpp -p - -i /usr/local/google/home/drussel/extern_src/avro-cpp-1.7.2/jsonschemas/circulardep -o circulardep.hh -n cd -U

circulardep_hh: CMakeFiles/circulardep_hh
circulardep_hh: circulardep.hh
circulardep_hh: CMakeFiles/circulardep_hh.dir/build.make
.PHONY : circulardep_hh

# Rule to build all files generated by this target.
CMakeFiles/circulardep_hh.dir/build: circulardep_hh
.PHONY : CMakeFiles/circulardep_hh.dir/build

CMakeFiles/circulardep_hh.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/circulardep_hh.dir/cmake_clean.cmake
.PHONY : CMakeFiles/circulardep_hh.dir/clean

CMakeFiles/circulardep_hh.dir/depend:
	cd /usr/local/google/home/drussel/extern_src/avro-cpp-1.7.2 && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /usr/local/google/home/drussel/extern_src/avro-cpp-1.7.2 /usr/local/google/home/drussel/extern_src/avro-cpp-1.7.2 /usr/local/google/home/drussel/extern_src/avro-cpp-1.7.2 /usr/local/google/home/drussel/extern_src/avro-cpp-1.7.2 /usr/local/google/home/drussel/extern_src/avro-cpp-1.7.2/CMakeFiles/circulardep_hh.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/circulardep_hh.dir/depend
