# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.5

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


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
CMAKE_COMMAND = /sw/rhel6-x64/devtools/cmake-3.5.2-gcc48/bin/cmake

# The command to remove a file.
RM = /sw/rhel6-x64/devtools/cmake-3.5.2-gcc48/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /pf/a/a270058/fesom-1.4

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /pf/a/a270058/fesom-1.4/build

# Include any dependencies generated for this target.
include fesom_cpl/CMakeFiles/fesom_C.dir/depend.make

# Include the progress variables for this target.
include fesom_cpl/CMakeFiles/fesom_C.dir/progress.make

# Include the compile flags for this target's objects.
include fesom_cpl/CMakeFiles/fesom_C.dir/flags.make

fesom_cpl/CMakeFiles/fesom_C.dir/fort_part.c.o: fesom_cpl/CMakeFiles/fesom_C.dir/flags.make
fesom_cpl/CMakeFiles/fesom_C.dir/fort_part.c.o: ../fesom_cpl/fort_part.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/pf/a/a270058/fesom-1.4/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object fesom_cpl/CMakeFiles/fesom_C.dir/fort_part.c.o"
	cd /pf/a/a270058/fesom-1.4/build/fesom_cpl && /opt/mpi/bullxmpi_mlx/1.2.8.3/bin/mpicc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/fesom_C.dir/fort_part.c.o   -c /pf/a/a270058/fesom-1.4/fesom_cpl/fort_part.c

fesom_cpl/CMakeFiles/fesom_C.dir/fort_part.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/fesom_C.dir/fort_part.c.i"
	cd /pf/a/a270058/fesom-1.4/build/fesom_cpl && /opt/mpi/bullxmpi_mlx/1.2.8.3/bin/mpicc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /pf/a/a270058/fesom-1.4/fesom_cpl/fort_part.c > CMakeFiles/fesom_C.dir/fort_part.c.i

fesom_cpl/CMakeFiles/fesom_C.dir/fort_part.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/fesom_C.dir/fort_part.c.s"
	cd /pf/a/a270058/fesom-1.4/build/fesom_cpl && /opt/mpi/bullxmpi_mlx/1.2.8.3/bin/mpicc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /pf/a/a270058/fesom-1.4/fesom_cpl/fort_part.c -o CMakeFiles/fesom_C.dir/fort_part.c.s

fesom_cpl/CMakeFiles/fesom_C.dir/fort_part.c.o.requires:

.PHONY : fesom_cpl/CMakeFiles/fesom_C.dir/fort_part.c.o.requires

fesom_cpl/CMakeFiles/fesom_C.dir/fort_part.c.o.provides: fesom_cpl/CMakeFiles/fesom_C.dir/fort_part.c.o.requires
	$(MAKE) -f fesom_cpl/CMakeFiles/fesom_C.dir/build.make fesom_cpl/CMakeFiles/fesom_C.dir/fort_part.c.o.provides.build
.PHONY : fesom_cpl/CMakeFiles/fesom_C.dir/fort_part.c.o.provides

fesom_cpl/CMakeFiles/fesom_C.dir/fort_part.c.o.provides.build: fesom_cpl/CMakeFiles/fesom_C.dir/fort_part.c.o


fesom_cpl/CMakeFiles/fesom_C.dir/psolve.c.o: fesom_cpl/CMakeFiles/fesom_C.dir/flags.make
fesom_cpl/CMakeFiles/fesom_C.dir/psolve.c.o: ../fesom_cpl/psolve.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/pf/a/a270058/fesom-1.4/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object fesom_cpl/CMakeFiles/fesom_C.dir/psolve.c.o"
	cd /pf/a/a270058/fesom-1.4/build/fesom_cpl && /opt/mpi/bullxmpi_mlx/1.2.8.3/bin/mpicc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/fesom_C.dir/psolve.c.o   -c /pf/a/a270058/fesom-1.4/fesom_cpl/psolve.c

fesom_cpl/CMakeFiles/fesom_C.dir/psolve.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/fesom_C.dir/psolve.c.i"
	cd /pf/a/a270058/fesom-1.4/build/fesom_cpl && /opt/mpi/bullxmpi_mlx/1.2.8.3/bin/mpicc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /pf/a/a270058/fesom-1.4/fesom_cpl/psolve.c > CMakeFiles/fesom_C.dir/psolve.c.i

fesom_cpl/CMakeFiles/fesom_C.dir/psolve.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/fesom_C.dir/psolve.c.s"
	cd /pf/a/a270058/fesom-1.4/build/fesom_cpl && /opt/mpi/bullxmpi_mlx/1.2.8.3/bin/mpicc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /pf/a/a270058/fesom-1.4/fesom_cpl/psolve.c -o CMakeFiles/fesom_C.dir/psolve.c.s

fesom_cpl/CMakeFiles/fesom_C.dir/psolve.c.o.requires:

.PHONY : fesom_cpl/CMakeFiles/fesom_C.dir/psolve.c.o.requires

fesom_cpl/CMakeFiles/fesom_C.dir/psolve.c.o.provides: fesom_cpl/CMakeFiles/fesom_C.dir/psolve.c.o.requires
	$(MAKE) -f fesom_cpl/CMakeFiles/fesom_C.dir/build.make fesom_cpl/CMakeFiles/fesom_C.dir/psolve.c.o.provides.build
.PHONY : fesom_cpl/CMakeFiles/fesom_C.dir/psolve.c.o.provides

fesom_cpl/CMakeFiles/fesom_C.dir/psolve.c.o.provides.build: fesom_cpl/CMakeFiles/fesom_C.dir/psolve.c.o


# Object files for target fesom_C
fesom_C_OBJECTS = \
"CMakeFiles/fesom_C.dir/fort_part.c.o" \
"CMakeFiles/fesom_C.dir/psolve.c.o"

# External object files for target fesom_C
fesom_C_EXTERNAL_OBJECTS =

fesom_cpl/libfesom_C.a: fesom_cpl/CMakeFiles/fesom_C.dir/fort_part.c.o
fesom_cpl/libfesom_C.a: fesom_cpl/CMakeFiles/fesom_C.dir/psolve.c.o
fesom_cpl/libfesom_C.a: fesom_cpl/CMakeFiles/fesom_C.dir/build.make
fesom_cpl/libfesom_C.a: fesom_cpl/CMakeFiles/fesom_C.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/pf/a/a270058/fesom-1.4/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking C static library libfesom_C.a"
	cd /pf/a/a270058/fesom-1.4/build/fesom_cpl && $(CMAKE_COMMAND) -P CMakeFiles/fesom_C.dir/cmake_clean_target.cmake
	cd /pf/a/a270058/fesom-1.4/build/fesom_cpl && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/fesom_C.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
fesom_cpl/CMakeFiles/fesom_C.dir/build: fesom_cpl/libfesom_C.a

.PHONY : fesom_cpl/CMakeFiles/fesom_C.dir/build

fesom_cpl/CMakeFiles/fesom_C.dir/requires: fesom_cpl/CMakeFiles/fesom_C.dir/fort_part.c.o.requires
fesom_cpl/CMakeFiles/fesom_C.dir/requires: fesom_cpl/CMakeFiles/fesom_C.dir/psolve.c.o.requires

.PHONY : fesom_cpl/CMakeFiles/fesom_C.dir/requires

fesom_cpl/CMakeFiles/fesom_C.dir/clean:
	cd /pf/a/a270058/fesom-1.4/build/fesom_cpl && $(CMAKE_COMMAND) -P CMakeFiles/fesom_C.dir/cmake_clean.cmake
.PHONY : fesom_cpl/CMakeFiles/fesom_C.dir/clean

fesom_cpl/CMakeFiles/fesom_C.dir/depend:
	cd /pf/a/a270058/fesom-1.4/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /pf/a/a270058/fesom-1.4 /pf/a/a270058/fesom-1.4/fesom_cpl /pf/a/a270058/fesom-1.4/build /pf/a/a270058/fesom-1.4/build/fesom_cpl /pf/a/a270058/fesom-1.4/build/fesom_cpl/CMakeFiles/fesom_C.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : fesom_cpl/CMakeFiles/fesom_C.dir/depend
