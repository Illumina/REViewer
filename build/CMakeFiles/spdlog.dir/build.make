# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.21

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.21.3/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.21.3/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/sclamons/Documents/software/REViewer

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/sclamons/Documents/software/REViewer/build

# Utility rule file for spdlog.

# Include any custom commands dependencies for this target.
include CMakeFiles/spdlog.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/spdlog.dir/progress.make

CMakeFiles/spdlog: CMakeFiles/spdlog-complete

CMakeFiles/spdlog-complete: spdlog-prefix/src/spdlog-stamp/spdlog-install
CMakeFiles/spdlog-complete: spdlog-prefix/src/spdlog-stamp/spdlog-mkdir
CMakeFiles/spdlog-complete: spdlog-prefix/src/spdlog-stamp/spdlog-download
CMakeFiles/spdlog-complete: spdlog-prefix/src/spdlog-stamp/spdlog-update
CMakeFiles/spdlog-complete: spdlog-prefix/src/spdlog-stamp/spdlog-patch
CMakeFiles/spdlog-complete: spdlog-prefix/src/spdlog-stamp/spdlog-configure
CMakeFiles/spdlog-complete: spdlog-prefix/src/spdlog-stamp/spdlog-build
CMakeFiles/spdlog-complete: spdlog-prefix/src/spdlog-stamp/spdlog-install
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/Users/sclamons/Documents/software/REViewer/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Completed 'spdlog'"
	/usr/local/Cellar/cmake/3.21.3/bin/cmake -E make_directory /Users/sclamons/Documents/software/REViewer/build/CMakeFiles
	/usr/local/Cellar/cmake/3.21.3/bin/cmake -E touch /Users/sclamons/Documents/software/REViewer/build/CMakeFiles/spdlog-complete
	/usr/local/Cellar/cmake/3.21.3/bin/cmake -E touch /Users/sclamons/Documents/software/REViewer/build/spdlog-prefix/src/spdlog-stamp/spdlog-done

spdlog-prefix/src/spdlog-stamp/spdlog-build: spdlog-prefix/src/spdlog-stamp/spdlog-configure
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/Users/sclamons/Documents/software/REViewer/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Performing build step for 'spdlog'"
	cd /Users/sclamons/Documents/software/REViewer/build/spdlog-prefix/src/spdlog-build && $(MAKE)
	cd /Users/sclamons/Documents/software/REViewer/build/spdlog-prefix/src/spdlog-build && /usr/local/Cellar/cmake/3.21.3/bin/cmake -E touch /Users/sclamons/Documents/software/REViewer/build/spdlog-prefix/src/spdlog-stamp/spdlog-build

spdlog-prefix/src/spdlog-stamp/spdlog-configure: spdlog-prefix/tmp/spdlog-cfgcmd.txt
spdlog-prefix/src/spdlog-stamp/spdlog-configure: spdlog-prefix/src/spdlog-stamp/spdlog-patch
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/Users/sclamons/Documents/software/REViewer/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Performing configure step for 'spdlog'"
	cd /Users/sclamons/Documents/software/REViewer/build/spdlog-prefix/src/spdlog-build && /usr/local/Cellar/cmake/3.21.3/bin/cmake -DCMAKE_INSTALL_PREFIX:PATH=/Users/sclamons/Documents/software/REViewer/build/install -DCMAKE_C_COMPILER=/Library/Developer/CommandLineTools/usr/bin/cc -DCMAKE_CXX_COMPILER=/Library/Developer/CommandLineTools/usr/bin/c++ "-GUnix Makefiles" /Users/sclamons/Documents/software/REViewer/build/spdlog-prefix/src/spdlog
	cd /Users/sclamons/Documents/software/REViewer/build/spdlog-prefix/src/spdlog-build && /usr/local/Cellar/cmake/3.21.3/bin/cmake -E touch /Users/sclamons/Documents/software/REViewer/build/spdlog-prefix/src/spdlog-stamp/spdlog-configure

spdlog-prefix/src/spdlog-stamp/spdlog-download: spdlog-prefix/src/spdlog-stamp/spdlog-urlinfo.txt
spdlog-prefix/src/spdlog-stamp/spdlog-download: spdlog-prefix/src/spdlog-stamp/spdlog-mkdir
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/Users/sclamons/Documents/software/REViewer/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Performing download step (download, verify and extract) for 'spdlog'"
	cd /Users/sclamons/Documents/software/REViewer/build/spdlog-prefix/src && /usr/local/Cellar/cmake/3.21.3/bin/cmake -P /Users/sclamons/Documents/software/REViewer/build/spdlog-prefix/src/spdlog-stamp/download-spdlog.cmake
	cd /Users/sclamons/Documents/software/REViewer/build/spdlog-prefix/src && /usr/local/Cellar/cmake/3.21.3/bin/cmake -P /Users/sclamons/Documents/software/REViewer/build/spdlog-prefix/src/spdlog-stamp/verify-spdlog.cmake
	cd /Users/sclamons/Documents/software/REViewer/build/spdlog-prefix/src && /usr/local/Cellar/cmake/3.21.3/bin/cmake -P /Users/sclamons/Documents/software/REViewer/build/spdlog-prefix/src/spdlog-stamp/extract-spdlog.cmake
	cd /Users/sclamons/Documents/software/REViewer/build/spdlog-prefix/src && /usr/local/Cellar/cmake/3.21.3/bin/cmake -E touch /Users/sclamons/Documents/software/REViewer/build/spdlog-prefix/src/spdlog-stamp/spdlog-download

spdlog-prefix/src/spdlog-stamp/spdlog-install: spdlog-prefix/src/spdlog-stamp/spdlog-build
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/Users/sclamons/Documents/software/REViewer/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Performing install step for 'spdlog'"
	cd /Users/sclamons/Documents/software/REViewer/build/spdlog-prefix/src/spdlog-build && $(MAKE) install
	cd /Users/sclamons/Documents/software/REViewer/build/spdlog-prefix/src/spdlog-build && /usr/local/Cellar/cmake/3.21.3/bin/cmake -E touch /Users/sclamons/Documents/software/REViewer/build/spdlog-prefix/src/spdlog-stamp/spdlog-install

spdlog-prefix/src/spdlog-stamp/spdlog-mkdir:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/Users/sclamons/Documents/software/REViewer/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Creating directories for 'spdlog'"
	/usr/local/Cellar/cmake/3.21.3/bin/cmake -E make_directory /Users/sclamons/Documents/software/REViewer/build/spdlog-prefix/src/spdlog
	/usr/local/Cellar/cmake/3.21.3/bin/cmake -E make_directory /Users/sclamons/Documents/software/REViewer/build/spdlog-prefix/src/spdlog-build
	/usr/local/Cellar/cmake/3.21.3/bin/cmake -E make_directory /Users/sclamons/Documents/software/REViewer/build/spdlog-prefix
	/usr/local/Cellar/cmake/3.21.3/bin/cmake -E make_directory /Users/sclamons/Documents/software/REViewer/build/spdlog-prefix/tmp
	/usr/local/Cellar/cmake/3.21.3/bin/cmake -E make_directory /Users/sclamons/Documents/software/REViewer/build/spdlog-prefix/src/spdlog-stamp
	/usr/local/Cellar/cmake/3.21.3/bin/cmake -E make_directory /Users/sclamons/Documents/software/REViewer/build/spdlog-prefix/src
	/usr/local/Cellar/cmake/3.21.3/bin/cmake -E make_directory /Users/sclamons/Documents/software/REViewer/build/spdlog-prefix/src/spdlog-stamp
	/usr/local/Cellar/cmake/3.21.3/bin/cmake -E touch /Users/sclamons/Documents/software/REViewer/build/spdlog-prefix/src/spdlog-stamp/spdlog-mkdir

spdlog-prefix/src/spdlog-stamp/spdlog-patch: spdlog-prefix/src/spdlog-stamp/spdlog-update
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/Users/sclamons/Documents/software/REViewer/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "No patch step for 'spdlog'"
	/usr/local/Cellar/cmake/3.21.3/bin/cmake -E echo_append
	/usr/local/Cellar/cmake/3.21.3/bin/cmake -E touch /Users/sclamons/Documents/software/REViewer/build/spdlog-prefix/src/spdlog-stamp/spdlog-patch

spdlog-prefix/src/spdlog-stamp/spdlog-update: spdlog-prefix/src/spdlog-stamp/spdlog-download
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/Users/sclamons/Documents/software/REViewer/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "No update step for 'spdlog'"
	/usr/local/Cellar/cmake/3.21.3/bin/cmake -E echo_append
	/usr/local/Cellar/cmake/3.21.3/bin/cmake -E touch /Users/sclamons/Documents/software/REViewer/build/spdlog-prefix/src/spdlog-stamp/spdlog-update

spdlog: CMakeFiles/spdlog
spdlog: CMakeFiles/spdlog-complete
spdlog: spdlog-prefix/src/spdlog-stamp/spdlog-build
spdlog: spdlog-prefix/src/spdlog-stamp/spdlog-configure
spdlog: spdlog-prefix/src/spdlog-stamp/spdlog-download
spdlog: spdlog-prefix/src/spdlog-stamp/spdlog-install
spdlog: spdlog-prefix/src/spdlog-stamp/spdlog-mkdir
spdlog: spdlog-prefix/src/spdlog-stamp/spdlog-patch
spdlog: spdlog-prefix/src/spdlog-stamp/spdlog-update
spdlog: CMakeFiles/spdlog.dir/build.make
.PHONY : spdlog

# Rule to build all files generated by this target.
CMakeFiles/spdlog.dir/build: spdlog
.PHONY : CMakeFiles/spdlog.dir/build

CMakeFiles/spdlog.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/spdlog.dir/cmake_clean.cmake
.PHONY : CMakeFiles/spdlog.dir/clean

CMakeFiles/spdlog.dir/depend:
	cd /Users/sclamons/Documents/software/REViewer/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/sclamons/Documents/software/REViewer /Users/sclamons/Documents/software/REViewer /Users/sclamons/Documents/software/REViewer/build /Users/sclamons/Documents/software/REViewer/build /Users/sclamons/Documents/software/REViewer/build/CMakeFiles/spdlog.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/spdlog.dir/depend

