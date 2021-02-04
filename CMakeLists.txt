cmake_minimum_required(VERSION 3.7)
project(REViewer)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(ExternalProject)

set(installDir ${CMAKE_CURRENT_BINARY_DIR}/install)


ExternalProject_Add(htslib
	BUILD_IN_SOURCE YES
	GIT_REPOSITORY "https://github.com/samtools/htslib.git"
	GIT_TAG "1.10.2"
	UPDATE_COMMAND ""
	CONFIGURE_COMMAND ""
	BUILD_COMMAND make
	INSTALL_COMMAND make install prefix=${installDir}
	LOG_DOWNLOAD YES
)


ExternalProject_Add(Boost
	BUILD_IN_SOURCE YES
	URL https://dl.bintray.com/boostorg/release/1.73.0/source/boost_1_73_0.tar.bz2
	UPDATE_COMMAND ""
	CONFIGURE_COMMAND ./bootstrap.sh --prefix=${installDir}/lib
	BUILD_COMMAND ./b2 install -j8   --prefix=${installDir} --with-filesystem --with-system --with-program_options link=static,shared
	INSTALL_COMMAND ""
)


ExternalProject_Add(spdlog
	GIT_REPOSITORY "https://github.com/gabime/spdlog.git"
	GIT_TAG "v1.6.1"
	UPDATE_COMMAND ""
	CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${installDir}
)


ExternalProject_Add(catch2
	GIT_REPOSITORY "https://github.com/catchorg/Catch2.git"
	GIT_TAG "v2.12.4"
	UPDATE_COMMAND ""
	CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${installDir}
	           -DBUILD_TESTING=OFF
)


ExternalProject_Add(reviewer
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/reviewer
	BUILD_ALWAYS YES
	CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${installDir}
		   -DCMAKE_PREFIX_PATH:PATH=${installDir}
)


ExternalProject_Add_StepDependencies(reviewer configure
	Boost spdlog htslib catch2)

