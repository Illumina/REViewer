cmake_minimum_required(VERSION 3.7)
project(REViewer VERSION 0.2.4 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if (NOT CMAKE_BUILD_TYPE)
    set(DEFAULT_CMAKE_BUILD_TYPE Release)
    set(CMAKE_BUILD_TYPE ${DEFAULT_CMAKE_BUILD_TYPE} CACHE STRING
            "Choose the type of build (default: ${DEFAULT_CMAKE_BUILD_TYPE})" FORCE)
endif ()

set(Boost_USE_STATIC_LIBS ON)
find_package(Boost 1.73 REQUIRED COMPONENTS program_options filesystem system)

find_package(spdlog REQUIRED)
find_package(ZLIB REQUIRED)
find_package(BZip2 REQUIRED)
find_package(LibLZMA REQUIRED)
find_package(CURL REQUIRED)
find_package(Catch2 REQUIRED)
find_package(Threads REQUIRED)
find_library(htslib libhts.a)
find_library(htslib hts)

add_subdirectory(thirdparty/graph-tools)

add_executable(REViewer
        app/REViewer.cpp
        app/Workflow.cpp app/Workflow.hh
        app/CatalogLoading.hh app/CatalogLoading.cpp
        app/ReferenceContigInfo.hh app/ReferenceContigInfo.cpp
        app/Reference.hh app/Reference.cpp
        app/GenomicRegion.hh app/GenomicRegion.cpp
        app/LocusSpecification.hh app/LocusSpecification.cpp
        app/VariantSpecification.hh app/VariantSpecification.cpp
        app/Read.hh app/Read.cpp
        app/LocusSpecDecoding.hh app/LocusSpecDecoding.cpp
        app/RegionGraph.hh app/RegionGraph.cpp
        app/GraphBlueprint.hh app/GraphBlueprint.cpp
        app/GenotypePaths.hh app/GenotypePaths.cpp
        app/Aligns.hh app/Aligns.cpp
        app/HtsHelpers.hh app/HtsHelpers.cpp
        app/Projection.hh app/Projection.cpp
        app/LanePlot.hh app/LanePlot.cpp
        app/GenerateSvg.hh app/GenerateSvg.cpp
        app/Origin.hh app/Origin.cpp
        app/Phasing.hh app/Phasing.cpp
        app/FragLenFilter.hh app/FragLenFilter.cpp)


target_include_directories(REViewer PUBLIC
        ${CMAKE_SOURCE_DIR}
        ${Boost_INCLUDE_DIRS}
        ${LIBLZMA_INCLUDE_DIRS}
        ${CURL_INCLUDE_DIRS}
        )

# Set static linking of gcc standard libraries to simplify binary distribution
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(STATIC_FLAGS -static-libgcc -static-libstdc++)
endif()

target_link_libraries(REViewer PUBLIC
        ${STATIC_FLAGS}
        graphtools
        ${htslib}
        ${Boost_LIBRARIES}
        ${LIBLZMA_LIBRARIES}
        ${CURL_LIBRARIES}
        ZLIB::ZLIB
        BZip2::BZip2
        Threads::Threads)

install(TARGETS REViewer RUNTIME DESTINATION bin)

