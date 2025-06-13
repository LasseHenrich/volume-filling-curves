cmake_minimum_required (VERSION 3.29.0)

project ("volume-filling-curves")

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  # using Visual Studio C++
  message("Using Visual Studio compiler flags")
  set(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} /W4")
  set(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} /MP") # parallel build
  SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} /wd\"4267\"")  # ignore conversion to smaller type (fires more aggressively than the gcc version, which is annoying)
  SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} /wd\"4244\"")  # ignore conversion to smaller type (fires more aggressively than the gcc version, which is annoying)
  SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} /wd\"4305\"")  # ignore truncation on initialization
  SET(CMAKE_CXX_FLAGS "${BASE_CXX_FLAGS} ${DISABLED_WARNINGS}")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MD")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MDd")

  add_definitions(/D "_CRT_SECURE_NO_WARNINGS")
  add_definitions(-DNOMINMAX)
  add_definitions(-D_USE_MATH_DEFINES)
else()
  # unrecognized
  message( FATAL_ERROR "Unrecognized compiler [${CMAKE_CXX_COMPILER_ID}]" )
endif()

set(VCPKG_ROOT "C:/Users/lasse/Documents/GitHub/vcpkg")

find_package(OpenVDB REQUIRED CONFIG)

# Deps
add_subdirectory(deps/geometry-central)
add_subdirectory(deps/polyscope)

# read the modules
add_subdirectory(modules)

# tinyad
include(FetchContent)
FetchContent_Declare(
  tinyad
  GIT_REPOSITORY https://github.com/patr-schm/tinyad.git
  GIT_TAG 75093e14ef0d7bb39657c5f3b2aba1251afaa38c
)

FetchContent_GetProperties(tinyad)
if(NOT tinyad_POPULATED)
  # Fetch the content using previously declared details
  FetchContent_Populate(tinyad)
  message(STATUS "tinyad_SOURCE_DIR: ${tinyad_SOURCE_DIR}")
  message(STATUS "tinyad_BINARY_DIR: ${tinyad_BINARY_DIR}")
  add_subdirectory(${tinyad_SOURCE_DIR} ${tinyad_BINARY_DIR})
endif()

# Add source to this project's executable.
add_executable (volume-filling-curves "volume-filling-curves.cpp")

target_include_directories(volume-filling-curves PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/modules/include)

# add the args.hxx project which we use for command line args
target_include_directories(volume-filling-curves PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/deps/polyscope/deps/args")
target_include_directories(volume-filling-curves PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/deps/libigl/include")
target_link_libraries(volume-filling-curves PUBLIC geometry-central polyscope modules TinyAD PRIVATE OpenVDB::openvdb)

target_include_directories(volume-filling-curves PRIVATE
    "${VCPKG_ROOT}/installed/x64-windows/include"
)

# tinyad needs C++17 (not sure if that means this project does too)
target_compile_features(volume-filling-curves PRIVATE cxx_std_17)
