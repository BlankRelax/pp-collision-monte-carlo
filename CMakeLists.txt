cmake_minimum_required(VERSION 3.14)
project(ppmc VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release)
endif()

if(MSVC)
  add_compile_options(/W4)
else()
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# The physics, as a library shared by the program and the tests.
add_library(ppmc_core
  src/PDF.cpp
  src/DrellYan.cpp
  src/QCDDijet.cpp
  src/Histogram.cpp)
target_include_directories(ppmc_core PUBLIC include)

# The command-line event generator.
add_executable(ppmc src/main.cpp)
target_link_libraries(ppmc PRIVATE ppmc_core)

# Physics self-tests: run with `ctest` from the build directory.
enable_testing()
add_executable(test_physics tests/test_physics.cpp)
target_link_libraries(test_physics PRIVATE ppmc_core)
add_test(NAME physics COMMAND test_physics)
