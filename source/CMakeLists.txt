cmake_minimum_required(VERSION 3.10)

project(nssi)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# This project is heavily based on the OpenMP library
find_package(OpenMP REQUIRED)

include_directories(include include/nlohmann include/Eigen)

# All the libraries are header-only, thus
# there's no necessasity to include something 
# else but my own source cxx files
file(GLOB sources "src/*.cpp")
add_executable(nssi ${sources})

# Link OpenMP
target_link_libraries(nssi PRIVATE OpenMP::OpenMP_CXX)