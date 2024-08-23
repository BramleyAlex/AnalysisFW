cmake_minimum_required(VERSION 3.25)

project(BDT)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(ROOT 6.16 CONFIG REQUIRED)

# We need this directory, and users of our library will need it too
# All users of this library will need at least C++17

# Add a shared module - modules are intended to be imported at runtime.
# - This is where you add the source files
add_library(TrainFactory SHARED TrainFactory.cpp)
#set_target_properties(TrainFactory PROPERTIES SUFFIX ".so")

# Set up the libraries and header search paths for this target
target_link_libraries(TrainFactory ROOT::Physics ROOT::Hist ROOT::TMVA ROOT::MathCore ROOT::Tree)

target_compile_features(TrainFactory PUBLIC cxx_std_17)

add_executable(Main main.cpp)
target_link_libraries(Main TrainFactory)