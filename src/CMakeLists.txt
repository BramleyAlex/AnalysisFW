find_package(ROOT 6.16 CONFIG REQUIRED)

# We need this directory, and users of our library will need it too
# All users of this library will need at least C++17

# Add a shared module - modules are intended to be imported at runtime.
# - This is where you add the source files
add_library(Analysis SHARED CLoop.cpp CLoopBackend.cpp Analysis.cpp Tools.cpp OutputTree.cpp HistogramContainer.cpp Analysis_HM.cpp Analysis_BDT.cpp)
set_target_properties(Analysis PROPERTIES SUFFIX ".so")

# Set up the libraries and header search paths for this target
target_link_libraries(Analysis ROOT::Physics ROOT::Hist ROOT::TMVA ROOT::MathCore ROOT::Tree)
target_include_directories(Analysis PUBLIC ../include)

target_compile_features(Analysis PUBLIC cxx_std_17)

install(TARGETS Analysis DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)