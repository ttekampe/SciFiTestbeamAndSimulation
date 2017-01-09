include_directories(${SciFiTestbeamAndSimulation_SOURCE_DIR}/include)

# Executables
add_executable(produceCorrectedFile produceCorrectedFile.cpp)
add_executable(clusterAnalysis clusterAnalysis.cpp)

# Linker information
target_link_libraries(produceCorrectedFile ${Boost_LIBRARIES} Calibration)
target_link_libraries(clusterAnalysis ${Boost_LIBRARIES} Calibration ClusterCreator ClusterMonitor Cluster EDouble)
