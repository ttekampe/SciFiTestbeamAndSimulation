include_directories(${SciFiTestbeamAndSimulation_SOURCE_DIR}/include)

# Executables
add_executable(clusterisation clusterisation.cpp)
add_executable(gains gains.cpp)
add_executable(produceCorrectedFile produceCorrectedFile.cpp)
add_executable(clusterAnalysis clusterAnalysis.cpp)

# Linker information
target_link_libraries(clusterisation ${ROOT_LIBRARIES} cluster_lib ClusterAlgorithms Calibration)
target_link_libraries(gains Calibration LHCbStyle ${Boost_LIBRARIES})
target_link_libraries(produceCorrectedFile ${Boost_LIBRARIES} Calibration)
target_link_libraries(clusterAnalysis ${Boost_LIBRARIES} ${ROOT_LIBRARIES} Calibration ClusterAlgorithms cluster_lib)