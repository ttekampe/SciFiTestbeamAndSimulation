include_directories(${SciFiTestbeamAndSimulation_SOURCE_DIR}/include)

# Executables
#add_executable(clusterisation clusterisation.cpp)
add_executable(gains gains.cpp)
add_executable(produceCorrectedFile produceCorrectedFile.cpp)
add_executable(clusterAnalysis clusterAnalysis.cpp)
add_executable(attScan attScan.cpp)


# Linker information
#target_link_libraries(clusterisation ${ROOT_LIBRARIES} Cluster ClusterAlgorithms Calibration)
target_link_libraries(gains Calibration LHCbStyle ${Boost_LIBRARIES})
target_link_libraries(produceCorrectedFile ${Boost_LIBRARIES} Calibration)
target_link_libraries(clusterAnalysis ${Boost_LIBRARIES} Calibration ClusterCreator ClusterMonitor Cluster EDouble)
target_link_libraries(attScan Calibration ClusterCreator Cluster LHCbStyle)
