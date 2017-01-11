include_directories(${SciFiTestbeamAndSimulation_SOURCE_DIR}/include)

# Executables
add_executable(produceCorrectedFile produceCorrectedFile.cpp)
add_executable(clusterAnalysis clusterAnalysis.cpp)
add_executable(clusterAnalysis2 clusterAnalysis2.cpp)
add_executable(test_cfgParser test_cfgParser.cpp)

# Linker information
target_link_libraries(produceCorrectedFile ${Boost_LIBRARIES} Calibration)
target_link_libraries(clusterAnalysis ${Boost_LIBRARIES} Calibration ClusterCreator ClusterMonitor Cluster EDouble)
target_link_libraries(clusterAnalysis2 ${Boost_LIBRARIES} Calibration ClusterCreator ClusterMonitor Cluster EDouble Alignment)
target_link_libraries(test_cfgParser ConfigParser)
