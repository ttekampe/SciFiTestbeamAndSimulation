include_directories(${SciFiTestbeamAndSimulation_SOURCE_DIR}/include)

# Executables
add_executable(clusterisation clusterisation.cpp)
add_executable(gains gains.cpp)
add_executable(produceCorrectedFile produceCorrectedFile.cpp)


# Linker information
target_link_libraries(clusterisation ${ROOT_LIBRARIES} cluster_lib ClusterAlgorithms Calibration)
target_link_libraries(gains ${ROOT_LIBRARIES} ${ROOFIT_LIBRARIES} LHCbStyle ${Boost_LIBRARIES} -lSpectrum)
target_link_libraries(produceCorrectedFile ${Boost_LIBRARIES} Calibration)