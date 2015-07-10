include_directories(${SciFiTestbeamAndSimulation_SOURCE_DIR}/include)

add_library(EDouble EDouble.cpp)
add_library(LHCbStyle lhcbStyle.cpp)
#add_library(ClusterAlgorithms ClusterAlgorithms.cpp)
add_library(Cluster Cluster.cpp)
add_library(Calibration Calibration.cpp)
add_library(ClusterCreator ClusterCreator.cpp)
add_library(ClusterMonitor ClusterMonitor.cpp)

target_link_libraries(Calibration ${ROOT_LIBRARIES} ${ROOFIT_LIBRARIES} LHCbStyle -lSpectrum)
target_link_libraries(ClusterCreator Cluster)
target_link_libraries(ClusterMonitor ${ROOT_LIBRARIES} Cluster)
