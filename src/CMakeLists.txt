include_directories(${SciFiTestbeamAndSimulation_SOURCE_DIR}/include)

add_library(Calibration Calibration.cpp)
add_library(Channel Channel.cpp)
add_library(Cluster Cluster.cpp)
add_library(ClusterCreator ClusterCreator.cpp)
add_library(ClusterMonitor ClusterMonitor.cpp)
add_library(ConfigParser ConfigParser.cpp)
add_library(EDouble EDouble.cpp)
add_library(LHCbStyle lhcbStyle.cpp)


target_link_libraries(Calibration Channel ${ROOT_LIBRARIES} ${ROOFIT_LIBRARIES} LHCbStyle -lSpectrum)
target_link_libraries(Cluster Channel)
target_link_libraries(ClusterMonitor ${ROOT_LIBRARIES} Cluster)
target_link_libraries(ClusterCreator Cluster)
