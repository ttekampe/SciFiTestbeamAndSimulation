include_directories(${SciFiTestbeamAndSimulation_SOURCE_DIR}/include)

add_library(EDouble EDouble.cpp)
add_library(LHCbStyle lhcbStyle.cpp)
add_library(ClusterAlgorithms ClusterAlgorithms.cpp)
add_library(cluster_lib Cluster.cpp)