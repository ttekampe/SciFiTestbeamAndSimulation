#ifndef CLUSTERALGORITHMS_H
#define CLUSTERALGORITHMS_H

//from here
#include "Cluster.h"

//from std
#include <iostream>
#include <fstream>
#include <vector>



void FindClustersInEvent(std::vector<Cluster*>& clusterVector,
                         const std::vector<Channel>& event,
                         const double neighbourThreshold,
                         const double seedThreshold,
                         const double sumThreshold,
                         const unsigned int maxClusterSize = {4},
                         bool debug = {false});


void FindClustersInEventBoole(std::vector<Cluster*>& clusterVector,
                                const std::vector<Channel>& event,
                                const double neighbourThreshold,
                                const double seedThreshold,
                                const double sumThreshold,
                                const unsigned int maxClusterSize = {4},
                                bool debug = {false});


#endif // CLUSTERALGORITHMS_H
