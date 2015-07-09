#ifndef CLUSTERCREATOR_H
#define CLUSTERCREATOR_H

//from std
#include <vector>

//from here
#include "Cluster.h"

typedef std::vector<Channel> Event;


class ClusterCreator{
private:
  std::vector<Cluster*> clusters;

public:
  ClusterCreator(){};
  ~ClusterCreator();
  void FindClustersInEventMax(
    const Event& dataVector
    ,const double neighbour_threshold
    ,const double seed_threshold
    ,const double sum_threshold
    );

  void FindClustersInEventBoole(
    const Event& event
    ,const double neighbourThreshold
    ,const double seedThreshold
    ,const double sumThreshold
    ,const int maxClusterSize
    ,bool debug
    );

  int getNumberOfClusters() {return clusters.size();}
  const std::vector<Cluster*> &getClusters() const {return clusters;}

};

#endif
