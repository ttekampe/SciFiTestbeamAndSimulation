#ifndef CLUSTERMONITOR_H
#define CLUSTERMONITOR_H

//from std
#include <vector>
#include <utility>

//from here
#include "Cluster.h"
#include "ClusterCreator.h"


class ClusterMonitor{

private:

public:
  ClusterMonitor(){};
  ~ClusterMonitor(){};

  void WriteToNtuple(const ClusterCreator& clCreator, std::string fileName);


};


#endif // CLUSTER_H
