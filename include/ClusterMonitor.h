#ifndef CLUSTERMONITOR_H
#define CLUSTERMONITOR_H

// from std
#include <map>
#include <utility>
#include <vector>

// from here
#include "Cluster.h"
#include "ClusterCreator.h"

class ClusterMonitor {
 private:
 public:
  ClusterMonitor() = default;
  ~ClusterMonitor() = default;
  typedef std::map<std::string, std::vector<double> > feature_map;

  void WriteToNtuple(const ClusterCreator& clCreator,
                     const std::string fileName,
                     const feature_map features = feature_map());
};

#endif  // CLUSTER_H
