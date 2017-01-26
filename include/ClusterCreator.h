#ifndef CLUSTERCREATOR_H
#define CLUSTERCREATOR_H

// from std
#include <vector>

// from here
#include "Cluster.h"

typedef std::vector<Channel> Event;

class ClusterCreator {
 private:
  std::vector<Cluster> clusters;

 public:
  ClusterCreator() = default;
  ~ClusterCreator() = default;

  std::vector<Cluster> FindClustersInEventMax(const Event& dataVector,
                                              const double neighbour_threshold,
                                              const double seed_threshold,
                                              const double sum_threshold);

  std::vector<Cluster> FindClustersInEventBoole(const Event& event,
                                                const double neighbourThreshold,
                                                const double seedThreshold,
                                                const double sumThreshold,
                                                const int maxClusterSize);

  std::size_t getNumberOfClusters() const;
  const std::vector<Cluster>& getClusters() const;
};

#endif
