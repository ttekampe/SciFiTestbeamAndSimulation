#ifndef CLUSTERMONITOR_H
#define CLUSTERMONITOR_H

// from std
#include <map>
#include <utility>
#include <vector>

// from here
#include "Cluster.h"
#include "ClusterCreator.h"

/**
 * \brief Stores cluster information into ntuples
 */

class ClusterMonitor {
 private:
 public:
  ClusterMonitor() = default;
  ~ClusterMonitor() = default;

  typedef std::map<std::string, std::vector<double> > feature_map;
  ///< A TTree is created based on this map, where the key is the branch name
  /// and the value a vector with the value for each event

  /**
   * \brief Stores the clusters collected by a ClusterCreator into a TTree
   * \param clCreator The ClusterCreator whose results shall be stored
   * \param fileName  The file in which the results are stored
   * \param features  Additional branches that shall go into the TTree, for
   * example track distances
   */
  void WriteToNtuple(const ClusterCreator& clCreator,
                     const std::string fileName,
                     const feature_map features = feature_map());
};

#endif  // CLUSTER_H
