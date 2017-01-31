#ifndef CLUSTERCREATOR_H
#define CLUSTERCREATOR_H

// from std
#include <vector>

// from here
#include "Cluster.h"

typedef std::vector<Channel> Event;

/**
 * \brief Searches events for clusters and stores them
 */

class ClusterCreator {
 private:
  std::vector<Cluster> clusters;  ///< Clusters that were found

 public:
  ClusterCreator() = default;
  ~ClusterCreator() = default;

  /**
   * \brief Search an event for clusters using Max Neuners algorithm
   * \param dataVector The event (a vector of channels) that is searched
   * \param seed_threshold the adc value a channel is required to have to start
   * a new cluster
   * \param sum_threshold The minimum charge a cluster must have to be kept
   */
  std::vector<Cluster> FindClustersInEventMax(const Event& dataVector,
                                              const double neighbour_threshold,
                                              const double seed_threshold,
                                              const double sum_threshold);

  /**
  * \brief Search an event for clusters using the algorithm from Boole
  * \param event The event (a vector of channels) that is searched
  * \param neighbourThreshold The adc value an event must reach to be added to
  * the cluster
  * \param seedThreshold the adc value a channel is required to have to start
  * a new cluster
  * \param sumThreshold The minimum charge a cluster must have to be kept
  * \param maxClusterSize The maximum size a cluster is allowed to reach
  */
  std::vector<Cluster> FindClustersInEventBoole(const Event& event,
                                                const double neighbourThreshold,
                                                const double seedThreshold,
                                                const double sumThreshold,
                                                const int maxClusterSize);

  /**
   * \return The number of clusters stored in this ClusterCreator
   */
  std::size_t getNumberOfClusters() const;

  /**
   * \return The clusters stored in this ClusterCreator
   */
  const std::vector<Cluster>& getClusters() const;
};

#endif
