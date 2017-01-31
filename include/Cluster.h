#ifndef CLUSTER_H
#define CLUSTER_H

// from std
#include <vector>
// from here
#include "Channel.h"

/**
 * \brief A cluster
 */
class Cluster {
 private:
  std::vector<Channel>
      RelatedChannels;  //< All adc channels that belong to this cluster

 public:
  Cluster() = default;
  ~Cluster() = default;
  //\brief copy construction
  Cluster(const Cluster& other);
  //\brief move construction
  Cluster(Cluster&& other);
  //\brief copy assignment
  Cluster& operator=(const Cluster& other);
  //\brief move assignment
  Cluster& operator=(Cluster&& other);

  ///\brief clears the list of channels belonging to this cluster
  void Clear();

  /**
   * \brief adds a channel to this cluster
   * \param uplNumber number of the uplink
   * \param chanNumber number of the adc within the uplink
   * \param adcValue adc value
   */
  void AddChannel(const unsigned int uplNumber, const unsigned int chanNumber,
                  const double adcValue);
  /**
   * \brief removes a channel from a cluster
   * \param chanNumber number (within the uplink) of the channel that shall be
   * removed
   */
  void RemoveChannel(const unsigned int chanNumber);

  /**
   * \brief resizes the cluster by removing the rightmost or leftmost channel,
   * depending on which adcValue is smaller, until the given size is reached
   * \param newsize cluster size after this function
   */
  void Resize(const unsigned int newsize);

  /**
   * \brief Calculates the mean cluster position weighted by the adcValue
   * \return charge weighted cluster mean in units of channel number
   */
  double GetChargeWeightedMean() const;

  /**
   * \brief Calculates the mean cluster position
   * \return cluster mean in units of channel number
   */
  double GetHitWeightedMean() const;

  /**
   * \brief Calculates the sum of adc values (collected charge) of a cluster
   * \return sum of adc values
   */
  double GetSumOfAdcValues() const;

  /**
   * \return The largest adc value in the clsuter
   */
  double GetMaximumAdcValue() const;

  /**
   * \return The size of the cluster
   */
  unsigned int GetClusterSize() const;

  /**
   * \return The number of the channel with smallest channel number (where the
   * cluster begins)
   */
  unsigned int GetMinChannel() const;

  /**
   * \return The number of the channel with highest channel number (where the
   * cluster ends)
   */
  unsigned int GetMaxChannel() const;

  /**
   * \return The number of the channel with highest adc value
   */
  unsigned int GetSeedChannelNumber() const;

  /**
   * Retruns the channels that belong to the cluster
   */
  const std::vector<Channel>& GetRelatedChannels() const;
};

#endif  // CLUSTER_H
