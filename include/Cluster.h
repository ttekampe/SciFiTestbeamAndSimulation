#ifndef CLUSTER_H
#define CLUSTER_H

#include <vector>

struct Channel{
  unsigned int Uplink;
  unsigned int ChannelNumber;
  double AdcValue;
};

class Cluster{
  private:
  std::vector<Channel> RelatedChannels;

  public:
  Cluster();
  ~Cluster();
  void Clear();
  void AddChannel(const unsigned int uplNumber, const unsigned int chanNumber, const double adcValue);
  void RemoveChannel(const unsigned int chanNumber);
  void Resize(const unsigned int newsize);
  double GetChargeWeightedMean();
  double GetHitWeightedMean();
  double GetSumOfAdcValues();
  double GetMaximumAdcValue();
  unsigned int GetClusterSize();
  unsigned int GetMinChannel();
  unsigned int GetMaxChannel();

};

#endif // CLUSTER_H
