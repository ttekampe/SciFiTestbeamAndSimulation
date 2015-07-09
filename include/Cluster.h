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
  Cluster(){}
  ~Cluster(){}
  void Clear() {RelatedChannels.clear();}
  void AddChannel(const unsigned int uplNumber, const unsigned int chanNumber, const double adcValue);
  void RemoveChannel(const unsigned int chanNumber);
  void Resize(const unsigned int newsize);
  double GetChargeWeightedMean() const;
  double GetHitWeightedMean() const;
  double GetSumOfAdcValues() const;
  double GetMaximumAdcValue() const;
  unsigned int GetClusterSize() const {return RelatedChannels.size();}
  unsigned int GetMinChannel() const;
  unsigned int GetMaxChannel() const;
  unsigned int GetSeedChannelNumber() const;

  const std::vector<Channel> &GetRelatedChannels() const {return RelatedChannels;}

};

#endif // CLUSTER_H
