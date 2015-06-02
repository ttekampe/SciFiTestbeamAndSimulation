#include "Cluster.h"

#include <iostream>
#include <vector>

Cluster::Cluster(){}

Cluster::~Cluster(){}

void Cluster::Clear(){
  RelatedChannels.clear();
}

void Cluster::AddChannel(const unsigned int uplNumber, const unsigned int chanNumber, const double adcValue){
  Channel c;
  c.Uplink = uplNumber;
  c.ChannelNumber = chanNumber;
  c.AdcValue = adcValue;
  RelatedChannels.push_back(c);
}

void Cluster::RemoveChannel(const unsigned int chanNumber){
  for(std::vector<Channel>::iterator chan = RelatedChannels.begin(); chan != RelatedChannels.end(); ++chan){
    if(chan->ChannelNumber == chanNumber){
      RelatedChannels.erase(chan);
      return;
    }
  }
  std::cerr << "Cluster::RemoveChannel: Channel number " << chanNumber << " not in cluster -> cannot remove." << std::endl;
}

void Cluster::Resize(const unsigned int newsize){
  while(RelatedChannels.size() > newsize){
    std::cout << "current cluster size is " << RelatedChannels.size() << " max size is " << newsize << " -> reducing" << std::endl;
    if(RelatedChannels.front().AdcValue > RelatedChannels.back().AdcValue) RelatedChannels.pop_back();
    else RelatedChannels.erase(RelatedChannels.begin());
  }
}

double Cluster::GetChargeWeightedMean(){
  double totalCharge{0.};
  double mean{0.};
  for(const auto& chan : RelatedChannels){
    mean += chan.ChannelNumber * chan.AdcValue;
    totalCharge += chan.AdcValue;
  }
  return mean/totalCharge;
}

double Cluster::GetHitWeightedMean(){
  double mean{0.};
  for(const auto& chan : RelatedChannels){
    mean += chan.ChannelNumber;
  }
  return mean/RelatedChannels.size();
}

double Cluster::GetSumOfAdcValues(){
  double sumOfAdcValues{0.};
    for(const auto& chan : RelatedChannels){
    sumOfAdcValues += chan.AdcValue;
  }
  return sumOfAdcValues;
}

double Cluster::GetMaximumAdcValue(){
  double maxAdcValue{0.};
  for(const auto& chan : RelatedChannels){
    if(chan.AdcValue > maxAdcValue) maxAdcValue = chan.AdcValue;
  }
  return maxAdcValue;
}

unsigned int Cluster::GetClusterSize(){
  return RelatedChannels.size();
}

unsigned int Cluster::GetMinChannel(){
  unsigned int minChannel{9999999};
  for(const auto& chan : RelatedChannels){
    if(chan.ChannelNumber < minChannel) minChannel = chan.ChannelNumber;
  }
  return minChannel;
}
unsigned int Cluster::GetMaxChannel(){
  unsigned int maxChannel{0};
  for(const auto& chan : RelatedChannels){
    if(chan.ChannelNumber > maxChannel) maxChannel = chan.ChannelNumber;
  }
  return maxChannel;
}
