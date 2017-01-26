// from std
#include <iostream>
#include <utility>
// from here
#include "Channel.h"

Channel::Channel(const Channel& other)  // copy construction
    : Uplink(other.Uplink),
      ChannelNumber(other.ChannelNumber),
      AdcValue(other.AdcValue) {}

Channel::Channel(Channel&& other)  // move construction
    : Uplink(std::move(other.Uplink)),
      ChannelNumber(std::move(other.ChannelNumber)),
      AdcValue(std::move(other.AdcValue)) {}

Channel& Channel::operator=(const Channel& other) {  // copy assignment
  Uplink = other.Uplink;
  ChannelNumber = other.ChannelNumber;
  AdcValue = other.AdcValue;

  return *this;
}

Channel& Channel::operator=(Channel&& other) {  // move assignment
  Uplink = std::move(other.Uplink);
  ChannelNumber = std::move(other.ChannelNumber);
  AdcValue = std::move(other.AdcValue);
  return *this;
}
