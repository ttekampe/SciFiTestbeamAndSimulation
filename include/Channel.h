#ifndef CHANNEL_H
#define CHANNEL_H

struct Channel {
  unsigned int Uplink;
  unsigned int ChannelNumber;
  double AdcValue;

  Channel() = default;
  ~Channel() = default;
  Channel(const Channel& other);             // copy construction
  Channel(Channel&& other);                  // move construction
  Channel& operator=(const Channel& other);  // copy assignment
  Channel& operator=(Channel&& other);       // move assignment
};
#endif
