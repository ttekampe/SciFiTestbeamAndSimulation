#ifndef CHANNEL_H
#define CHANNEL_H

/**
 * \brief An ADC channel
 */
struct Channel {
  unsigned int Uplink;         //< Number of the uplink the adc belongs to
  unsigned int ChannelNumber;  //< Number of the adc within its uplink
  double AdcValue;             //< Value of this adc channel

  Channel() = default;
  ~Channel() = default;
  //\brief copy construction
  Channel(const Channel& other);

  //\bief move construction
  Channel(Channel&& other);

  //\bief copy assignment
  Channel& operator=(const Channel& other);

  //\bief move assignment
  Channel& operator=(Channel&& other);
};
#endif
