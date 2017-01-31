#ifndef CONFIGPARSER_H
#define CONFIGPARSER_H

#include <iostream>
#include <map>
#include <string>
#include <vector>

struct FibMatInfo {
  std::string name;
  int uplink_min;
  int uplink_max;
  int adc_min{0};
  int adc_max{128};
  double x_offset;
  double z_position;
};

/**
 * \brief Read the configuration needed to perform a cluster analysis from a txt
 * file
 */
class ConfigParser {
 public:
  ConfigParser(std::string file_name);
  ~ConfigParser() {}
  // std::string operator[](std::string) const;
  std::vector<std::map<std::string, std::string>> GetConfMap() const {
    return data;
  }

  std::vector<FibMatInfo> GetConf() const;

  void print() const;

 private:
  std::vector<std::map<std::string, std::string>> data;
};

#endif
