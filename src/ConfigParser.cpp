#include "ConfigParser.h"
#include <exception>
#include <fstream>
#include <sstream>

ConfigParser::ConfigParser(const std::string file_name) {
  std::ifstream cfg_file(file_name);
  std::map<std::string, std::string> current_conf;
  for (std::string line; std::getline(cfg_file, line);) {

    // ignore lines that are marked as comments
    if (line[0] == '#') {
      continue;
    }

    std::string prefix("name");
    if (line.substr(0, prefix.size()) == prefix) {
      // begin a new fibre mat if a new name is set
      if (!current_conf.empty()) {
        data.push_back(current_conf);
      }
      current_conf.clear();
    }

    std::istringstream iss(line);
    std::string key, eq, val;

    if (!(iss >> key >> eq >> val >> std::ws) || eq != "=" ||
        iss.get() != EOF) {
      std::cerr << "Warning: The following line of your cfg file could not be "
                   "interpreted and will be ignored:\n";
      std::cerr << line << "\n";
    }

    current_conf[key] = val;
  }
  if (!current_conf.empty()) {
    data.push_back(current_conf);
  }
}

std::vector<FibMatInfo> ConfigParser::GetConf() const {
  std::vector<FibMatInfo> rv;
  for (const auto &mat : data) {
    FibMatInfo current_mat;
    for (const auto &c : mat) {
      if (c.first == "name") {
        current_mat.name = c.second;
      } else if (c.first == "uplink_min") {
        current_mat.uplink_min = std::stoi(c.second);
      } else if (c.first == "uplink_max") {
        current_mat.uplink_max = std::stoi(c.second);
      } else if (c.first == "adc_min") {
        current_mat.adc_min = std::stoi(c.second);
      } else if (c.first == "adc_max") {
        current_mat.adc_max = std::stoi(c.second);
      } else if (c.first == "z_position") {
        current_mat.z_position = std::stod(c.second);
      } else {
        std::cerr << "Unknown option on config file:\t" << c.first << " = "
                  << c.second << "\n";
      }
    }
    rv.push_back(current_mat);
  }
  return rv;
}

void ConfigParser::print() const {
  std::cout << data.size() << " fibre mats are set up\n";
  std::cout << "Data is configured as follows:\n";
  for (const auto &cfgMap : data) {
    for (const auto &element : cfgMap) {
      std::cout << element.first << "\t" << element.second << "\n";
    }
  }
}

// std::string ConfigParser::operator[][](const unsigned int idx,
//                                        const std::string key) const {
//   if (data[idx].find("f") == data.end()) {
//     // the key does not exist in the map
//     std::cerr << "ConfigParser::operator[] the key you are looking for
//     does "
//                  "not exist!\n";
//     std::cerr << key << "\n";
//     throw std::invalid_argument("Key is not present in config");
//   } else {
//     return data[idx].at(key);
//   }
// }
