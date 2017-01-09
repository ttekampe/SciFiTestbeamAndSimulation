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
      data.push_back(current_conf);
      current_conf.clear();
    }

    std::istringstream iss(line);
    std::string key, eq, val;

    if (!(iss >> key >> eq >> val >> std::ws) || eq != "=" ||
        iss.get() != EOF) {
      std::cerr << "Warning: The following line of you cfg file could not be "
                   "interpreted:\n";
      std::cerr << line << "\n";
    }

    current_conf[key] = val;
  }
  if (!current_conf.empty()) {
    data.push_back(current_conf);
  }
}

// std::string ConfigParser::operator[][](const unsigned int idx,
//                                        const std::string key) const {
//   if (data[idx].find("f") == data.end()) {
//     // the key does not exist in the map
//     std::cerr << "ConfigParser::operator[] the key you are looking for does "
//                  "not exist!\n";
//     std::cerr << key << "\n";
//     throw std::invalid_argument("Key is not present in config");
//   } else {
//     return data[idx].at(key);
//   }
// }
