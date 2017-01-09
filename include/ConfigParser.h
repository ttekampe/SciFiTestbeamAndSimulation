#include <iostream>
#include <map>
#include <string>
#include <vector>

class ConfigParser {
public:
  ConfigParser(std::string file_name);
  ~ConfigParser() {}
  // std::string operator[](std::string) const;
  std::vector<std::map<std::string, std::string>> GetConfMap() const {
    return data;
  }

private:
  std::vector<std::map<std::string, std::string>> data;
};
