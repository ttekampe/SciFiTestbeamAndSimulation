//from stl
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

//from ROOT
#include "TFile.h"
#include "TTree.h"

//from here
#include "Calibration.h"

//from BOOST
#include "boost/program_options.hpp"

namespace po = boost::program_options;





struct config{
  std::string file2correct;
  std::string ledFileName;
  std::string darkFileName;
  bool debug;
};




int parseOptions(config &c, int argc, char *argv[]){

  // declare options
  po::options_description desc("Allowed options") ;
  desc.add_options()
    ("help", "show this help")
    ("file2correct,f", po::value<std::string>(&c.file2correct), "test beam data file containing the led run")
    ;

  // actually do the parsing
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  // show help and exit
  if ((argc == 0) || (vm.count("help"))) {
    std::cout << desc << "\n";
    return 1;
  }

  return 0;
}

int main(int argc, char *argv[]){

  config c;
  if (parseOptions(c, argc, argv)!= 0){
    std::cerr << "Can not parse options!" << std::endl;
    return 0;
  }
  const unsigned int runNumber = runNumberFromFilename(c.file2correct);
  calibrationRunNumbers calNum = lookUpCalibrationFiles(runNumber);

  std::cout << "Run number: " << runNumber << "\tdark calib: " << calNum.dark << "\tled: " << calNum.led << std::endl;

  c.darkFileName = "btsoftware_" + std::to_string(calNum.dark) + "_calib_dark_ntuple.root";
  c.ledFileName = "btsoftware_" + std::to_string(calNum.led) + "_calib_led_ntuple.root";


  TFile darkFile(c.darkFileName.c_str(), "READ");
  TFile dataFile(c.file2correct.c_str(), "READ");
  TFile ledFile(c.ledFileName.c_str(), "READ");

  TTree* darkTree = dynamic_cast<TTree*>(darkFile.Get("rawData"));
  TTree* dataTree = dynamic_cast<TTree*>(dataFile.Get("rawData"));
  TTree* ledTree = dynamic_cast<TTree*>(ledFile.Get("rawData"));

  if(!darkTree||!ledTree||!dataTree){
    std::cerr << "could not open one of the trees!" << std::endl;
    return 0;
  }

  TString gainFileName = c.ledFileName;
  gainFileName.ReplaceAll(".root", "_gain.txt");
  std::ifstream gainFile("/home/tobi/SciFi/results/gains/" + gainFileName);
  if(!gainFile){
    produceGains(ledTree, 3, 4, 1, 128, 4, c.ledFileName, "/home/tobi/SciFi/results/gains/");
    gainFile.open("/home/tobi/SciFi/results/gains/" + gainFileName);
  }


  return 0;
}
