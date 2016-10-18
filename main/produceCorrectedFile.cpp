//from stl
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <unistd.h>//cwd
//from ROOT
#include "TFile.h"
#include "TTree.h"

//from here
#include "Calibration.h"

//from BOOST
#include "boost/program_options.hpp"
namespace po = boost::program_options;
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>

struct config{
  std::string file2correct;
  std::string ledFileName;
  std::string darkFileName;
  unsigned int uplinkMin;
  unsigned int uplinkMax;
  bool debug;
  std::string outputLocation;//AD 18-10-16
};
//todo. Change all output locations to the config struct val
//todo. make directories if they don't exist already
//todo. 
//AD 18-10-16. Check for existence of file. See http://stackoverflow.com/questions/12774207/fastest-way-to-check-if-a-file-exist-using-standard-c-c11-c
inline bool file_exists (const std::string& name) {
  std::ifstream f(name.c_str());
  return f.good();
}

int parseOptions(config &c, int argc, char *argv[]){

  // declare options
  po::options_description desc("Allowed options") ;
  desc.add_options()
    ("help,h", "show this help")
    ("file2correct,f", po::value<std::string>(&c.file2correct), "test beam data file")
    ("umax,u", po::value<unsigned int>(&c.uplinkMax)->default_value(4), "uplink number to start at [1, ... , 8]")
    ("umin,l", po::value<unsigned int>(&c.uplinkMin)->default_value(3), "uplink number to stop at [1, ... , 8]")
    ("outputlocation,o",po::value<std::string>(&c.outputLocation), "output location")
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
  if(c.uplinkMax<c.uplinkMin){
      std::cerr << "error in configuration: uplink min is bigger that uplink max\n";
      return  0;
  }
  //fix issue with multiple underscores. AD 18-10-2016
  const unsigned int runNumber = runNumberFromFilename(removePath(c.file2correct).Data());
  //fix issue with search for runNumbers.txt. AD 18-10-2016
  std::string exec_path = boost::filesystem::system_complete(argv[0]).string();
  std::string path_head = exec_path.substr(0, exec_path.find("SciFiTestbeamAndSimulation"));
  std::string location =  path_head+"SciFiTestbeamAndSimulation/RunNumbers/runNumbers.txt";
  std::cout<<"looking for location"<<location<<std::endl;
  //check if file exists
  if(!file_exists(location)){
    std::cout<<"Something terribly wrong with finding the run numbers"<<std::endl;
    return 0;
  }
  calibrationRunNumbers calNum = lookUpCalibrationFiles(runNumber, location);
  
  std::cout << "Run number: " << runNumber << "\tdark calib: " << calNum.dark << "\tled: " << calNum.led << std::endl;

  std::string data_path = boost::filesystem::system_complete(c.file2correct).string();
  data_path = data_path.substr(0,data_path.find(removePath(c.file2correct).Data()));
  std::cout << "Finding data at " << data_path << std::endl;
  
  c.darkFileName = data_path + "btsoftware_" + std::to_string(calNum.dark) + "_calib_dark_ntuple.root";
  c.ledFileName =  data_path + "/btsoftware_" + std::to_string(calNum.led) + "_calib_led_ntuple.root";

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
  gainFileName = removePath(gainFileName);
  gainFileName.ReplaceAll(".root", "_gain.txt");
  gainFileName = "/home/tobi/SciFi/results/gains/" + gainFileName;


  std::cout << " looking for file " << gainFileName << std::endl;
  std::ifstream gainFile(gainFileName);
  if(!gainFile){
    std::cout << "could not find gainfile -> producing gains" << std::endl;
    produceGains(ledTree, c.uplinkMin, c.uplinkMax, 1, 128, 4, c.ledFileName, "/home/tobi/SciFi/results/gains/");
  }
  else std::cout << "found gain file." << std::endl;
  gainFile.close();

  std::cout << "reading gains" << std::endl;
  std::map<unsigned int, std::map<unsigned int, double>> gains = readGains(gainFileName.Data());


  std::cout << "reading pedestals" << std::endl;
  std::map<unsigned int, std::map<unsigned int, double> > pedestals = getPedestals(c.darkFileName, c.uplinkMin, c.uplinkMax, 128);

  for(unsigned int i=c.uplinkMin; i<=c.uplinkMax; ++i){
    for(unsigned int j = 1; j<=128; ++j){
      std::cout << "uplink: " << i << "\tadc: " << j << "\t pesdestal: " << pedestals[i][j]<< "\t gain: " << gains[i][j] << std::endl;
    }
  }


  TString newFileName = "/corrected/" + removePath(c.file2correct);
  newFileName.ReplaceAll(".root", "_corrected.root");


//  TString newFileName = "/data/old_testbeam/corrected" + c.file2correct;
//  newFileName.ReplaceAll(".root", "_corrected.root");


  std::cout << "correcting file" << std::endl;
  correctFile(dataTree, gains, pedestals, c.uplinkMin, c.uplinkMax, 128, newFileName);

  return 0;
}
