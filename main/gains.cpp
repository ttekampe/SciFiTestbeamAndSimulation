//from std
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <algorithm>

//from ROOT
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TSpectrum.h"
#include "TH1D.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooMinuit.h"
#include "RooAbsReal.h"

//from here
#include "lhcbStyle.h"
#include "Calibration.h"

//from BOOST
#include "boost/program_options.hpp"

namespace po = boost::program_options;

struct config{
  unsigned int nGaussians;
  unsigned int uplinkNumber;
  std::string inputFile;
  bool debug;
};



int parseOptions(config &c, int argc, char *argv[]){

  // declare options
  po::options_description desc("Allowed options") ;
  desc.add_options()
    ("help", "show this help")
    ("ngaussians,n", po::value<unsigned int>(&c.nGaussians)->default_value(4), "number aus gaussians to fit to the distribution")
    ("inputfile,f", po::value<std::string>(&c.inputFile), "test beam data file containing the led run")
    ("uplink,u", po::value<unsigned int>(&c.uplinkNumber), "number of the uplink")
    ("debug,d", po::bool_switch(&c.debug), "debug output")
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

  lhcb::lhcbStyle();
  if(!c.debug){
    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
    RooMsgService::instance().setSilentMode(true);
  }


  const unsigned int runNumber = runNumberFromFilename(c.inputFile);
  calibrationRunNumbers calNum = lookUpCalibrationFiles(runNumber, "/data/testbeam/data/runNumbers.txt");

  std::cout << "Run number: " << runNumber << "\tdark calib: " << calNum.dark << "\tled: " << calNum.led << std::endl;

  std::string ledFileName = "btsoftware_" + std::to_string(calNum.led) + "_calib_led_ntuple.root";

  TFile inputFile( ("/data/testbeam/" + ledFileName).c_str(), "READ");
  TTree* inputTree = dynamic_cast<TTree*>(inputFile.Get("rawData"));

  produceGains(inputTree, 3, 4, 1, 128, c.nGaussians, ledFileName, "/home/tobi/SciFi/results/gains/");

  TString saveName = ledFileName;
  saveName.Remove(0, saveName.Last('/')+1);
  saveName.ReplaceAll(".root", "_gain.txt");
  std::string gainFileName("/home/tobi/SciFi/results/gains/" + saveName);

  readGains(gainFileName);
//  GetGain(inputTree, 3, 4, 1, 128, c.nGaussians, c.inputFile);

  return 0;
}
