// from std
#include <cmath>
#include <iostream>
#include <vector>

// from ROOT
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TMath.h"
#include "TLeaf.h"
// from here
#include "Calibration.h"
#include "EDouble.h"

// from BOOST
#include "boost/filesystem.hpp"
#include "boost/program_options.hpp"

namespace po = boost::program_options;
struct config {
  std::string file2analyze;
  unsigned int chanMin;
  unsigned int chanMax;
  std::string outputDir;
};

int parseOptions(config &c, int argc, char *argv[]) {

  // declare options
  po::options_description desc("Allowed options");
  desc.add_options()("help,h", "show this help")
    ("file,f",po::value<std::string>(&c.file2analyze),"corrected test beam data file")
    ("chanMin,l",po::value<unsigned int>(&c.chanMin)->default_value(1),"lowest channel (default 1)")
    ("chanMax,u",po::value<unsigned int>(&c.chanMax)->default_value(128),"highest channel (default 128)")
    ("odir,o",po::value<std::string>(&c.outputDir)->default_value(std::string(std::getenv("HOME")) + "/SciFi"),"directory the output is written to");

  // actually do the parsing
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  // show help and exit
  if ((argc == 1) || (vm.count("help"))) {
    std::cout << desc << "\n";
    return 1;
  }

  return 0;
}


void drawProfile(TTree* tree, const unsigned int chanMin, const unsigned int chanMax,TString filename){
  /*
    format of the ntuple:
    Uplink_1_adc_1  
    Uplink_1_adc_2  
    Uplink_1_adc_3  
   */
  std::vector<TH1F> hists;
  //init histograms
  for(unsigned int uplink = 1; uplink<=4;++uplink){//loop over uplinks
    int u1 = uplink+(uplink-1);
    int u2 = u1+1;
    
    TH1F tmp(Form("BeamProfileUplinks_%d_%d",u1,u2),
	     Form("BeamProfileUplinks_%d_%d",u1,u2),
	     2*chanMax-chanMin+1,chanMin,2*chanMax);//both uplinks in one hist
    hists.push_back(tmp);
  }
  //fill histograms
  for(unsigned int uplink = 1; uplink <=4; ++uplink){
    int u1 = uplink+(uplink-1);
    int u2 = u1+1;    
    for(int channel = chanMin; channel <=chanMax;++channel){//loop over channels per uplink
      double u1val = tree->GetLeaf(Form("Uplink_%d_adc_%d",u1,channel))->GetValue();
      hists.at(uplink-1).SetBinContent(channel,u1val);
      double u2val = tree->GetLeaf(Form("Uplink_%d_adc_%d",u2,channel))->GetValue();
      hists.at(uplink-1).SetBinContent(2*channel,u2val);
    }//end loop on channels
  }//end loop on uplinks
  TFile f1(filename,"RECREATE");
  f1.cd();
  for(auto hist: hists){
    hist.Write();
  }
  f1.Close();
}

int main(int argc, char *argv[]) {
  std::cout<<"Trivial script to plot the beam profile for the given file"<<std::endl;
  config c;
  if (parseOptions(c, argc, argv) != 0) {
    std::cerr << "Can not parse options!" << std::endl;
    return 0;
  }

  TString savename(c.outputDir);
  savename+=Form("%d_beam_profile.root",  runNumberFromFilename(removePath(c.file2analyze).Data()));
  std::cout<<"will draw beam profile and save to "<<savename.Data()<<std::endl;
  //drawProfile(tree,c.chanMin,c.chanMax,savename);
  
  return 0;
}
