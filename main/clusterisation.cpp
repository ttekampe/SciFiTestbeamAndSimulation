//from here
#include "Cluster.h"
#include "ClusterAlgorithms.h"
#include "Calibration.h"

//from std
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

//from ROOT
#include "TFile.h"
#include "TTree.h"
#include "TLeafObject.h"

struct Config{
  std::string FileName;
  std::string PedestalFileName;
  std::string GainFileName;
  bool debug;
};



int main(){
  Config cfg;

cfg.FileName = "/data/old_testbeam/corrected/btsoftware_1414634337_datarun_ntuple_corrected.root";
//  cfg.FileName = "/data/old_testbeam/btsoftware_1414634337_datarun_ntuple.root";
  cfg.GainFileName = "/home/tobi/SciFi/results/gains/btsoftware_1414634319_calib_led_ntuple_gain.txt";
  cfg.PedestalFileName = "/data/old_testbeam/btsoftware_1414634301_calib_dark_ntuple.root";

//  cfg.FileName = "/fhgfs/users/ttekampe/SciFi/testbeamData/old/btsoftware_1414634337_datarun_ntuple.root";
//  cfg.PedestalFileName = "/fhgfs/users/ttekampe/SciFi/testbeamData/old/btsoftware_1414634301_calib_dark_ntuple.root";
//  cfg.GainFileName = "/fhgfs/users/ttekampe/SciFi/testbeamData/old/led_gain/btsoftware_1414634319_calib_led_gain_uplink34.txt";
  cfg.debug = false;

//  std::cout << "reading pedestals" << std::endl;
//  std::map<unsigned int, std::map<unsigned int, double>> pedestals = getPedestals(cfg.PedestalFileName, 34, 34, 128);
//
//
//  std::cout << "reading gains" << std::endl;
//  std::map<unsigned int, std::map<unsigned int, double>> gains = readGains(cfg.GainFileName);


  TFile inputFile(cfg.FileName.c_str(), "READ");
  if(!inputFile.IsOpen()){
    std::cerr << "could not open file " << cfg.FileName << std::endl;
    return 0;
  }
  TTree* inputTree = dynamic_cast<TTree*>(inputFile.Get("rawData"));
  if(inputTree == nullptr){
    std::cerr << "inputtree is nullptr" << std::endl;
    return 0;
  }
  std::cout << "start parsing the TTree into vectors" << std::endl;
//  std::vector<std::vector<Channel>*>* dataVector = parseRootTree(inputTree, 34, 34, 128, pedestals, gains);
  std::vector<std::vector<Channel>*>* dataVector = parseCorrectedRootTree(inputTree, 34, 34, 128);
  std::cout << "done. parsed " << dataVector->size() << " events." << std::endl;
  std::vector<Cluster*> myClusters;



  std::cout << "start to search for clusters" << std::endl;
  long double sumOfAdcVals = 0.;
  for(unsigned int i = 0; i < dataVector->size(); i++){
  //int i = 5757;
//    std::cout << "searching for clusters in event numer " << i << " of " << dataVector->size() << std::endl;
    std::vector<Channel> tmpVec = *(dataVector->at(i));
    for(const auto& chan : *(dataVector->at(i))){
      sumOfAdcVals += chan.AdcValue;
    }

    //for(const auto& chan : tmpVec){
    //  std::cout << "channel number " << chan.ChannelNumber << " adc value " << chan.AdcValue << std::endl;

    //}
                                              //neighbor, seed, sum, maxsize
    FindClustersInEvent(myClusters, tmpVec, 1.5, 2.5, 4.5, 100, cfg.debug);
  }
  std::cout << "sum of all adc values: " << sumOfAdcVals << "\n";
  //std::vector<Channel> tmpVec = *(dataVector->at(39));
  //FindClustersInEvent(myClusters, tmpVec, 1.5, 2.5, 4.5);

  std::cout << "done" << std::endl;

  std::cout << "found " << myClusters.size() << " clusters" << std::endl;

  double meanClusterSize = 0.;
  for(const auto& c : myClusters){
      meanClusterSize += c->GetClusterSize();
    //std::cout << "cluster ranges from channel number " << c->GetMinChannel() << " to " << c->GetMaxChannel() << std::endl
    //           << "hit weighted mean: " << c->GetHitWeightedMean() << " charge weighted mean: " << c->GetChargeWeightedMean() << std::endl;
  }
  meanClusterSize /= myClusters.size();
  std::cout << "Mean cluster size: " << meanClusterSize << std::endl;
  return 0;
}
