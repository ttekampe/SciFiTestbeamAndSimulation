//from here
#include "Cluster.h"
#include "ClusterAlgorithms.h"
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


std::vector<double> GetGains(std::string fileName, int nAdcs){
  std::ifstream inputFile(fileName);
  std::string line;
  int chanNumber;
  double gain;

  std::vector<double> gains(nAdcs);
  int i = 0;
  while (std::getline(inputFile, line) && i < nAdcs + 4) {
    if (i >= 4) { //skip first 4 lines with textfile discription
      std::cout << line << std::endl;
      std::istringstream ss(line);
      ss >> chanNumber >> gain;
      gains[i-4] = gain;
    }
    ++i;
  }
  return gains;
}


std::vector<double> GetPedestals(std::string fileName, int uplink, int nAdcs){
  TFile pedestalFile(fileName.c_str(), "READ");
  if(!pedestalFile.IsOpen()){
    std::cerr << "unable to open pedestal file " << fileName << std::endl;
  }
  TTree* pedestalTree = dynamic_cast<TTree*>(pedestalFile.Get("rawData"));
  if(pedestalTree == nullptr){
    std::cerr << "pedestal tree is nullptr" << std::endl;
  }
  std::string branchName = "Uplink_" + std::to_string(uplink) + "_adc_";

  std::vector<double> pedestals(nAdcs, 0.);
  float* adcVals = new float[nAdcs];
  for(int i = 0; i<nAdcs; ++i){
    pedestalTree->SetBranchAddress((branchName + std::to_string(i+1)).c_str(), &adcVals[i]);
  }
  for(int i = 0; i < pedestalTree->GetEntriesFast(); ++i){
    pedestalTree->GetEntry(i);
    for(int j = 0; j < nAdcs; ++j){
      pedestals[j] += adcVals[j];
    }
  }
  for(int j = 0; j < nAdcs; ++j){
    pedestals[j] /= pedestalTree->GetEntriesFast();
  }
  pedestalFile.Close();
  delete[] adcVals;
  return pedestals;
}

std::vector<std::vector<Channel>*>* parseRootTree(TTree* dataTree, int uplink, int nAdcs, std::vector<double>& pedestals, std::vector<double>& gains){
  std::cout << "Reading " << nAdcs << " channels from uplink " << uplink << std::endl;
  std::string branchName = "Uplink_" + std::to_string(uplink) + "_adc_";
  std::vector<std::vector<Channel>*>* dataVector = new std::vector<std::vector<Channel>*>(dataTree->GetEntriesFast());

  float* adcVals = new float[nAdcs];
  for(int i = 0; i<nAdcs; ++i){
    dataTree->SetBranchAddress((branchName + std::to_string(i+1)).c_str(), &adcVals[i]);
  }

  for(int i = 0; i < dataTree->GetEntriesFast(); ++i){
    dataTree->GetEntry(i);
    std::vector<Channel>* event = new std::vector<Channel>(nAdcs);
    for(int j = 0; j < nAdcs; ++j){
      Channel c;
      c.ChannelNumber = j;
      //std::cout << "adcValue = (" << adcVals[j] << " - " << pedestals[j] <<  " )/ " << gains[j] << std::endl;
      c.AdcValue = (adcVals[j] - pedestals[j]) / gains[j];
      event->at(j) = c;
    }
    dataVector->at(i) = event;
  }
  delete adcVals;
  return dataVector;
}

struct CalibrationFiles{
  std::string pedestal;
  std::string gain;
};

int main(){
  Config cfg;
  cfg.FileName = "/fhgfs/users/ttekampe/SciFi/testbeamData/old/btsoftware_1414634337_datarun_ntuple.root";
  cfg.PedestalFileName = "/fhgfs/users/ttekampe/SciFi/testbeamData/old/btsoftware_1414634301_calib_dark_ntuple.root";
  cfg.GainFileName = "/fhgfs/users/ttekampe/SciFi/testbeamData/old/led_gain/btsoftware_1414634319_calib_led_gain_uplink34.txt";
  cfg.debug = false;

  std::cout << "reading pedestals" << std::endl;
  std::vector<double> pedestals = GetPedestals(cfg.PedestalFileName, 34, 128);


  std::cout << "reading gains" << std::endl;
  std::vector<double> gains = GetGains(cfg.GainFileName, 128);

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
  std::vector<std::vector<Channel>*>* dataVector = parseRootTree(inputTree, 34, 128, pedestals, gains);
  std::cout << "done" << std::endl;
  std::vector<Cluster*> myClusters;



  std::cout << "start to search for clusters" << std::endl;

  for(unsigned int i = 0; i < dataVector->size(); i++){
  //int i = 5757;
//    std::cout << "searching for clusters in event numer " << i << " of " << dataVector->size() << std::endl;
    std::vector<Channel> tmpVec = *(dataVector->at(i));

    //for(const auto& chan : tmpVec){
    //  std::cout << "channel number " << chan.ChannelNumber << " adc value " << chan.AdcValue << std::endl;

    //}

    FindClustersInEventFTStyle(myClusters, tmpVec, 1.5, 2.5, 4.5, 100, cfg.debug);
  }

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
