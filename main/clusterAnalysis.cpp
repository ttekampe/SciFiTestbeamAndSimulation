//from std
#include <vector>
#include <cmath>
#include <iostream>

//from ROOT
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2I.h"

//from here
//#include "Cluster.h"
//#include "ClusterAlgorithms.h"
#include "Calibration.h"
#include "ClusterCreator.h"
#include "ClusterMonitor.h"

//from BOOST
#include "boost/program_options.hpp"

namespace po = boost::program_options;


struct config{
  std::string file2analyse;
  bool debug;
  bool simulation;
  std::string clusterAlg;
};


void addClusterShapeToHist(const Cluster& cl, std::vector<std::vector<double> >& adcValuesRelativeToSeed){
  int clusterSeed = cl.GetSeedChannelNumber();
  for(const auto& chan : cl.GetRelatedChannels()){
    unsigned int position = std::abs(int(chan.ChannelNumber) - clusterSeed);
    if (cl.GetClusterSize() <= adcValuesRelativeToSeed.size()){
      adcValuesRelativeToSeed.at(position).push_back(chan.AdcValue);
    }
  }
}

std::pair<double , double> getMean(std::vector<double> data){
  if(data.empty()){
    return std::make_pair(0., 0.);
  }
  double mean{0.};
  for(const auto& entry : data){
    mean += entry;
  }
  mean /= data.size();

  double stdDev{0.};
  for(const auto& entry : data){
    stdDev += (mean - entry) * (mean - entry);
  }
  stdDev = std::sqrt(stdDev / (data.size() - 1.) );

  return std::make_pair(mean, stdDev);
}


int parseOptions(config &c, int argc, char *argv[]){

  // declare options
  po::options_description desc("Allowed options") ;
  desc.add_options()
    ("help", "show this help")
    ("file,f", po::value<std::string>(&c.file2analyse), "corrected test beam data file")
    ("simulation,s", po::bool_switch(&c.simulation), "Simulated input?")
    ("clusteralg,c", po::value<std::string>(&c.clusterAlg)->default_value("b"), "clustering algorithm: b for Boole or m for Maxs")
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

  TFile inputFile(c.file2analyse.c_str(), "READ");
  TTree* inputTree;
  if(c.simulation){
    inputTree = dynamic_cast<TTree*>( inputFile.Get("layer_0") );
  }
  else{
    inputTree = dynamic_cast<TTree*>( inputFile.Get("rawData") );
  }

  if(!inputTree){
    std::cerr << "Tree not found!" << std::endl;
    return 0;
  }

  std::vector<std::vector<Channel>*>* data;

  if(c.simulation){
    data = parseCorrectedRootTree(inputTree, 1, 4, 128, true);
  }
  else{
    data = parseCorrectedRootTree(inputTree, 3, 4, 128, false);
  }


  ClusterCreator clCreator;
  for (const auto& event : *data){

    if(c.clusterAlg == "b") clCreator.FindClustersInEventBoole(*event, 1.5, 2.5, 4.0, 100, false);
    if(c.clusterAlg == "m") clCreator.FindClustersInEventMax(*event, 1.5, 2.5, 4.0);
                                                  //neighbor, seed, sum, maxsize, debug in simu 3, 5, 8
  }

  std::cout << "Found " << clCreator.getNumberOfClusters() << " clusters in " << inputTree->GetEntriesFast() << " events!\n";
  std::string tag = "";
  if (c.clusterAlg == "m") tag = "_max";

  ClusterMonitor clMonitor;
  clMonitor.WriteToNtuple(clCreator, ("/home/tobi/SciFi/results/clusters/" + removePath(c.file2analyse).ReplaceAll(".root", "_clusterAnalyis" + tag +".root")).Data() );


//  std::vector<Cluster*> clusterVec;
//
//  for (const auto& event : *data){
//                                                  //neighbor, seed, sum, maxsize simu 3, 5, 8
//    FindClustersInEventBoole(clusterVec, *event, 1.5, 2.5, 4.5, 100, false);
////    FindClustersInEventFTStyle(clusterVec, *event, 0.5, 1.5, 1.5, 100, false);
//
//
//  }
//
//  std::cout << "Found " << clusterVec.size() << " clusters in " << inputTree->GetEntriesFast() << " events!\n";
//
//  TFile results("/home/tobi/SciFi/results/clusters/" + removePath(c.file2analyse).ReplaceAll(".root", "_clusterAnalyis.root"), "RECREATE");
//  TTree* resultTree = new TTree("clusterAnalysis", "");
//
//  double d_chargeWeightedMean = 0.0;
//  resultTree->Branch("chargeWeightedMean", &d_chargeWeightedMean, "chargeWeightedMean/D");
//
//  double d_hitWeightedMean = 0.0;
//  resultTree->Branch("hitWeightedMean", &d_hitWeightedMean, "hitWeightedMean/D");
//
//  double d_sumCharge = 0.0;
//  resultTree->Branch("sumCharge", &d_sumCharge, "sumCharge/D");
//
//  double d_maxCharge = 0.0;
//  resultTree->Branch("maxCharge", &d_maxCharge, "maxCharge/D");
//
//  int i_clusterSize = 0;
//  resultTree->Branch("clusterSize", &i_clusterSize, "clusterSize/I");
//
////  double d_clusterShape = 0.0;
////  resultTree->Branch("clusterShape", &d_clusterShape, "clusterShape/D");
//  TH1D clusterShape("clusterShape", "", 10, 0, 10);
//
//  std::vector<std::vector<double> > adcValuesRelativeToSeed(10);
//
//
//  for(const auto& clust : clusterVec){
//    addClusterShapeToHist(*clust, adcValuesRelativeToSeed);
//    d_chargeWeightedMean = clust->GetChargeWeightedMean();
//    d_hitWeightedMean = clust->GetHitWeightedMean();
//    d_sumCharge = clust->GetSumOfAdcValues();
//    d_maxCharge = clust->GetMaximumAdcValue();
//    i_clusterSize = clust->GetClusterSize();
//    resultTree->Fill();
//  }
//  std::cout << "Filling cluster shape histogram..\n";
//  std::pair<double, double> current(0., 0.);
//  for(unsigned int i = 0; i<adcValuesRelativeToSeed.size(); ++i){
//    current = getMean(adcValuesRelativeToSeed.at(i));
//    clusterShape.SetBinContent(i+1, current.first);
//    clusterShape.SetBinError(i+1, current.second);
//  }
//
//  clusterShape.Write();
//  resultTree->Write();
//  results.Close();


  return 0;
}
