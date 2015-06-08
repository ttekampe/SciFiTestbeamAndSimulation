//from std
#include <vector>

//from ROOT
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2I.h"

//from here
#include "Cluster.h"
#include "ClusterAlgorithms.h"
#include "Calibration.h"

//from BOOST
#include "boost/program_options.hpp"

namespace po = boost::program_options;


struct config{
  std::string file2analyse;
  bool debug;
};



int parseOptions(config &c, int argc, char *argv[]){

  // declare options
  po::options_description desc("Allowed options") ;
  desc.add_options()
    ("help", "show this help")
    ("file,f", po::value<std::string>(&c.file2analyse), "corrected test beam data file")
    ("debug,2", po::bool_switch(&c.debug), "debug output")
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
  TTree* inputTree = dynamic_cast<TTree*>( inputFile.Get("rawData") );
  if(!inputTree){
    std::cerr << "Tree not found!" << std::endl;
    return 0;
  }

  std::vector<std::vector<Channel>*>* data = parseCorrectedRootTree(inputTree, 3, 4, 128);

  std::vector<Cluster*> clusterVec;

  for (const auto& event : *data){
                                                  //neighbor, seed, sum, maxsize
    FindClustersInEventFTStyle(clusterVec, *event, 1.5, 2.5, 4.5, 100, false);
//    FindClustersInEvent(clusterVec, *event, 1.5, 2.5, 4.5, 4, false);
  }

  std::cout << "Found " << clusterVec.size() << " clusters in " << inputTree->GetEntriesFast() << " events!\n";

  TFile results("/home/tobi/SciFi/results/clusters/" + removePath(c.file2analyse).ReplaceAll(".root", "_clusterAnalyis.root"), "RECREATE");
  TTree* resultTree = new TTree("clusterAnalysis", "");

  double d_chargeWeightedMean = 0.0;
  resultTree->Branch("chargeWeightedMean", &d_chargeWeightedMean, "chargeWeightedMean/D");

  double d_hitWeightedMean = 0.0;
  resultTree->Branch("hitWeightedMean", &d_hitWeightedMean, "hitWeightedMean/D");

  double d_sumCharge = 0.0;
  resultTree->Branch("sumCharge", &d_sumCharge, "sumCharge/D");

  double d_maxCharge = 0.0;
  resultTree->Branch("maxCharge", &d_maxCharge, "maxCharge/D");

  int i_clusterSize = 0.0;
  resultTree->Branch("clusterSize", &i_clusterSize, "clusterSize/I");


//  TH1D chargeWeightedMeanHist("chargeWeightedMeanHist", "", 6, 0, 6)
//  TH1D hitWeightedMeanHist
//  TH1D sumChargeHist
//  TH1D maxChargeHist
//  TH1I clusterSizeHist

  for(const auto& clust : clusterVec){
    d_chargeWeightedMean = clust->GetChargeWeightedMean();
    d_hitWeightedMean = clust->GetHitWeightedMean();
    d_sumCharge = clust->GetSumOfAdcValues();
    d_maxCharge = clust->GetMaximumAdcValue();
    i_clusterSize = clust->GetClusterSize();
    resultTree->Fill();
  }

  resultTree->Write();
  results.Close();


  return 0;
}
