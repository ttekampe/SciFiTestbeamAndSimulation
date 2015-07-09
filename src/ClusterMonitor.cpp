//from std
#include <vector>
#include <utility>

//from ROOT
#include "TFile.h"
#include "TTree.h"

//from here
#include "Cluster.h"
#include "ClusterMonitor.h"




void ClusterMonitor::WriteToNtuple(const ClusterCreator& clCreator, std::string fileName){
  TFile results(fileName.c_str(), "RECREATE");
  TTree* resultTree = new TTree("clusterAnalysis", "");

  double d_chargeWeightedMean{0.0};
  resultTree->Branch("chargeWeightedMean", &d_chargeWeightedMean, "chargeWeightedMean/D");

  double d_hitWeightedMean{0.0};
  resultTree->Branch("hitWeightedMean", &d_hitWeightedMean, "hitWeightedMean/D");

  double d_sumCharge{0.0};
  resultTree->Branch("sumCharge", &d_sumCharge, "sumCharge/D");

  double d_maxCharge{0.0};
  resultTree->Branch("maxCharge", &d_maxCharge, "maxCharge/D");

  int i_clusterSize{0};
  resultTree->Branch("clusterSize", &i_clusterSize, "clusterSize/I");


  for(const auto& clust : clCreator.getClusters()){
    d_chargeWeightedMean = clust->GetChargeWeightedMean();
    d_hitWeightedMean = clust->GetHitWeightedMean();
    d_sumCharge = clust->GetSumOfAdcValues();
    d_maxCharge = clust->GetMaximumAdcValue();
    i_clusterSize = clust->GetClusterSize();
    resultTree->Fill();
  }

  resultTree->Write();
  results.Close();


}
