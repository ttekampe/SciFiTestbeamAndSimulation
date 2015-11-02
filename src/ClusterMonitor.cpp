//from std
#include <iostream>
#include <vector>
#include <utility>

//from ROOT
#include "TFile.h"
#include "TTree.h"

//from here
#include "Cluster.h"
#include "ClusterMonitor.h"




void ClusterMonitor::WriteToNtuple(const ClusterCreator& clCreator, const std::string fileName, const feature_map features){
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

  std::map<std::string, double*> ptrMap;
  for(const auto& feature : features){
    if(feature.second.size() != clCreator.getNumberOfClusters()){
      std::cerr << "Cannot add feature " << feature.first << " because the vector size does not match the number of found clusters\n";
      std::cerr << clCreator.getNumberOfClusters() << " clusters and " << feature.second.size() << " values for the feature!\n"; 
      continue;
    }
    std::cout << "Adding feature " << feature.first << " to result tree\n";
    ptrMap[feature.first] = new double();
    resultTree->Branch(feature.first.c_str(), ptrMap[feature.first], (feature.first+"/D").c_str() );
  }

  int index = 0;
  for(const auto& clust : clCreator.getClusters()){
    d_chargeWeightedMean = clust->GetChargeWeightedMean();
    d_hitWeightedMean = clust->GetHitWeightedMean();
    d_sumCharge = clust->GetSumOfAdcValues();
    d_maxCharge = clust->GetMaximumAdcValue();
    i_clusterSize = clust->GetClusterSize();

    for(auto& ptr : ptrMap){
      *(ptr.second) = features.at(ptr.first).at(index);
    }

    resultTree->Fill();
    ++index;
  }

  resultTree->Write();
  results.Close();


}
