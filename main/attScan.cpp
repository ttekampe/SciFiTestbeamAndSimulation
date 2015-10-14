//from std
#include <vector>
#include <cmath>
#include <string>
#include <iostream>

//from ROOT
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TAxis.h"
#include "TF1.h"


//from here
#include "Calibration.h"
#include "ClusterCreator.h"


#include "lhcbStyle.h"

int main(){
  lhcb::lhcbStyle();

  std::vector<std::string> files = {
    "305"
    ,"355"
    ,"455"
    ,"555"
    ,"655"
    ,"755"
    ,"855"
    ,"955"
    ,"1055"
    ,"1155"
    ,"1255"
    ,"1355"
    ,"1455"
    ,"1555"
    ,"1655"
    ,"1755"
    ,"1855"
    ,"1955"
    ,"2055"
    ,"2155"
    ,"2255"
  };

  std::vector<double> lightyields;
  std::vector<double> efficiencies;
  std::vector<double> lightyieldsErr;
  std::vector<double> efficienciesErr;

  for(const auto& f : files){
    std::string fname = "/data/testbeam/simulation/attScan/testbeam_simulation_position_" + f + "_at_0deg.root";
    TFile inputFile(fname.c_str(), "READ");
    TTree* inputTree;
    inputTree = dynamic_cast<TTree*>( inputFile.Get("layer_0") );


    if(!inputTree){
      std::cerr << "Tree not found!" << std::endl;
      return 0;
    }

    std::vector<std::vector<Channel>*>* data;
    data = parseCorrectedRootTree(inputTree, 1, 4, 128, true);



    ClusterCreator clCreator;
    int currentNumberOfClusters{0};
    int missedEvents{0};
    for (const auto& event : *data){
      clCreator.FindClustersInEventBoole(*event, 1.5, 2.5, 4.0, 100, false);

      if (currentNumberOfClusters == clCreator.getNumberOfClusters()) ++missedEvents;
      currentNumberOfClusters = clCreator.getNumberOfClusters();
                                                    //neighbor, seed, sum, maxsize, debug in simu 3, 5, 8
    }

    double meanLightYield{0};
    double skippedEvents{0};
    for(const auto& cl : clCreator.getClusters()){
      if(cl->GetClusterSize()>5 || cl->GetChargeWeightedMean()> 97. || cl->GetChargeWeightedMean()< 92.){
        ++skippedEvents;
        continue;
      }
      meanLightYield += cl->GetSumOfAdcValues();
    }
    meanLightYield /= (double)clCreator.getNumberOfClusters();

    double stdDevLightYield{0};
    for(const auto& cl : clCreator.getClusters()){
      if(cl->GetClusterSize()>5 || cl->GetChargeWeightedMean()> 97. || cl->GetChargeWeightedMean()< 92.) continue;
      stdDevLightYield += (cl->GetSumOfAdcValues() - meanLightYield) * (cl->GetSumOfAdcValues() - meanLightYield);
    }
    stdDevLightYield *= 1./( (double)clCreator.getNumberOfClusters() - 1. );
    stdDevLightYield = TMath::Sqrt(stdDevLightYield);

    std::cout << "Mean light yield: " << meanLightYield << " +/- " << stdDevLightYield << "\n";
    lightyields.push_back(meanLightYield);
    lightyieldsErr.push_back(stdDevLightYield/TMath::Sqrt((double)clCreator.getNumberOfClusters()));


    double nEvents = inputTree->GetEntriesFast() - skippedEvents;
    std::cout << "Found " << clCreator.getNumberOfClusters() << " clusters in " << nEvents << " events!\n";
    std::cout << "Missed " << missedEvents << " events\n";
    double event2OneOrMoreClusterEff = 1. - (double)missedEvents/nEvents;
    double event2OneOrMoreClusterEffErr = TMath::Sqrt( (1-event2OneOrMoreClusterEff)*event2OneOrMoreClusterEff/nEvents );
    std::cout << "Rate of MCEvent producing one or more clusters: " << event2OneOrMoreClusterEff << " +/- " << event2OneOrMoreClusterEffErr << "\n";

    efficiencies.push_back(event2OneOrMoreClusterEff);
    efficienciesErr.push_back(event2OneOrMoreClusterEffErr);
  }

  for(const auto& f : files){
    std::string fname = "/data/testbeam/simulation/attScan/testbeam_simulation_position_" + f + "_at_0deg.root";
    TFile inputFile(fname.c_str(), "READ");
    TTree* inputTree;
    inputTree = dynamic_cast<TTree*>( inputFile.Get("layer_0") );


    if(!inputTree){
      std::cerr << "Tree not found!" << std::endl;
      return 0;
    }

    std::vector<std::vector<Channel>*>* data;
    data = parseCorrectedRootTree(inputTree, 1, 4, 128, true, 0.45);



    ClusterCreator clCreator;
    int currentNumberOfClusters{0};
    int missedEvents{0};
    for (const auto& event : *data){
      clCreator.FindClustersInEventBoole(*event, 1.5, 2.5, 4.0, 100, false);

      if (currentNumberOfClusters == clCreator.getNumberOfClusters()) ++missedEvents;
      currentNumberOfClusters = clCreator.getNumberOfClusters();
                                                    //neighbor, seed, sum, maxsize, debug in simu 3, 5, 8
    }

    double meanLightYield{0};
    double skippedEvents{0};
    for(const auto& cl : clCreator.getClusters()){
    //  if(cl->GetClusterSize()>5 || cl->GetChargeWeightedMean()> 97. || cl->GetChargeWeightedMean()< 92.){
    //    ++skippedEvents;
    //    continue;
    //  }
      meanLightYield += cl->GetSumOfAdcValues();
    }
    meanLightYield /= (double)clCreator.getNumberOfClusters();

    double stdDevLightYield{0};
    for(const auto& cl : clCreator.getClusters()){
      if(cl->GetClusterSize()>5 || cl->GetChargeWeightedMean()> 97. || cl->GetChargeWeightedMean()< 92.) continue;
      stdDevLightYield += (cl->GetSumOfAdcValues() - meanLightYield) * (cl->GetSumOfAdcValues() - meanLightYield);
    }
    stdDevLightYield *= 1./( (double)clCreator.getNumberOfClusters() - 1. );
    stdDevLightYield = TMath::Sqrt(stdDevLightYield);

    std::cout << "Mean light yield: " << meanLightYield << " +/- " << stdDevLightYield << "\n";
    lightyields.push_back(meanLightYield);
    lightyieldsErr.push_back(stdDevLightYield/TMath::Sqrt((double)clCreator.getNumberOfClusters()));


    double nEvents = inputTree->GetEntriesFast() - skippedEvents;
    std::cout << "Found " << clCreator.getNumberOfClusters() << " clusters in " << nEvents << " events!\n";
    std::cout << "Missed " << missedEvents << " events\n";
    double event2OneOrMoreClusterEff = 1. - (double)missedEvents/nEvents;
    double event2OneOrMoreClusterEffErr = TMath::Sqrt( (1-event2OneOrMoreClusterEff)*event2OneOrMoreClusterEff/nEvents );
    std::cout << "Rate of MCEvent producing one or more clusters: " << event2OneOrMoreClusterEff << " +/- " << event2OneOrMoreClusterEffErr << "\n";

    efficiencies.push_back(event2OneOrMoreClusterEff);
    efficienciesErr.push_back(event2OneOrMoreClusterEffErr);
  }


  //std::vector<double> zeroes(lightyields.size(), 0.);
  TGraphErrors graph(lightyields.size(),lightyields.data(), efficiencies.data(), lightyieldsErr.data(), efficienciesErr.data());
  //TGraphErrors graph(lightyields.size(),lightyields.data(), efficiencies.data(), zeroes.data(), efficienciesErr.data());

  TF1* func = new TF1("func", "[0]*(1-[1]*exp(-x/[2]))", 0, 22);
  func->SetParameter(0,0.9);
  func->SetParLimits(0,0.,1.);

  func->SetParameter(1,3.5);
  func->SetParLimits(1,0,100);

  func->SetParameter(2,2.7);
  func->SetParLimits(2,1,20);

  TF1* romansFunc = new TF1("romansFunc", "0.99*(1-3.4964*exp(-x/2.7035))", 0, 22);
  romansFunc->SetLineColor(kRed);

  graph.Fit(func, "R");

  TCanvas can;
  graph.Draw("ALP");
  graph.GetXaxis()->SetTitle("Light yield in p.e.");
  graph.GetYaxis()->SetTitle("#varepsilon");
  graph.Draw("AP");
  graph.Print("v");

  romansFunc->Draw("same");

  can.SaveAs("lightyieldVsEff_artReduced.pdf");

  delete  func;

}
