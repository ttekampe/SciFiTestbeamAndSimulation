//from std
#include <vector>
#include <cmath>
#include <iostream>
#include <map>
#include <algorithm> 
#include <fstream>

//from ROOT
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2I.h"
#include "TMath.h"
#include "TGraph.h"
#include "TF1.h"
#include "TCanvas.h"

//from RooFit
#include "RooFit.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooFitResult.h"

//from here
//#include "Cluster.h"
//#include "ClusterAlgorithms.h"
#include "Calibration.h"
#include "ClusterCreator.h"
#include "ClusterMonitor.h"
#include "EDouble.h"

//from BOOST
#include "boost/program_options.hpp"

namespace po = boost::program_options;

struct config{
  std::vector<std::string> files2analyse;
  bool debug;
  bool simulation;
  std::string clusterAlg;
  std::string tag;
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
  ("file,f", po::value<std::vector<std::string>>(&c.files2analyse)->multitoken(), "corrected test beam data file")
  ("simulation,s", po::bool_switch(&c.simulation), "Simulated input?")
  // ("clusteralg,c", po::value<std::string>(&c.clusterAlg)->default_value("b"), "clustering algorithm: b for Boole or m for Maxs")
  ("tag,t", po::value<std::string>(&c.tag)->default_value(""), "tag that is added to the output file name")
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

std::string getPositionFromFileName(std::string fileName){
  std::string position;
  if(fileName.find("btsoftware") != std::string::npos){ // data
    std::map<std::string, std::string> runNumbersAndPositions;
    //Pos Angle   Runnumber     Frac of events with 1 or more clusters
    runNumbersAndPositions["1431786652"] =  "a";    
    runNumbersAndPositions["1432091294"] =  "a";    
    runNumbersAndPositions["1432264510"] =  "a";    
    runNumbersAndPositions["1432350729"] =  "a";    
    runNumbersAndPositions["1432341050"] =  "a";    
    runNumbersAndPositions["1432169457"] =  "c";   
    runNumbersAndPositions["1432089102"] =  "c";   
    runNumbersAndPositions["1432187205"] =  "c";   
    runNumbersAndPositions["1432358809"] =  "c";   

    //AttScan
    runNumbersAndPositions["1431883436"] = "35.5";
    runNumbersAndPositions["1431880915"] = "45.5";
    runNumbersAndPositions["1431853383"] = "55.5";
    runNumbersAndPositions["1431873726"] = "65.5";
    runNumbersAndPositions["1431876961"] = "75.5";
    runNumbersAndPositions["1431847243"] = "85.5";
    runNumbersAndPositions["1431871033"] = "95.5";
    runNumbersAndPositions["1431857717"] = "105.5";
    runNumbersAndPositions["1431868949"] = "115.5";
    runNumbersAndPositions["1431866804"] = "135.5";
    runNumbersAndPositions["1431886395"] = "145.5";
    runNumbersAndPositions["1431850788"] = "155.5";
    runNumbersAndPositions["1431893934"] = "165.5";
    runNumbersAndPositions["1431864557"] = "175.5";
    runNumbersAndPositions["1431845031"] = "185.5";
    runNumbersAndPositions["1431896347"] = "195.5";
    runNumbersAndPositions["1431855605"] = "205.5";
    runNumbersAndPositions["1431899004"] = "215.5";

    std::string currentRunNumber = std::to_string(runNumberFromFilename(fileName));
    position = runNumbersAndPositions.find(currentRunNumber) != runNumbersAndPositions.end() ? runNumbersAndPositions[currentRunNumber] : "";

  }
  else{ //simulation
    std::string lookFor = "testbeam_simulation_position_";
    auto posStartAt = fileName.find(lookFor) + lookFor.size();
    auto posEndAt = fileName.find("_", posStartAt);
    position = fileName.substr(posStartAt, posEndAt-posStartAt);
    std::cout << "Current position is " << position << "\n";
  }
  return position;
}



std::pair<EDouble, EDouble> analyse(std::string file2analyse, const config& c){
  
  TFile inputFile(file2analyse.c_str(), "READ");
  if(!inputFile.IsOpen()){
    std::cout << "Could not open " << file2analyse << "\n";
    return std::make_pair(EDouble(0, 0), EDouble(0, 0));
  }

  std::cout << "Analysis " << file2analyse << "\n";

  TTree* inputTree;
  if(c.simulation){
    inputTree = dynamic_cast<TTree*>( inputFile.Get("layer_0") );
  }
  else{
    inputTree = dynamic_cast<TTree*>( inputFile.Get("rawData") );
  }

  if(!inputTree){
    std::cerr << "Tree not found!" << std::endl;
    return std::make_pair(EDouble(0, 0), EDouble(0, 0));
  }

  /*
  In the 2015 testbeam 4 fibre mats were tested:
  HD2 Uplink 5 and 6
  HD1 Uplink 7 and 8
  CERN4 Uplink 1 and 2
  SLAYER3(DUT) Uplink 3 and 4
  */

  std::map<std::string, std::vector<std::vector<Channel>*>* > data;

  std::vector<double> xPositions;
  std::string offsetFileName = ( "/home/ttekampe/SciFi/results/moduleOffset/" + removePath(file2analyse).ReplaceAll(".root", ".txt") ).Data();
  bool produceOffsetFile{false};
  std::vector<double> xOffsets;
  std::vector<double> track_distances;

  if(!c.simulation){
    std::cout << "Checking for offsetFile...\n";
     std::ifstream offsetFile(offsetFileName);
    if(!offsetFile){
      std::cout << "Offsetfile does not exists yet -> producing...\n";
      produceOffsetFile = true;
    }
    else{
      std::string line;
      std::getline(offsetFile, line);
      std::istringstream ss(line);
      double offset;
      if(ss >> offset) xOffsets.push_back(offset);
      std::cout << "Found xoffset: " << xOffsets[0] << "\n";
    }
  }

  std::vector<double> zPositions = {0., 247.0*1000, 469.0*1000}; // HD2, Slayer, CERN

  if(c.simulation){
    data["simulation"] = parseCorrectedRootTree(inputTree, 1, 4, 128);
  }
  else{
    data["cern"] = parseCorrectedRootTree(inputTree, 1, 2, 128);
    data["slayer"] = parseCorrectedRootTree(inputTree, 3, 4, 128);
    data["HD2"] = parseCorrectedRootTree(inputTree, 5, 6, 128);
  }

  std::map<std::string, ClusterCreator> clCreators;
  if(!c.simulation){
    clCreators["cern"] = ClusterCreator();
    clCreators["slayer"] = ClusterCreator();
    clCreators["HD2"] = ClusterCreator();
  }

  clCreators["simulation"] = ClusterCreator(); // in the case of data this one stores the mathed clusters

  std::map<std::string, std::vector<Cluster*>> clustersInModule;
  clustersInModule["cern"];
  clustersInModule["slayer"];
  clustersInModule["HD2"];

  ClusterCreator cl2Analyse;
  unsigned int currentNumberOfClusters{0};
  double missedEvents{0};
  double foundEvents{0};
  if(c.simulation){
    for (const auto& event : *data["simulation"]){

      // if(c.clusterAlg == "b") 
      clCreators["simulation"].FindClustersInEventBoole(*event, 1.5, 2.5, 4.0, 100, false);
      // if(c.clusterAlg == "m") clCreators["simulation"].FindClustersInEventMax(*event, 1.5, 2.5, 4.0);
      if (currentNumberOfClusters == clCreators["simulation"].getNumberOfClusters()) ++missedEvents;
      currentNumberOfClusters = clCreators["simulation"].getNumberOfClusters();
                                                      //neighbor, seed, sum, maxsize, debug in simu 3, 5, 8
    }
  }
  else{
    for(unsigned int i = 0; i<inputTree->GetEntriesFast(); ++i){
      clustersInModule["cern"] = clCreators["cern"].FindClustersInEventBoole(*(data["cern"]->at(i)), 1.5, 2.5, 4.0, 100, false);
      clustersInModule["slayer"] = clCreators["slayer"].FindClustersInEventBoole(*(data["slayer"]->at(i)), 1.5, 2.5, 4.0, 100, false);
      clustersInModule["HD2"] = clCreators["HD2"].FindClustersInEventBoole(*(data["HD2"]->at(i)), 1.5, 2.5, 4.0, 100, false);
      
      if(  clustersInModule["cern"].size() == 1
        && clustersInModule["slayer"].size() == 1
        && clustersInModule["HD2"].size() == 1
        ){

        clustersInModule["slayer"] = clCreators["simulation"].FindClustersInEventBoole(*(data["slayer"]->at(i)), 1.5, 2.5, 4.0, 100, false);


        xPositions = {
          clustersInModule["HD2"][0]->GetChargeWeightedMean() * 250.,
          clustersInModule["slayer"][0]->GetChargeWeightedMean() * 250.,
          clustersInModule["cern"][0]->GetChargeWeightedMean() * 250.

        };

        double dx = xPositions[2] - xPositions[0];
        double dz = zPositions[2] - zPositions[0];

        double slope =  dx / dz;
        double constant = xPositions[0];

        //std::cout << "Track is " << constant << " + " << slope << " * x\n";

        double trackAtSlayer = constant + slope * zPositions[1];
        //std::cout << "Track at slayer should be " << trackAtSlayer << " and is " << xPositions[1] << "\n";
        if(produceOffsetFile) xOffsets.push_back(trackAtSlayer - xPositions[1]);
        else{
          track_distances.push_back(trackAtSlayer - xPositions[1] - xOffsets[0]);
        }


      }

      //std::cout << "\n";
      if(!clustersInModule["cern"].empty() && !clustersInModule["HD2"].empty()) ++missedEvents;

      if(!clustersInModule["cern"].empty() && !clustersInModule["slayer"].empty() && !clustersInModule["HD2"].empty()){ // if cluster in all modules store the slayer one for analysis
        //clCreators["simulation"].FindClustersInEventBoole((*data["slayer"]->at(i)), 1.5, 2.5, 4.0, 100, false);
        ++foundEvents;

      }

    }

  }

    if(produceOffsetFile){
      //fit gaussian to offsets
      TCanvas can_offset;

      std::cout << "Collected " << xOffsets.size() << " xoffsets\n";

      double sum = std::accumulate(xOffsets.begin(), xOffsets.end(), 0.0);
      double mean = sum / xOffsets.size();

      double sq_sum = std::inner_product(xOffsets.begin(), xOffsets.end(), xOffsets.begin(), 0.0);
      double stdev = std::sqrt(sq_sum / xOffsets.size() - mean * mean);

      std::cout << "mean: " << mean << " stdev: " << stdev << "\n";

      RooRealVar rv_xoffset("rv_xoffset", "", mean - 300 , mean + 300);

      RooDataSet dataset_xoffset("dataset_xoffset", "", RooArgSet(rv_xoffset));
      for(const auto xoffset : xOffsets){
        if(xoffset > mean + 300 || xoffset < mean - 300) continue;

        //std::cout << "Adding " << xoffset << " to offset dataset\n";
        rv_xoffset.setVal(xoffset);
        dataset_xoffset.add(RooArgSet(rv_xoffset));
      }
      dataset_xoffset.Print();
      RooRealVar rv_mean("rv_mean", "", mean, mean - 100 , mean + 100);
      RooRealVar rv_sigma("rv_sigma", "", 50, 0,  100);  
      RooGaussian gauss("gauss", "", rv_xoffset, rv_mean, rv_sigma);

      RooFitResult* fr = gauss.fitTo(dataset_xoffset, RooFit::Save(true), RooFit::Minimizer("Minuit2", "minimize"));

      RooPlot* plot = rv_xoffset.frame();
      dataset_xoffset.plotOn(plot);
      gauss.plotOn(plot);
      plot->Draw();

      can_offset.SaveAs("/home/ttekampe/SciFi/results/moduleOffset/" + removePath(file2analyse).ReplaceAll(".root", ".pdf") );

      fr->Print("v");

      std::ofstream offsetFile( ("/home/ttekampe/SciFi/results/moduleOffset/" + removePath(file2analyse).ReplaceAll(".root", ".txt") ).Data() );
      offsetFile << rv_mean.getVal() << "\n";
      offsetFile.close();

      // restart and run with offsets
      std::map<std::string, std::vector<std::vector<Channel>*>* > data;
      for (auto& module : data){
        for(unsigned int entryIndex = 0; entryIndex < module.second->size(); ++entryIndex){
          delete module.second->at(entryIndex);
        }
      }

      delete inputTree;
      inputFile.Close();
      return analyse(file2analyse, c);
    }

    double meanLightYield{0};
    for(const auto& cl : clCreators["simulation"].getClusters()){
      meanLightYield += cl->GetSumOfAdcValues();
    }
    meanLightYield /= (double)clCreators["simulation"].getNumberOfClusters();

    double stdDevLightYield{0};
    for(const auto& cl : clCreators["simulation"].getClusters()){
      stdDevLightYield += (cl->GetSumOfAdcValues() - meanLightYield) * (cl->GetSumOfAdcValues() - meanLightYield);
    }

    stdDevLightYield *= 1./((double)clCreators["simulation"].getNumberOfClusters() - 1.);
    stdDevLightYield = TMath::Sqrt(stdDevLightYield) / TMath::Sqrt((double)clCreators["simulation"].getNumberOfClusters());

    std::cout << "Mean light yield: " << meanLightYield << " +/- " << stdDevLightYield << "\n";

    EDouble light;
    EDouble eff;

    if(c.simulation){
      std::cout << "Found " << clCreators["simulation"].getNumberOfClusters() << " clusters in " << inputTree->GetEntriesFast() << " events!\n";
      std::cout << "Missed " << missedEvents << " events\n";
      double event2OneOrMoreClusterEff = 1. - (double)missedEvents/data["simulation"]->size();
      double event2OneOrMoreClusterEffErr = TMath::Sqrt( (1-event2OneOrMoreClusterEff)*event2OneOrMoreClusterEff/(double)data["simulation"]->size() );
      std::cout << "Rate of MCEvent producing one or more clusters: " << event2OneOrMoreClusterEff << " +/- "
      <<  event2OneOrMoreClusterEffErr << "\n";

      light = EDouble(meanLightYield, stdDevLightYield);
      eff = EDouble(event2OneOrMoreClusterEff, event2OneOrMoreClusterEffErr);
      //hitEffFile << getPositionFromFileName(file2analyse) << "," << meanLightYield << "," << stdDevLightYield << ","
      //     << event2OneOrMoreClusterEff << "," <<  event2OneOrMoreClusterEffErr << "\n";
    }
    else{
      std::cout << "Found " << missedEvents << " events that caused one or more cluster in all three but the slayer modules.\n";
      std::cout << "Found " << foundEvents << " events that caused one or more clusters in all modules\n";
      double efficiency = foundEvents / missedEvents;
      double effErr = TMath::Sqrt(efficiency*(1.-efficiency) / missedEvents);
      std::cout << "This is the fraction of " << foundEvents / missedEvents << "\n";

      // hitEffFile << "data" << "," << meanLightYield << "," << stdDevLightYield << ","
      //      << efficiency << "," <<  effErr << "\n";

      light = EDouble(meanLightYield, stdDevLightYield);
      eff = EDouble(efficiency, effErr);


    }


    // if (c.clusterAlg == "m") c.tag += "_max";

    std::map<std::string, std::vector<double>> features;
    features["distance_from_track"] = track_distances;

    ClusterMonitor clMonitor;
    clMonitor.WriteToNtuple(clCreators["simulation"], 
      ("/home/ttekampe/SciFi/results/clusters/" + removePath(file2analyse).ReplaceAll(".root", "_clusterAnalyis" + c.tag +".root")).Data(),
      features );
      for (auto& module : data){
       for(unsigned int entryIndex = 0; entryIndex < module.second->size(); ++entryIndex){
          delete module.second->at(entryIndex);
        }
      }
    delete inputTree;
    inputFile.Close();


    return std::make_pair(light, eff);
  }



int main(int argc, char *argv[]){

  config c;
  if (parseOptions(c, argc, argv)!= 0){
    std::cerr << "Can not parse options!" << std::endl;
    return 0;
  }

  std::ofstream hitEffFile("/home/ttekampe/SciFi/results/hitEfficiency.txt");
  hitEffFile << "position\tlightyield\tlightyieldErr\tefficiency\tefficiencyErr\n";
  
  std::pair<EDouble, EDouble> lightAndEff;

  for(const auto& file2analyse : c.files2analyse){

    lightAndEff = analyse(file2analyse, c);
    hitEffFile << getPositionFromFileName(file2analyse) << "," << lightAndEff.first.GetVal() << "," << lightAndEff.first.GetError() << ","
               << lightAndEff.second.GetVal() << "," <<  lightAndEff.second.GetError() << "\n";

  }

  return 0;
}
