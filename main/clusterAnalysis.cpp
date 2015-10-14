//from std
#include <vector>
#include <cmath>
#include <iostream>
#include <map>

//from ROOT
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2I.h"
#include "TMath.h"

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
    ("clusteralg,c", po::value<std::string>(&c.clusterAlg)->default_value("b"), "clustering algorithm: b for Boole or m for Maxs")
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

int main(int argc, char *argv[]){

  config c;
  if (parseOptions(c, argc, argv)!= 0){
    std::cerr << "Can not parse options!" << std::endl;
    return 0;
  }


  for(const auto& file2analyse : c.files2analyse){

    TFile inputFile(file2analyse.c_str(), "READ");
    if(!inputFile.IsOpen()){
      std::cout << "Could not open " << file2analyse << "\n";
      continue;
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
      return 0;
    }

    /*
    In the 2015 testbeam 4 fibre mats were tested:
    HD2 Uplink 5 and 6
    HD1 Uplink 7 and 8
    CERN4 Uplink 1 and 2
    SLAYER3(DUT) Uplink 3 and 4
    */

    std::map<std::string, std::vector<std::vector<Channel>*>* > data;

    if(c.simulation){
      data["simulation"] = parseCorrectedRootTree(inputTree, 1, 4, 128, true);
    }
    else{
      data["cern"] = parseCorrectedRootTree(inputTree, 1, 2, 128, false);
      data["slayer"] = parseCorrectedRootTree(inputTree, 3, 4, 128, false);
      data["HD2"] = parseCorrectedRootTree(inputTree, 5, 6, 128, false);
      data["HD1"] = parseCorrectedRootTree(inputTree, 7, 8, 128, false);
    }

    std::map<std::string, ClusterCreator> clCreators;
    if(!c.simulation){
      clCreators["cern"] = ClusterCreator();
      clCreators["slayer"] = ClusterCreator();
      clCreators["HD2"] = ClusterCreator();
      clCreators["HD1"] = ClusterCreator();
    }

    clCreators["simulation"] = ClusterCreator(); // in the case of data this one stores the mathed clusters

    std::map<std::string, bool> clusterInModule;
    clusterInModule["cern"] = false;
    clusterInModule["slayer"] = false;
    clusterInModule["HD2"] = false;
    clusterInModule["HD1"] = false;

    ClusterCreator cl2Analyse;
    int currentNumberOfClusters{0};
    double missedEvents{0};
    double foundEvents{0};
    if(c.simulation){
        for (const auto& event : *data["simulation"]){

          if(c.clusterAlg == "b") clCreators["simulation"].FindClustersInEventBoole(*event, 1.5, 2.5, 4.0, 100, false);
          if(c.clusterAlg == "m") clCreators["simulation"].FindClustersInEventMax(*event, 1.5, 2.5, 4.0);
          if (currentNumberOfClusters == clCreators["simulation"].getNumberOfClusters()) ++missedEvents;
          currentNumberOfClusters = clCreators["simulation"].getNumberOfClusters();
                                                        //neighbor, seed, sum, maxsize, debug in simu 3, 5, 8
        }
    }
    else{
        for(unsigned int i = 0; i<inputTree->GetEntriesFast(); ++i){
            clusterInModule["cern"] = clCreators["cern"].FindClustersInEventBoole(*(data["cern"]->at(i)), 1.5, 2.5, 4.0, 100, false);
            clusterInModule["slayer"] = clCreators["slayer"].FindClustersInEventBoole(*(data["slayer"]->at(i)), 1.5, 2.5, 4.0, 100, false);
            clusterInModule["HD2"] = clCreators["HD2"].FindClustersInEventBoole(*(data["HD2"]->at(i)), 1.5, 2.5, 4.0, 100, false);
            clusterInModule["HD1"] = clCreators["HD1"].FindClustersInEventBoole(*(data["HD1"]->at(i)), 1.5, 2.5, 4.0, 100, false);

            if(clusterInModule["cern"] && clusterInModule["HD2"] && clusterInModule["HD1"]) ++missedEvents;

            if(clusterInModule["cern"] && clusterInModule["slayer"] && clusterInModule["HD2"] && clusterInModule["HD1"]){ // if cluster in all modules store the slayer one for analysis
                clCreators["simulation"].FindClustersInEventBoole((*data["slayer"]->at(i)), 1.5, 2.5, 4.0, 100, false);
                ++foundEvents;
            }

        }
    }

    if(c.simulation){
        std::cout << "Found " << clCreators["simulation"].getNumberOfClusters() << " clusters in " << inputTree->GetEntriesFast() << " events!\n";
        std::cout << "Missed " << missedEvents << " events\n";
        double event2OneOrMoreClusterEff = 1. - (double)missedEvents/data["simulation"]->size();
        std::cout << "Rate of MCEvent producing one or more clusters: " << event2OneOrMoreClusterEff << " +/- "
                  << TMath::Sqrt( (1-event2OneOrMoreClusterEff)*event2OneOrMoreClusterEff/(double)data["simulation"]->size() )  << "\n";
    }
    else{
      std::cout << "Found " << missedEvents << " events that caused one or more cluster in all three but the slayer modules.\n";
      std::cout << "Found " << foundEvents << " events that caused one or more clusters in all modules\n";
      std::cout << "This is the fraction of " << foundEvents / missedEvents << "\n";
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
    stdDevLightYield = TMath::Sqrt(stdDevLightYield);

    std::cout << "Mean light yield: " << meanLightYield << " +/- " << stdDevLightYield << "\n";


    if (c.clusterAlg == "m") c.tag += "_max";

    ClusterMonitor clMonitor;
    clMonitor.WriteToNtuple(clCreators["simulation"], ("/home/tobi/SciFi/results/clusters/" + removePath(file2analyse).ReplaceAll(".root", "_clusterAnalyis" + c.tag +".root")).Data() );

  }

  return 0;
}
