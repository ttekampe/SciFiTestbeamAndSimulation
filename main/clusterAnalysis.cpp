// from std
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <vector>

// from ROOT
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2I.h"
#include "TMath.h"
#include "TTree.h"

// from RooFit
#include "RooDataSet.h"
#include "RooFit.h"
#include "RooFitResult.h"
#include "RooGaussian.h"
#include "RooPlot.h"
#include "RooRealVar.h"

// from here
//#include "Cluster.h"
//#include "ClusterAlgorithms.h"
#include "Calibration.h"
#include "ClusterCreator.h"
#include "ClusterMonitor.h"
#include "EDouble.h"

// from BOOST
#include "boost/filesystem.hpp"
#include "boost/program_options.hpp"

namespace po = boost::program_options;

struct config {
  std::vector<std::string> files2analyse;
  bool debug;
  bool simulation;
  std::string clusterAlg;
  std::string tag;
  std::string outputDir;
};

int parseOptions(config &c, int argc, char *argv[]) {
  // declare options
  po::options_description desc("Allowed options");
  desc.add_options()("help", "show this help")(
      "file,f",
      po::value<std::vector<std::string>>(&c.files2analyse)->multitoken(),
      "corrected test beam data file")(
      "simulation,s", po::bool_switch(&c.simulation), "Simulated input?")
      // ("clusteralg,c",
      // po::value<std::string>(&c.clusterAlg)->default_value("b"), "clustering
      // algorithm: b for Boole or m for Maxs")
      ("tag,t", po::value<std::string>(&c.tag)->default_value(""),
       "tag that is added to the output file name")(
          "debug,d", po::bool_switch(&c.debug), "debug output")(
          "odir,o",
          po::value<std::string>(&c.outputDir)
              ->default_value(std::string(std::getenv("HOME")) + "/SciFi"),
          "directory the output is written to");

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

std::string getPositionFromFileName(std::string fileName) {
  // maps testbeam run numbers to the position in testbeam jargon:
  // a: beam hits module close to mirror, b: middle, c: close to scipm
  // in the case of the attenuation scan, the measured distance is returned
  // see https://twiki.cern.ch/twiki/bin/view/LHCb/SciFiTrackerTestBeam2015
  std::string position;
  if (fileName.find("btsoftware") != std::string::npos) {  // data
    std::map<std::string, std::string> runNumbersAndPositions;
    runNumbersAndPositions["1431786652"] = "a";
    runNumbersAndPositions["1432091294"] = "a";
    runNumbersAndPositions["1432264510"] = "a";
    runNumbersAndPositions["1432350729"] = "a";
    runNumbersAndPositions["1432341050"] = "a";
    runNumbersAndPositions["1432169457"] = "c";
    runNumbersAndPositions["1432089102"] = "c";
    runNumbersAndPositions["1432187205"] = "c";
    runNumbersAndPositions["1432358809"] = "c";

    // AttScan
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

    std::string currentRunNumber =
        std::to_string(runNumberFromFilename(fileName));
    position = runNumbersAndPositions.find(currentRunNumber) !=
                       runNumbersAndPositions.end()
                   ? runNumbersAndPositions[currentRunNumber]
                   : "";

  } else {  // simulation
    std::string lookFor = "testbeam_simulation_position_";
    auto posStartAt = fileName.find(lookFor) + lookFor.size();
    auto posEndAt = fileName.find("_", posStartAt);
    position = fileName.substr(posStartAt, posEndAt - posStartAt);
    std::cout << "Current position is " << position << "\n";
  }
  return position;
}

std::pair<EDouble, EDouble> analyse(std::string file2analyse, const config &c,
                                    bool recursive_call = false) {
  TFile inputFile(file2analyse.c_str(), "READ");
  if (!inputFile.IsOpen()) {
    std::cout << "Could not open " << file2analyse << "\n";
    return std::make_pair(EDouble(0, 0), EDouble(0, 0));
  }

  std::cout << "Analysis " << file2analyse << "\n";

  TTree *inputTree;
  if (c.simulation) {
    inputTree = dynamic_cast<TTree *>(inputFile.Get("layer_0"));
  } else {
    inputTree = dynamic_cast<TTree *>(inputFile.Get("rawData"));
  }

  if (!inputTree) {
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

  std::map<std::string, std::vector<std::vector<Channel>>> data;

  // offsets between the modules perpendicular to the beam
  std::vector<double> xPositions;
  boost::filesystem::create_directories(c.outputDir + "/results/moduleOffset/");
  std::string offsetFileName =
      (c.outputDir + "/results/moduleOffset/" +
       removePath(file2analyse).ReplaceAll(".root", ".txt"))
          .Data();
  bool produceOffsetFile{false};
  std::vector<double> xOffsets;
  std::vector<double> track_distances;
  std::vector<double> track_dxdz;

  // for data the offset needs to be known to measure the distance
  // between the "track" and a found cluster
  if (!c.simulation) {
    std::cout << "Checking for offsetFile...\n";
    std::ifstream offsetFile(offsetFileName);
    if (!offsetFile) {
      if (recursive_call) {
        std::cerr
            << "Already runnig a second time, still no offset file\n"
            << "Maybe you do not have write permissions at your output dir\n"
            << "Giving up";
        throw;
      }
      std::cout << "Offsetfile does not exists yet -> producing...\n";
      produceOffsetFile = true;
    } else {
      std::string line;
      std::getline(offsetFile, line);
      std::istringstream ss(line);
      double offset;
      if (ss >> offset) xOffsets.push_back(offset);
      std::cout << "Found xoffset: " << xOffsets[0] << "\n";
    }
  }

  // distance between modules in beam direction in the order HD2, Slayer, CERN
  std::vector<double> zPositions = {0., 247.0 * 1000, 469.0 * 1000};

  if (c.simulation) {
    data["simulation"] = parseRootTree(inputTree, 1, 4, 128);
  } else {
    data["cern"] = parseRootTree(inputTree, 1, 2, 128);
    data["slayer"] = parseRootTree(inputTree, 3, 4, 128);
    data["HD2"] = parseRootTree(inputTree, 5, 6, 128);
  }

  std::map<std::string, ClusterCreator> clCreators;
  if (!c.simulation) {
    clCreators["cern"] = ClusterCreator();
    clCreators["slayer"] = ClusterCreator();
    clCreators["HD2"] = ClusterCreator();
  }

  // in the case of data this one stores the matched clusters
  clCreators["simulation"] = ClusterCreator();

  std::map<std::string, std::vector<Cluster>> clustersInModule;
  clustersInModule["cern"];
  clustersInModule["slayer"];
  clustersInModule["HD2"];

  ClusterCreator cl2Analyse;
  unsigned int currentNumberOfClusters{0};
  double missedEvents{0};
  double foundEvents{0};
  if (c.simulation) {
    for (const auto &event : data["simulation"]) {
      clCreators["simulation"].FindClustersInEventBoole(event, 1.5, 2.5, 4.0,
                                                        100);
      if (currentNumberOfClusters ==
          clCreators["simulation"].getNumberOfClusters())
        ++missedEvents;
      currentNumberOfClusters = clCreators["simulation"].getNumberOfClusters();
    }
  } else {
    for (unsigned int i = 0; i < inputTree->GetEntriesFast(); ++i) {
      clustersInModule["cern"] = clCreators["cern"].FindClustersInEventBoole(
          (data["cern"].at(i)), 1.5, 2.5, 4.0, 100);
      clustersInModule["slayer"] =
          clCreators["slayer"].FindClustersInEventBoole((data["slayer"].at(i)),
                                                        1.5, 2.5, 4.0, 100);
      clustersInModule["HD2"] = clCreators["HD2"].FindClustersInEventBoole(
          (data["HD2"].at(i)), 1.5, 2.5, 4.0, 100);

      if (clustersInModule["cern"].size() == 1 &&
          clustersInModule["slayer"].size() == 1 &&
          clustersInModule["HD2"].size() ==
              1) {  // if there is a cluster in all
                    // modules, keep the one in the
                    // slayer module

        // why am I calling this again? Doesnt matter it's fast
        clustersInModule["slayer"] =
            clCreators["simulation"].FindClustersInEventBoole(
                (data["slayer"].at(i)), 1.5, 2.5, 4.0, 100);

        // get the x positions of the clusters to do some primitive tracking
        xPositions = {
            clustersInModule["HD2"][0].GetChargeWeightedMean() * 250.,
            clustersInModule["slayer"][0].GetChargeWeightedMean() * 250.,
            clustersInModule["cern"][0].GetChargeWeightedMean() * 250.};

        // just calculate the equation for a strright line through the cluster
        // positions in the HD2 and the CERN modules
        double dx = xPositions[2] - xPositions[0];
        double dz = zPositions[2] - zPositions[0];

        double slope = dx / dz;
        double constant = xPositions[0];

        // calculate the position of the track at the slayer module
        double trackAtSlayer = constant + slope * zPositions[1];
        // store the (not yet calibrated) offset in the vector for the offset
        // file or if it already exists store the distance between the track
        // and the cluster
        if (produceOffsetFile)
          xOffsets.push_back(trackAtSlayer - xPositions[1]);
        else {
          track_distances.push_back(trackAtSlayer - xPositions[1] -
                                    xOffsets[0]);
          track_dxdz.push_back(slope);
        }
      }

      // count found and missed events for the efficiency calculation
      if (!clustersInModule["cern"].empty() &&
          !clustersInModule["HD2"].empty()) {
        ++missedEvents;
      }
      if (!clustersInModule["cern"].empty() &&
          !clustersInModule["slayer"].empty() &&
          !clustersInModule["HD2"].empty()) {
        ++foundEvents;
      }
    }
  }

  if (produceOffsetFile) {
    // fit gaussian to offsets
    TCanvas can_offset;

    std::cout << "Collected " << xOffsets.size() << " xoffsets\n";

    double sum = std::accumulate(xOffsets.begin(), xOffsets.end(), 0.0);
    double mean = sum / xOffsets.size();

    double sq_sum = std::inner_product(xOffsets.begin(), xOffsets.end(),
                                       xOffsets.begin(), 0.0);
    double stdev = std::sqrt(sq_sum / xOffsets.size() - mean * mean);

    std::cout << "mean: " << mean << " stdev: " << stdev << "\n";

    RooRealVar rv_xoffset("rv_xoffset", "", mean - 300, mean + 300);

    RooDataSet dataset_xoffset("dataset_xoffset", "", RooArgSet(rv_xoffset));
    for (const auto xoffset : xOffsets) {
      // veto absurd offset values (noise / wrongly associated clusters)
      if (xoffset > mean + 300 || xoffset < mean - 300) continue;
      rv_xoffset.setVal(xoffset);
      dataset_xoffset.add(RooArgSet(rv_xoffset));
    }
    dataset_xoffset.Print();
    RooRealVar rv_mean("rv_mean", "", mean, mean - 100, mean + 100);
    RooRealVar rv_sigma("rv_sigma", "", 50, 0, 100);
    RooGaussian gauss("gauss", "", rv_xoffset, rv_mean, rv_sigma);

    RooFitResult *fr = gauss.fitTo(dataset_xoffset, RooFit::Save(true),
                                   RooFit::Minimizer("Minuit2", "minimize"));

    RooPlot *plot = rv_xoffset.frame();
    dataset_xoffset.plotOn(plot);
    gauss.plotOn(plot);
    plot->Draw();

    boost::filesystem::create_directories(c.outputDir +
                                          "/results/moduleOffset/");
    can_offset.SaveAs(c.outputDir + "/results/moduleOffset/" +
                      removePath(file2analyse).ReplaceAll(".root", ".pdf"));

    fr->Print("v");

    std::ofstream offsetFile(
        (c.outputDir + "/results/moduleOffset/" +
         removePath(file2analyse).ReplaceAll(".root", ".txt"))
            .Data());
    // store the found offset in the file
    offsetFile << rv_mean.getVal() << "\n";
    offsetFile.close();

    // restart the procedure, this time the offsetfile will exist

    delete inputTree;
    inputFile.Close();
    return analyse(file2analyse, c, true);
  }

  double sumLightYield = std::accumulate(
      clCreators["simulation"].getClusters().begin(),
      clCreators["simulation"].getClusters().end(), 0.,
      [](double ly, const Cluster &cl) { return ly + cl.GetSumOfAdcValues(); });

  double sumLightYieldSq = std::accumulate(
      clCreators["simulation"].getClusters().begin(),
      clCreators["simulation"].getClusters().end(), 0.,
      [](double ly, const Cluster &cl) {
        return (ly + cl.GetSumOfAdcValues() * cl.GetSumOfAdcValues());
      });

  double meanLightYield =
      sumLightYield / (double)clCreators["simulation"].getNumberOfClusters();
  double stdDevLightYield = std::sqrt(
      sumLightYieldSq / (double)clCreators["simulation"].getNumberOfClusters() -
      meanLightYield * meanLightYield);

  stdDevLightYield /=
      std::sqrt((double)clCreators["simulation"].getNumberOfClusters());

  std::cout << "Mean light yield: " << meanLightYield << " +/- "
            << stdDevLightYield << "\n";

  EDouble light;
  EDouble eff;

  if (c.simulation) {
    std::cout << "Found " << clCreators["simulation"].getNumberOfClusters()
              << " clusters in " << inputTree->GetEntriesFast() << " events!\n";
    std::cout << "Missed " << missedEvents << " events\n";
    double event2OneOrMoreClusterEff =
        1. - (double)missedEvents / data["simulation"].size();
    double event2OneOrMoreClusterEffErr = TMath::Sqrt(
        (1 - event2OneOrMoreClusterEff) * event2OneOrMoreClusterEff /
        (double)data["simulation"].size());
    std::cout << "Rate of MCEvent producing one or more clusters: "
              << event2OneOrMoreClusterEff << " +/- "
              << event2OneOrMoreClusterEffErr << "\n";

    light = EDouble(meanLightYield, stdDevLightYield);
    eff = EDouble(event2OneOrMoreClusterEff, event2OneOrMoreClusterEffErr);
    // hitEffFile << getPositionFromFileName(file2analyse) << "," <<
    // meanLightYield << "," << stdDevLightYield << ","
    //     << event2OneOrMoreClusterEff << "," <<  event2OneOrMoreClusterEffErr
    //     << "\n";
  } else {
    std::cout << "Found " << missedEvents << " events that caused one or more "
                                             "cluster in all three but the "
                                             "slayer modules.\n";
    std::cout << "Found " << foundEvents
              << " events that caused one or more clusters in all modules\n";
    double efficiency = foundEvents / missedEvents;
    double effErr = TMath::Sqrt(efficiency * (1. - efficiency) / missedEvents);
    std::cout << "This is the fraction of " << foundEvents / missedEvents
              << "\n";

    // hitEffFile << "data" << "," << meanLightYield << "," << stdDevLightYield
    // << ","
    //      << efficiency << "," <<  effErr << "\n";

    light = EDouble(meanLightYield, stdDevLightYield);
    eff = EDouble(efficiency, effErr);
  }

  // if (c.clusterAlg == "m") c.tag += "_max";

  std::map<std::string, std::vector<double>> features;
  features["distance_from_track"] = track_distances;
  features["track_dxdz"] = track_dxdz;

  ClusterMonitor clMonitor;
  boost::filesystem::create_directories(c.outputDir + "/results/clusters/");
  clMonitor.WriteToNtuple(
      clCreators["simulation"],
      (c.outputDir + "/results/clusters/" +
       removePath(file2analyse)
           .ReplaceAll(".root", "_clusterAnalyis" + c.tag + ".root"))
          .Data(),
      features);
  delete inputTree;
  inputFile.Close();

  return std::make_pair(light, eff);
}

int main(int argc, char *argv[]) {
  config c;
  if (parseOptions(c, argc, argv) != 0) {
    std::cerr << "Can not parse options!" << std::endl;
    return 0;
  }

  std::ofstream hitEffFile(c.outputDir + "/results/hitEfficiency.txt");
  hitEffFile
      << "position\tlightyield\tlightyieldErr\tefficiency\tefficiencyErr\n";

  std::pair<EDouble, EDouble> lightAndEff;

  for (const auto &file2analyse : c.files2analyse) {
    lightAndEff = analyse(file2analyse, c);
    hitEffFile << getPositionFromFileName(file2analyse) << ","
               << lightAndEff.first.GetVal() << ","
               << lightAndEff.first.GetError() << ","
               << lightAndEff.second.GetVal() << ","
               << lightAndEff.second.GetError() << "\n";
  }

  return 0;
}
