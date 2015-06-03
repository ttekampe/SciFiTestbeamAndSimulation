#ifndef CALIBRATION_H
#define CALIBRATION_H

//from std
#include <map>
#include <vector>
#include <string>

//from ROOT
#include "TTree.h"

//from here
#include "Cluster.h"


struct CalibrationFiles{
  std::string pedestal;
  std::string gain;
};

struct calibrationRunNumbers{
  unsigned int dark;
  unsigned int led;
};

std::vector<std::vector<Channel>*>* parseRootTree(TTree* dataTree,
                                                  unsigned int uplinkMin,
                                                  unsigned int uplinkMax,
                                                  unsigned int nAdcs,
                                                  const std::map<unsigned int, std::map<unsigned int, double>>& pedestals,
                                                  const std::map<unsigned int, std::map<unsigned int, double>>& gains);



void produceGains(TTree* t,
                  const unsigned int uplinkMin,
                  const unsigned int uplinkMax,
                  const unsigned int adcIDmin,
                  const unsigned int adcIDmax,
                  const unsigned int maxGaussians,
                  TString fileName,
                  std::string savePath);


std::map<unsigned int, std::map<unsigned int, double>> readGains(std::string fileName);

calibrationRunNumbers lookUpCalibrationFiles(const unsigned int runNumber);

unsigned int runNumberFromFilename(std::string filename);

std::map<unsigned int, std::map<unsigned int, double> > getPedestals(const std::string fileName,
                                                                     const unsigned int uplinkMin,
                                                                     const unsigned int uplinkMax,
                                                                     const unsigned int nAdcs);

#endif
