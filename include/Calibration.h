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

typedef std::vector<Channel> Event;

struct CalibrationFiles{
  std::string pedestal;
  std::string gain;
};

struct calibrationRunNumbers{
  unsigned int dark;
  unsigned int led;
};

std::vector<Event*>* parseRootTree(TTree* dataTree,
                                  unsigned int uplinkMin,
                                  unsigned int uplinkMax,
                                  unsigned int nAdcs,
                                  const std::map<unsigned int, std::map<unsigned int, double>>& pedestals,
                                  const std::map<unsigned int, std::map<unsigned int, double>>& gains);

std::vector<Event*>* parseCorrectedRootTree(TTree* dataTree,
                                            unsigned int uplinkMin,
                                            unsigned int uplinkMax,
                                            unsigned int nAdcs,
                                            double factor = 1.);


void produceGains(TTree* t,
                  const unsigned int uplinkMin,
                  const unsigned int uplinkMax,
                  const unsigned int adcIDmin,
                  const unsigned int adcIDmax,
                  const unsigned int maxGaussians,
                  TString fileName,
                  std::string savePath);


std::map<unsigned int, std::map<unsigned int, double>> readGains(std::string fileName);

calibrationRunNumbers lookUpCalibrationFiles(const unsigned int runNumber, const std::string catalogueFileName);

unsigned int runNumberFromFilename(std::string filename);

std::map<unsigned int, std::map<unsigned int, double> > getPedestals(const std::string fileName,
                                                                     const unsigned int uplinkMin,
                                                                     const unsigned int uplinkMax,
                                                                     const unsigned int nAdcs);

void correctFile(TTree* tree2correct,
                 const std::map<unsigned int, std::map<unsigned int, double>>& gains,
                 const std::map<unsigned int, std::map<unsigned int, double> >& pedestals,
                 const unsigned int uplinkMin,
                 const unsigned int uplinkMax,
                 const unsigned int nAdcs,
                 TString newFileName);

TString removePath(TString str);

#endif
