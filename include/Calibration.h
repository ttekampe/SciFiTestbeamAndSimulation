#ifndef CALIBRATION_H
#define CALIBRATION_H

// from std
#include <map>
#include <string>
#include <vector>

// from ROOT
#include "TTree.h"

// from here
#include "Cluster.h"

typedef std::vector<Channel> Event;

struct CalibrationFiles {
  std::string pedestal;
  std::string gain;
};

struct calibrationRunNumbers {
  unsigned int dark;
  unsigned int led;
};

typedef std::map<unsigned int, double> map_uint_double;
typedef std::map<unsigned int, map_uint_double> map_uint_map_uint_double;

std::vector<Event> parseRootTree(
    TTree* dataTree, unsigned int uplinkMin, unsigned int uplinkMax,
    unsigned int nAdcs, double factor = 1.,
    const map_uint_map_uint_double& pedestals = map_uint_map_uint_double(),
    const map_uint_map_uint_double& gains = map_uint_map_uint_double());

// std::vector<Event> parseCorrectedRootTree(TTree* dataTree,
//                                           unsigned int uplinkMin,
//                                           unsigned int uplinkMax,
//                                           unsigned int nAdcs,
//                                           double factor = 1.);

void produceGains(TTree* t, const unsigned int uplinkMin,
                  const unsigned int uplinkMax, const unsigned int adcIDmin,
                  const unsigned int adcIDmax, const unsigned int maxGaussians,
                  TString fileName, std::string savePath);

map_uint_map_uint_double readGains(std::string fileName);

calibrationRunNumbers lookUpCalibrationFiles(
    const unsigned int runNumber, const std::string catalogueFileName);

unsigned int runNumberFromFilename(std::string filename);

map_uint_map_uint_double getPedestals(const std::string fileName,
                                      const unsigned int uplinkMin,
                                      const unsigned int uplinkMax,
                                      const unsigned int nAdcs);

void correctFile(TTree* tree2correct, const map_uint_map_uint_double& gains,
                 const map_uint_map_uint_double& pedestals,
                 const unsigned int uplinkMin, const unsigned int uplinkMax,
                 const unsigned int nAdcs, TString newFileName);

TString removePath(TString str);

#endif
