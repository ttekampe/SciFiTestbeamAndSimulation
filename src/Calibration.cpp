// from std
#include <algorithm>
#include <fstream>
#include <iostream>
#include <memory>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// from ROOT
#include "RooAbsReal.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooGaussian.h"
#include "RooMinuit.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TSpectrum.h"
#include "TString.h"
#include "TTree.h"

// from here
#include "Calibration.h"
#include "Cluster.h"
#include "lhcbStyle.h"

TString removePath(TString str) {
  Ssiz_t posOfSlash = str.Last('/');
  if (posOfSlash != kNPOS) {
    str = str.Replace(0, posOfSlash + 1, "");
  }
  return str;
}

std::vector<Event> parseRootTree(TTree *dataTree, unsigned int uplinkMin,
                                 unsigned int uplinkMax, unsigned int nAdcs,
                                 double factor,
                                 const map_uint_map_uint_double &pedestals,
                                 const map_uint_map_uint_double &gains) {
  std::vector<Event> dataVector(dataTree->GetEntries());

  if (factor != 1.)
    std::cout << "Parsing corrected root tree and correct the adc value for "
                 "the factor of "
              << factor << "!\n";

  bool perform_correction{false};
  if (!pedestals.empty() && !gains.empty()) {
    std::cout << "Correcting for gain and pedestal\n";
    perform_correction = true;
  }

  const unsigned int nUplinks = uplinkMax - uplinkMin + 1;
  std::vector<std::vector<float>> adcVals(nUplinks, std::vector<float>(nAdcs));

  dataTree->SetBranchStatus("*", 0);
  for (unsigned int uplink = uplinkMin; uplink <= uplinkMax; ++uplink) {
    std::string branchName = "Uplink_" + std::to_string(uplink) + "_adc_";
    for (unsigned int adc = 0; adc < nAdcs; ++adc) {
      dataTree->SetBranchStatus((branchName + std::to_string(adc + 1)).c_str(),
                                1);
      dataTree->SetBranchAddress((branchName + std::to_string(adc + 1)).c_str(),
                                 &adcVals[uplink - uplinkMin][adc]);
    }
  }
  std::cout << "channels per event: " << (nAdcs * (uplinkMax - uplinkMin + 1))
            << "\n";
  for (unsigned int i = 0; i < dataTree->GetEntriesFast(); ++i) {
    dataTree->GetEntry(i);
    Event event(nAdcs * (uplinkMax - uplinkMin + 1));
    for (unsigned int uplink = uplinkMin; uplink <= uplinkMax; ++uplink) {
      for (unsigned int adc = 0; adc < nAdcs; ++adc) {
        Channel c;
        c.Uplink = uplink;
        c.ChannelNumber = adc + 1;
        if (adcVals[uplink - uplinkMin][adc] > 0) {
          c.AdcValue = adcVals[uplink - uplinkMin][adc] * factor;
        } else
          c.AdcValue = 0;
        if (perform_correction) {
          c.AdcValue = (c.AdcValue - pedestals.at(uplink).at(adc + 1)) /
                       gains.at(uplink).at(adc + 1);
        }
        event[adc + (uplink - uplinkMin) * nAdcs] = std::move(c);
      }
    }
    dataVector[i] = std::move(event);
  }

  //  dataTree->ResetBranchAddresses();

  return dataVector;
}

void produceGains(TTree *t, const unsigned int uplinkMin,
                  const unsigned int uplinkMax, const unsigned int adcIDmin,
                  const unsigned int adcIDmax, const unsigned int maxGaussians,
                  TString fileName, std::string savePath) {
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
  RooMsgService::instance().setSilentMode(true);
  lhcb::lhcbStyle();
  TString saveName = fileName;
  saveName.Remove(0, saveName.Last('/') + 1);
  saveName.ReplaceAll(".root", "_gain.txt");
  std::ofstream gainFile(savePath + saveName);
  gainFile << "uplinkNumber\tadcNumber\tgain\tgainErr" << std::endl;

  TString canvasName = saveName;
  canvasName.ReplaceAll(".txt", ".pdf");

  TCanvas can("", "");
  can.SaveAs((savePath + canvasName + "[").Data());

  double numFits{0.};
  double failedFits{0.};

  TString branchNameTemplate = "Uplink_";
  for (unsigned int uplinkNumber = uplinkMin; uplinkNumber <= uplinkMax;
       ++uplinkNumber) {
    TString branchNameTemplate2 =
        branchNameTemplate + std::to_string(uplinkNumber) + "_adc_";

    for (unsigned int adcNumber = adcIDmin; adcNumber <= adcIDmax;
         ++adcNumber) {
      unsigned int nGaussians = maxGaussians;
      TString branchName = branchNameTemplate2 + std::to_string(adcNumber);
      t->SetBranchStatus("*", 0);
      t->SetBranchStatus(branchName, 1);
      TSpectrum spec(nGaussians);
      TH1D hist("hist", "hist", 120, 400., 1000.);
      t->Draw(branchName + ">>hist");
      spec.Search(&hist);
      nGaussians = spec.GetNPeaks();
      auto peaks = spec.GetPositionX();

      auto max_peak = *std::max_element(peaks, peaks + nGaussians);
      auto min_peak = *std::min_element(peaks, peaks + nGaussians);

      RooRealVar adcCount(branchName, branchName, min_peak - 30.,
                          max_peak + 30.);
      RooRealVar mean("mean", "", min_peak, 420., 650.);
      RooRealVar sigma1("sigma1", "", 7., 0., 20.);
      RooRealVar sigma2("sigma2", "", 4., 0., 20.);
      RooRealVar gain("gain", "", (max_peak - min_peak) / (nGaussians - 1.),
                      45., 75.);
      RooDataSet dataSet("dataSet", "", t, RooArgSet(adcCount));

      std::vector<std::shared_ptr<RooGaussian>> gaussians(nGaussians, nullptr);
      std::vector<std::shared_ptr<RooFormulaVar>> means(nGaussians, nullptr);
      std::vector<std::shared_ptr<RooFormulaVar>> sigmas(nGaussians, nullptr);
      std::vector<std::shared_ptr<RooRealVar>> fractions(nGaussians - 1,
                                                         nullptr);

      RooArgList gaussianList("gaussianList");
      RooArgList fractionList("fractionList");
      for (unsigned int j = 0; j < nGaussians; ++j) {
        means[j] = std::make_shared<RooFormulaVar>(
            TString("mean_" + std::to_string(j)).Data(), "",
            TString("mean+gain*" + std::to_string(j)).Data(),
            RooArgList(mean, gain));
        sigmas[j] = std::make_shared<RooFormulaVar>(
            TString("sigma_" + std::to_string(j)).Data(), "",
            TString("sqrt(sigma1*sigma1+" + std::to_string(j) +
                    "*sigma2*sigma2)")
                .Data(),
            RooArgList(sigma1, sigma2));
        gaussians[j] = std::make_shared<RooGaussian>(
            TString("gaussian_" + std::to_string(j)).Data(), "", adcCount,
            *means[j], *sigmas[j]);
        gaussianList.add(*gaussians[j]);
        if (j > 0) {
          fractions[j - 1] = std::make_shared<RooRealVar>(
              TString("fraction_" + std::to_string(j)).Data(), "", 1. / 2. / j,
              0., 1.);
          fractionList.add(*fractions[j - 1]);
        }
      }
      RooAddPdf pdf("pdf", "", gaussianList, fractionList);

      // RooFitResult* fs = pdf.fitTo(dataSet, RooFit::Save(true));
      // fs->Print("v");
      ++numFits;
      RooAbsReal *nll = pdf.createNLL(dataSet);

      RooMinuit minu(*nll);
      minu.setStrategy(2);  // ensure minimum true, errors correct

      // call MIGRAD -- minimises the likelihood
      int migradStatusCode = minu.migrad();  // catch status code

      // could also get an intermediary RooFitResult if you want
      //      RooFitResult *migradFitResult = m.save();
      //      int migradStatusCode2 = migradFitResult->status();

      // call HESSE -- calculates error matrix
      int hesseStatusCode = minu.hesse();

      RooFitResult *fs =
          minu.save();  // equivalent result to that from fitTo call

      // check for success
      int covarianceQuality = fs->covQual();

      //    int minosStatusCode = minu.minos();

      bool isAGoodFit =
          ((hesseStatusCode == 0) && (migradStatusCode == 0) &&
           (covarianceQuality == 3));  // && (minosStatusCode==0) );

      if (!isAGoodFit || gain.getError() == 0.0) {
        ++failedFits;
        std::cout << "Fit for uplink " << uplinkNumber << " and adc "
                  << adcNumber << " failed!" << std::endl;
        fs->Print("v");
        gainFile << uplinkNumber << "\t" << adcNumber << "\t" << 0 << "\t" << 0
                 << std::endl;
      } else
        gainFile << uplinkNumber << "\t" << adcNumber << "\t" << gain.getVal()
                 << "\t" << gain.getError() << std::endl;

      //      fs->Print("v");

      RooPlot *myFrame = adcCount.frame();
      dataSet.plotOn(myFrame);
      pdf.plotOn(myFrame);
      myFrame->Draw();
      can.SaveAs((savePath + canvasName).Data());
    }
  }
  can.SaveAs((savePath + canvasName + "]").Data());

  std::cout << failedFits / numFits * 100. << "% of all fits failed."
            << std::endl;
}

map_uint_map_uint_double readGains(std::string fileName) {
  std::ifstream inputFile(fileName);
  std::string line;
  int uplinkNumber{0};
  int adcNumber{0};
  double gain{0};

  map_uint_map_uint_double gains;
  while (std::getline(inputFile, line)) {
    std::istringstream ss(line);
    if (ss >> uplinkNumber >> adcNumber >> gain) {
      std::cout << uplinkNumber << "\t" << adcNumber << "\t" << gain
                << std::endl;

      gains[uplinkNumber][adcNumber] = gain;
    }
  }
  return gains;
}

calibrationRunNumbers lookUpCalibrationFiles(
    const unsigned int runNumber, const std::string catalogueFileName) {
  //  std::ifstream fileCatalogue("/data/testbeam/data/runNumbers.txt");
  std::ifstream fileCatalogue(catalogueFileName);
  unsigned int currentRunNumber = 0;
  unsigned int currentLedNumber = 0;
  unsigned int currentDarkNumber = 0;
  std::string line = "";
  calibrationRunNumbers rv;
  while (std::getline(fileCatalogue, line)) {
    std::istringstream ss(line);
    ss >> currentDarkNumber >> currentLedNumber >> currentRunNumber;
    if (currentRunNumber == runNumber) {
      rv.dark = currentDarkNumber;
      rv.led = currentLedNumber;
      return rv;
    }
  }
  throw std::runtime_error(
      "calibration::lookUpCalibrationFiles: Unable to "
      "find a match in run number catalogue.");
}

unsigned int runNumberFromFilename(std::string filename) {
  filename = removePath(filename);
  // In May15 testbeam, unixtime was used as runNumber and contained 10 digits
  std::regex re("\\d{10,}");
  std::smatch match;
  std::regex_search(filename, match, re);
  return std::stoi(match[0]);
}

map_uint_map_uint_double getPedestals(const std::string fileName,
                                      const unsigned int uplinkMin,
                                      const unsigned int uplinkMax,
                                      const unsigned int nAdcs) {
  TFile pedestalFile(fileName.c_str(), "READ");
  if (!pedestalFile.IsOpen()) {
    std::cerr << "unable to open pedestal file " << fileName << std::endl;
  }
  TTree *pedestalTree = dynamic_cast<TTree *>(pedestalFile.Get("rawData"));
  if (pedestalTree == nullptr) {
    std::cerr << "pedestal tree is nullptr" << std::endl;
  }
  const unsigned int nUplinks = uplinkMax - uplinkMin + 1;
  std::cout << "creating double array of form double[" << nUplinks << "]["
            << nAdcs << "]\n";

  std::vector<std::vector<float>> rawAdcVals(nUplinks,
                                             std::vector<float>(nAdcs));

  std::cout << "setting branch addresses" << std::endl;
  std::string branchNameTemplate1 = "Uplink_";
  for (unsigned int uplink = uplinkMin; uplink <= uplinkMax; ++uplink) {
    std::string branchNameTemplate2 =
        branchNameTemplate1 + std::to_string(uplink) + "_adc_";
    for (unsigned int adc = 1; adc <= nAdcs; ++adc) {
      std::string branchName = branchNameTemplate2 + std::to_string(adc);
      std::cout << "setting adress for branch " << branchName << " to double["
                << uplink - uplinkMin << "][" << adc - 1 << "]" << std::endl;
      pedestalTree->SetBranchAddress(branchName.c_str(),
                                     &rawAdcVals[uplink - uplinkMin][adc - 1]);
    }
  }
  std::cout << "loop over file" << std::endl;
  map_uint_map_uint_double pedestals;
  for (unsigned int i = 0; i < pedestalTree->GetEntriesFast(); ++i) {
    pedestalTree->GetEntry(i);
    for (unsigned int uplink = uplinkMin; uplink <= uplinkMax; ++uplink) {
      for (unsigned int adc = 1; adc <= nAdcs; ++adc) {
        pedestals[uplink][adc] += rawAdcVals[uplink - uplinkMin][adc - 1];
      }
    }
  }

  for (unsigned int uplink = uplinkMin; uplink <= uplinkMax; ++uplink) {
    for (unsigned int adcNum = 1; adcNum <= nAdcs; ++adcNum) {
      pedestals[uplink][adcNum] /= pedestalTree->GetEntriesFast();
    }
  }
  pedestalTree->ResetBranchAddresses();
  pedestalFile.Close();

  return pedestals;
}

void correctFile(TTree *tree2correct, const map_uint_map_uint_double &gains,
                 const map_uint_map_uint_double &pedestals,
                 const unsigned int uplinkMin, const unsigned int uplinkMax,
                 const unsigned int nAdcs, TString newFileName) {
  //  newFileName = removePath(newFileName);
  //  TFile correctedFile(("/data/testbeam/data/corrected/" +
  //  newFileName).Data(), "RECREATE");

  TFile correctedFile(newFileName.Data(), "RECREATE");
  auto correctedTree = new TTree("rawData", "rawData");

  // setup fill mechanism for new TTree
  const unsigned int nUplinks = uplinkMax - uplinkMin + 1;
  std::vector<std::vector<float>> rawAdcVals(nUplinks,
                                             std::vector<float>(nAdcs));
  std::vector<std::vector<float>> correctedAdcVals(nUplinks,
                                                   std::vector<float>(nAdcs));

  std::string branchNameTemplate1 = "Uplink_";
  for (unsigned int uplink = uplinkMin; uplink <= uplinkMax; ++uplink) {
    std::string branchNameTemplate2 =
        branchNameTemplate1 + std::to_string(uplink) + "_adc_";
    for (unsigned int adc = 0; adc < nAdcs; ++adc) {
      std::string branchName = branchNameTemplate2 + std::to_string(adc + 1);
      tree2correct->SetBranchAddress(branchName.c_str(),
                                     &rawAdcVals[uplink - uplinkMin][adc]);
      correctedTree->Branch(branchName.c_str(),
                            &correctedAdcVals[uplink - uplinkMin][adc],
                            (branchName + "/F").c_str());
    }
  }
  for (unsigned int i = 0; i < tree2correct->GetEntriesFast(); ++i) {
    tree2correct->GetEntry(i);
    for (unsigned int uplink = uplinkMin; uplink <= uplinkMax; ++uplink) {
      for (unsigned int adc = 1; adc <= nAdcs; ++adc) {
        //      std::cout << uplink << "\t" << adc << "\n";
        if (gains.at(uplink).at(adc))
          correctedAdcVals[uplink - uplinkMin][adc - 1] =
              (rawAdcVals[uplink - uplinkMin][adc - 1] -
               pedestals.at(uplink).at(adc)) /
              gains.at(uplink).at(adc);  // if gain is not 0
        else
          correctedAdcVals[uplink - uplinkMin][adc - 1] =
              0.0;  // gain fit failed!
        //      std::cout << correctedAdcVals[uplink - uplinkMin][adc-1] <<
        //      std::endl;
      }
    }
    correctedTree->Fill();
  }
  correctedTree->ResetBranchAddresses();
  correctedTree->Write();
  correctedFile.Close();
}
