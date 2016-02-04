//from std
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>

//from ROOT
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TSpectrum.h"
#include "TH1D.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooMinuit.h"
#include "RooAbsReal.h"

//from here
#include "lhcbStyle.h"
#include "Cluster.h"
#include "Calibration.h"



TString removePath(TString str){
  Ssiz_t posOfSlash = str.Last('/');
  if(posOfSlash != kNPOS){
    str = str.Replace(0, posOfSlash+1, "");
  }
  return str;
}

std::vector<Event*>* parseRootTree(TTree* dataTree,
                                    unsigned int uplinkMin,
                                    unsigned int uplinkMax,
                                    unsigned int nAdcs,
                                    const std::map<unsigned int, std::map<unsigned int, double>>& pedestals,
                                    const std::map<unsigned int, std::map<unsigned int, double>>& gains){

  std::vector<Event*>* dataVector = new std::vector<Event*>(dataTree->GetEntriesFast());

  const unsigned int nUplinks = uplinkMax-uplinkMin+1;
  float** adcVals = new float*[nUplinks];
  for(unsigned int i = 0; i< nAdcs; ++i){
    adcVals[i] = new float[nAdcs];
  }

  dataTree->SetBranchStatus("*", 0);
  for(unsigned int uplink=uplinkMin; uplink<=uplinkMax; ++uplink){
    std::string branchName = "Uplink_" + std::to_string(uplink) + "_adc_";
    for(unsigned int adc = 0; adc<nAdcs; ++adc){
      dataTree->SetBranchStatus((branchName + std::to_string(adc+1)).c_str(), 1);
      dataTree->SetBranchAddress((branchName + std::to_string(adc+1)).c_str(), &adcVals[uplink-uplinkMin][adc]);
    }
  }
    for(unsigned int i = 0; i < dataTree->GetEntriesFast(); ++i){
      dataTree->GetEntry(i);
      Event* event = new Event(nAdcs * (uplinkMax-uplinkMin+1) );
      for(unsigned int uplink=uplinkMin; uplink<=uplinkMax; ++uplink){
        for(unsigned int adc = 0; adc < nAdcs; ++adc){
          Channel c;
          c.Uplink = uplink;
          c.ChannelNumber = adc+1;
          c.AdcValue = (adcVals[uplink-uplinkMin][adc] - pedestals.at(uplink).at(adc+1)) / gains.at(uplink).at(adc+1);
          event->at(adc + (uplink - uplinkMin)*nAdcs ) = c;
        }
      }
      dataVector->at(i) = event;
    }

  dataTree->ResetBranchAddresses();
  for(unsigned int i = 0; i< nUplinks; ++i){
    delete[] adcVals[i];
  }
  delete[] adcVals;

  return dataVector;
}

std::vector<Event*>* parseCorrectedRootTree(TTree* dataTree,
                                            unsigned int uplinkMin,
                                            unsigned int uplinkMax,
                                            unsigned int nAdcs,
                                            double factor
                                          ){

  std::vector<Event*>* dataVector = new std::vector<Event*>(dataTree->GetEntriesFast());

  if(factor != 1.) std::cout << "Parsing corrected root tree and correct the adc value for the factor of " << factor << "!\n";

  const unsigned int nUplinks = uplinkMax-uplinkMin+1;
  float** adcVals = new float*[nUplinks];
  for(unsigned int i = 0; i< nUplinks; ++i){
    adcVals[i] = new float[nAdcs];
  }

//  std::vector<std::vector<float> > adcVals(nUplinks, std::vector<float>(nAdcs));

  dataTree->SetBranchStatus("*", 0);
  for(unsigned int uplink=uplinkMin; uplink<=uplinkMax; ++uplink){
    std::string branchName = "Uplink_" + std::to_string(uplink) + "_adc_";
    for(unsigned int adc = 0; adc<nAdcs; ++adc){
      dataTree->SetBranchStatus((branchName + std::to_string(adc+1)).c_str(), 1);
      dataTree->SetBranchAddress((branchName + std::to_string(adc+1)).c_str(), &adcVals[uplink-uplinkMin][adc]);
    }
  }
  std::cout << "channels per event: " << ( nAdcs * (uplinkMax-uplinkMin+1) ) << "\n";
    for(unsigned int i = 0; i < dataTree->GetEntriesFast(); ++i){
      dataTree->GetEntry(i);
      Event* event = new Event(nAdcs * (uplinkMax-uplinkMin+1) );
      for(unsigned int uplink=uplinkMin; uplink<=uplinkMax; ++uplink){
        for(unsigned int adc = 0; adc < nAdcs; ++adc){
//          std::cout << "stroring adcVals[" << uplink-uplinkMin << "][" << adc << "] at " <<  (adc + (uplink - uplinkMin)*nAdcs) << "\n";
          Channel c;
          c.Uplink = uplink;
          c.ChannelNumber = adc+1;
          if(adcVals[uplink-uplinkMin][adc]>0){
            c.AdcValue = adcVals[uplink-uplinkMin][adc] * factor;
          }
          else
            c.AdcValue = 0;
          //if(adcVals[uplink-uplinkMin][adc]>10 || adcVals[uplink-uplinkMin][adc] <-10){
          //}
//          std::cout << uplink << "\t" << adc << "\t" << adcVals[uplink-uplinkMin][adc] << "\n";
          event->at(adc + (uplink - uplinkMin)*nAdcs ) = c;
        }
      }
      dataVector->at(i) = event;
    }

//  dataTree->ResetBranchAddresses();
  for(unsigned int i = 0; i< nUplinks; ++i){
    delete[] adcVals[i];
  }
  delete[] adcVals;

  return dataVector;
}



void produceGains(TTree* t,
                  const unsigned int uplinkMin,
                  const unsigned int uplinkMax,
                  const unsigned int adcIDmin,
                  const unsigned int adcIDmax,
                  const unsigned int maxGaussians,
                  TString fileName,
                  std::string savePath){
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
  RooMsgService::instance().setSilentMode(true);
  lhcb::lhcbStyle();
  TString saveName = fileName;
  saveName.Remove(0, saveName.Last('/')+1);
  saveName.ReplaceAll(".root", "_gain.txt");
  std::ofstream gainFile(savePath + saveName);
  gainFile << "uplinkNumber\tadcNumber\tgain\tgainErr" << std::endl;

  TString canvasName = saveName;
  canvasName.ReplaceAll(".txt", ".pdf");

  TCanvas can("", "");
  can.SaveAs( (savePath + canvasName + "[").Data() );

  double numFits{0.};
  double failedFits{0.};

  TString branchNameTemplate = "Uplink_";
  for(unsigned int uplinkNumber = uplinkMin; uplinkNumber<=uplinkMax; ++uplinkNumber){
    TString branchNameTemplate2 = branchNameTemplate + std::to_string(uplinkNumber) + "_adc_";

    for(unsigned int adcNumber = adcIDmin; adcNumber<=adcIDmax; ++adcNumber){
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

      auto max_peak = *std::max_element(peaks, peaks+nGaussians);
      auto min_peak = *std::min_element(peaks, peaks+nGaussians);

//      std::cout << "Found peaks at ";
//      for(unsigned int k = 0; k<nGaussians; ++k){
//        std::cout << peaks [k] << "\t";
//      }
//      std::cout << std::endl;

      RooRealVar adcCount(branchName, branchName, min_peak-30., max_peak+30.);
      RooRealVar mean("mean", "", min_peak, 420., 650.);
      RooRealVar sigma1("sigma1", "", 7., 0., 20.);
      RooRealVar sigma2("sigma2", "", 4., 0., 20.);
      RooRealVar gain("gain", "", (max_peak-min_peak)/(nGaussians-1.), 45., 75.);
      RooDataSet dataSet("dataSet", "", t, RooArgSet(adcCount));

      std::vector<RooGaussian*> gaussians(nGaussians, nullptr);
      std::vector<RooFormulaVar*> means(nGaussians, nullptr);
      std::vector<RooFormulaVar*> sigmas(nGaussians, nullptr);
      std::vector<RooRealVar*> fractions(nGaussians-1, nullptr);

      RooArgList gaussianList("gaussianList");
      RooArgList fractionList("fractionList");
      for(unsigned int j=0; j<nGaussians; ++j){
        means[j] = new RooFormulaVar( TString("mean_" + std::to_string(j) ).Data(), "", TString("mean+gain*" + std::to_string(j) ).Data(), RooArgList(mean, gain));
        sigmas[j] = new RooFormulaVar( TString("sigma_" + std::to_string(j) ).Data(), "", TString("sqrt(sigma1*sigma1+" + std::to_string(j) + "*sigma2*sigma2)").Data(), RooArgList(sigma1, sigma2));
        gaussians[j] = new RooGaussian( TString("gaussian_" + std::to_string(j) ).Data(), "", adcCount, *means[j], *sigmas[j]);
        gaussianList.add(*gaussians[j]);
        if(j>0){
          fractions[j-1] = new RooRealVar( TString("fraction_" + std::to_string(j) ).Data(), "", 1./2./j, 0., 1.);
          fractionList.add(*fractions[j-1]);
        }
      }
      RooAddPdf pdf("pdf", "", gaussianList, fractionList);

      //RooFitResult* fs = pdf.fitTo(dataSet, RooFit::Save(true));
      //fs->Print("v");
      ++numFits;
      RooAbsReal* nll = pdf.createNLL(dataSet);

      RooMinuit minu(*nll);
      minu.setStrategy(2); // ensure minimum true, errors correct

      // call MIGRAD -- minimises the likelihood
      int migradStatusCode = minu.migrad(); // catch status code

      // could also get an intermediary RooFitResult if you want
    //      RooFitResult *migradFitResult = m.save();
    //      int migradStatusCode2 = migradFitResult->status();

      // call HESSE -- calculates error matrix
      int hesseStatusCode = minu.hesse();

      RooFitResult *fs = minu.save();  // equivalent result to that from fitTo call


      // check for success
      int covarianceQuality = fs->covQual();

  //    int minosStatusCode = minu.minos();


      bool isAGoodFit = ((hesseStatusCode==0) && (migradStatusCode==0) && (covarianceQuality==3));// && (minosStatusCode==0) );


      if(!isAGoodFit || gain.getError() == 0.0){
        ++failedFits;
        std::cout << "Fit for uplink " << uplinkNumber << " and adc " << adcNumber << " failed!" << std::endl;
        fs->Print("v");
        gainFile << uplinkNumber << "\t" << adcNumber << "\t" << 0 << "\t" << 0 << std::endl;
      }
      else gainFile << uplinkNumber << "\t" << adcNumber << "\t" << gain.getVal() << "\t" << gain.getError() << std::endl;

//      fs->Print("v");

      RooPlot* myFrame = adcCount.frame();
      dataSet.plotOn(myFrame);
      pdf.plotOn(myFrame);
      myFrame->Draw();
      can.SaveAs( (savePath + canvasName).Data());

      for(unsigned int i=0; i<nGaussians; ++i){
        delete means[i];
        delete gaussians[i];
        if(i>0){
          delete fractions[i-1];
        }
      }
    }
  }
  can.SaveAs( (savePath + canvasName + "]").Data() );

  std::cout << failedFits/numFits*100. << "% of all fits failed." << std::endl;
}

std::map<unsigned int, std::map<unsigned int, double>> readGains(std::string fileName){
  std::ifstream inputFile(fileName);
  std::string line;
  int uplinkNumber{0};
  int adcNumber{0};
  double gain{0};

//  std::map<unsigned int, double> means;
//  std::map<unsigned int, double> goodGains;


  std::map<unsigned int, std::map<unsigned int, double>> gains;
  while (std::getline(inputFile, line)) {
      std::istringstream ss(line);
      if(ss >> uplinkNumber >> adcNumber >> gain){
        std::cout << uplinkNumber << "\t" << adcNumber << "\t" << gain << std::endl;

//        if(means.find( uplinkNumber ) == means.end()){
//          means[uplinkNumber] = 0.;
//          goodGains[uplinkNumber] = 0.;
//        }

        gains[uplinkNumber][adcNumber] = gain;

//        if(!gain == 0){
//          means[uplinkNumber] += gain;
//          ++goodGains[uplinkNumber];
//        }
      }
  }
  // replace gains from failed fits by the mean of all gains of the same uplink maybe not a good idea!
//  for(auto& uplink : gains){
//    means[uplink.first] /= goodGains[uplink.first];
//    for(auto& adc : uplink.second){
//      if (adc.second == 0){
//        std::cout << "replacing " << adc.second << " by " << means[uplink.first] << std::endl;
//        adc.second = means[uplink.first];
//      }
//    }
//  }
  return gains;
}

calibrationRunNumbers lookUpCalibrationFiles(const unsigned int runNumber, const std::string catalogueFileName){
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
      if(currentRunNumber == runNumber){
        rv.dark = currentDarkNumber;
        rv.led = currentLedNumber;
        break;
      }
    }
  return rv;
  }

unsigned int runNumberFromFilename(std::string filename){
  if(size_t posOfSlash = filename.rfind("/") != std::string::npos){
    filename = filename.substr(posOfSlash+1);
  }
  return std::stoi(filename.substr( filename.find("_")+1, filename.find("_", filename.find("_")+1)-filename.find("_")-1 ));
}

std::map<unsigned int, std::map<unsigned int, double> > getPedestals(const std::string fileName,
                                                                     const unsigned int uplinkMin,
                                                                     const unsigned int uplinkMax,
                                                                     const unsigned int nAdcs){
  TFile pedestalFile(fileName.c_str(), "READ");
  if(!pedestalFile.IsOpen()){
    std::cerr << "unable to open pedestal file " << fileName << std::endl;
  }
  TTree* pedestalTree = dynamic_cast<TTree*>(pedestalFile.Get("rawData"));
  if(pedestalTree == nullptr){
    std::cerr << "pedestal tree is nullptr" << std::endl;
  }
  const unsigned int nUplinks = uplinkMax-uplinkMin+1;
  std::cout << "creating double array of form double[" << nUplinks << "][" << nAdcs << "]\n";
  float** rawAdcVals = new float*[nUplinks];
  for(unsigned int i = 0; i< nUplinks; ++i){
    rawAdcVals[i] = new float[nAdcs];
  }
  std::cout << "setting branch addresses" << std::endl;
  std::string branchNameTemplate1 = "Uplink_";
  for(unsigned int uplink = uplinkMin; uplink<=uplinkMax; ++uplink){
    std::string branchNameTemplate2 = branchNameTemplate1 + std::to_string(uplink) + "_adc_";
    for(unsigned int adc = 1; adc<=nAdcs; ++adc){
      std::string branchName = branchNameTemplate2 + std::to_string(adc);
      std::cout << "setting adress for branch " << branchName << " to double[" << uplink - uplinkMin << "][" << adc-1 << "]"  <<  std::endl;
      pedestalTree->SetBranchAddress(branchName.c_str(), &rawAdcVals[uplink - uplinkMin][adc-1]);
    }
  }
  std::cout << "loop over file" << std::endl;
  std::map<unsigned int, std::map<unsigned int, double> > pedestals;
  for(unsigned int i =0; i<pedestalTree->GetEntriesFast(); ++i){
    pedestalTree->GetEntry(i);
    for(unsigned int uplink = uplinkMin; uplink<=uplinkMax; ++uplink){
      for(unsigned int adc = 1; adc <= nAdcs; ++adc){
        pedestals[uplink][adc] += rawAdcVals[uplink - uplinkMin][adc-1];
      }
    }
  }

  for(unsigned int uplink = uplinkMin; uplink<=uplinkMax; ++uplink){
    for(unsigned int adcNum = 1; adcNum <= nAdcs; ++adcNum){
      pedestals[uplink][adcNum] /= pedestalTree->GetEntriesFast();
    }
  }
  pedestalTree->ResetBranchAddresses();
  pedestalFile.Close();
  for(unsigned int i = 0; i< nUplinks; ++i){
    delete[] rawAdcVals[i];
  }
  delete[] rawAdcVals;

  return pedestals;
}


void correctFile(TTree* tree2correct,
                 const std::map<unsigned int, std::map<unsigned int, double> > &gains,
                 const std::map<unsigned int, std::map<unsigned int, double> > &pedestals,
                 const unsigned int uplinkMin,
                 const unsigned int uplinkMax,
                 const unsigned int nAdcs,
                 TString newFileName){

//  newFileName = removePath(newFileName);
//  TFile correctedFile(("/data/testbeam/data/corrected/" + newFileName).Data(), "RECREATE");

  TFile correctedFile(newFileName.Data(), "RECREATE");
  TTree* correctedTree = new TTree("rawData", "rawData");

  //setup fill mechanism for new TTree
  const unsigned int nUplinks = uplinkMax-uplinkMin+1;
  float** correctedAdcVals = new float*[nUplinks];
  for(unsigned int i = 0; i< nUplinks; ++i){
    correctedAdcVals[i] = new float[nAdcs];
  }
  float** rawAdcVals = new float*[nUplinks];
  for(unsigned int i = 0; i< nUplinks; ++i){
    rawAdcVals[i] = new float[nAdcs];
  }

  std::string branchNameTemplate1 = "Uplink_";
  for(unsigned int uplink = uplinkMin; uplink<=uplinkMax; ++uplink){
    std::string branchNameTemplate2 = branchNameTemplate1 + std::to_string(uplink) + "_adc_";
    for(unsigned int adc = 0; adc<nAdcs; ++adc){
      std::string branchName = branchNameTemplate2 + std::to_string(adc+1);
      tree2correct->SetBranchAddress(branchName.c_str(), &rawAdcVals[uplink - uplinkMin][adc]);
      correctedTree->Branch(branchName.c_str(), &correctedAdcVals[uplink - uplinkMin][adc], (branchName + "/F").c_str() );
    }
  }
  for(unsigned int i=0; i<tree2correct->GetEntriesFast(); ++i){
    tree2correct->GetEntry(i);
    for(unsigned int uplink = uplinkMin; uplink<=uplinkMax; ++uplink){
      for(unsigned int adc = 1; adc<=nAdcs; ++adc){
//      std::cout << uplink << "\t" << adc << "\n";
      if (gains.at(uplink).at(adc)) correctedAdcVals[uplink - uplinkMin][adc-1] = ( rawAdcVals[uplink - uplinkMin][adc-1] - pedestals.at(uplink).at(adc) ) / gains.at(uplink).at(adc); //if gain is not 0
      else correctedAdcVals[uplink - uplinkMin][adc-1] = 0.0; // gain fit failed!
//      std::cout << correctedAdcVals[uplink - uplinkMin][adc-1] << std::endl;
      }
    }
    correctedTree->Fill();
  }
  correctedTree->ResetBranchAddresses();
  correctedTree->Write();
  correctedFile.Close();

  for(unsigned int i = 0; i< nUplinks; ++i){
    delete[] correctedAdcVals[i];
    delete[] rawAdcVals[i];
  }
  delete[] correctedAdcVals;
  delete[] rawAdcVals;
}
