//from std
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
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

//from BOOST
#include "boost/program_options.hpp"

namespace po = boost::program_options;

struct config{
  unsigned int nGaussians;
  unsigned int uplinkNumber;
  std::string inputFile;
  bool debug;
};

void GetGain(TTree* t, const unsigned int uplinkMin, const unsigned int uplinkMax, const unsigned int adcIDmin, const unsigned int adcIDmax, std::string saveName = {""}, const unsigned int maxGaussians = {4}){
  std::ofstream gainFile(saveName);
  gainFile << "uplinkNumber\tadcNumber\tgain\tgainErr" << std::endl;

  TString canvasName = saveName;
  canvasName.ReplaceAll(".root", ".pdf");

  TCanvas can("", "");
  can.SaveAs( (canvasName + "[").Data() );

  double numFits{0.};
  double failedFits{0.};

  TString branchNameTemplate = "Uplink_";
  for(unsigned int uplinkNumber = uplinkMin; uplinkNumber<=uplinkMax; ++uplinkNumber){
    TString branchNameTemplate2 = branchNameTemplate + std::to_string(uplinkNumber) + "_";

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
      Float_t* peaks = spec.GetPositionX();

      Float_t max_peak = *std::max_element(peaks, peaks+nGaussians);
      Float_t min_peak = *std::min_element(peaks, peaks+nGaussians);

//      std::cout << "Found peaks at ";
//      for(unsigned int k = 0; k<nGaussians; ++k){
//        std::cout << peaks [k] << "\t";
//      }
//      std::cout << std::endl;

      RooRealVar adcCount(branchName, "adcCount", min_peak-30., max_peak+30.);
      RooRealVar mean("mean", "", min_peak, 420., 650.);
      RooRealVar sigma("sigma", "", 10., 0., 20.);
      RooRealVar gain("gain", "", (max_peak-min_peak)/(nGaussians-1.), 45., 75.);
      RooDataSet dataSet("dataSet", "", t, RooArgSet(adcCount));

      std::vector<RooGaussian*> gaussians(nGaussians, nullptr);
      std::vector<RooFormulaVar*> means(nGaussians, nullptr);
      std::vector<RooRealVar*> fractions(nGaussians-1, nullptr);

      RooArgList gaussianList("gaussianList");
      RooArgList fractionList("fractionList");
      for(unsigned int j=0; j<nGaussians; ++j){
        means[j] = new RooFormulaVar( TString("mean_" + std::to_string(j) ).Data(), "", TString("mean+gain*" + std::to_string(j) ).Data(), RooArgList(mean, gain));
        gaussians[j] = new RooGaussian( TString("gaussian_" + std::to_string(j) ).Data(), "", adcCount, *means[j], sigma);
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
      fs->Print("v");

      if(!isAGoodFit){
        ++failedFits;
        std::cout << "Fit for uplink " << uplinkNumber << " and adc " << adcNumber << " failed!" << std::endl;
      }

      gainFile << uplinkNumber << "\t" << adcNumber << "\t" << gain.getVal() << "\t" << gain.getError() << std::endl;

      RooPlot* myFrame = adcCount.frame();
      dataSet.plotOn(myFrame);
      pdf.plotOn(myFrame);
      myFrame->Draw();
      can.SaveAs( (canvasName).Data() );

      for(unsigned int i=0; i<nGaussians; ++i){
        delete means[i];
        delete gaussians[i];
        if(i>0){
          delete fractions[i-1];
        }
      }
    }
  }
  can.SaveAs( (canvasName + "]").Data() );

  std::cout << failedFits/numFits*100. << "% of all fits failed." << std::endl;
}

void searchBranches(TTree* inputTree, unsigned int& uplinkMin, unsigned int& uplinkMax, unsigned int& adcIDmin, unsigned int& adcIDmax){
  TString branchNameTemplate = "Uplink_";
  for(unsigned int uplinkNumber = uplinkMin; uplinkNumber<=uplinkMax; ++uplinkNumber){
    TString branchNameTemplate2 = branchNameTemplate + std::to_string(uplinkNumber) + "_";

    for(unsigned int adcNumber = adcIDmin; adcNumber<=adcIDmax; ++adcNumber){
      unsigned int nGaussians = maxGaussians;
      TString branchName = branchNameTemplate2 + std::to_string(adcNumber);
}

int parseOptions(config &c, int argc, char *argv[]){

  // declare options
  po::options_description desc("Allowed options") ;
  desc.add_options()
    ("help", "show this help")
    ("ngaussians,n", po::value<unsigned int>(&c.nGaussians)->default_value(4), "number aus gaussians to fit to the distribution")
    ("inputfile,f", po::value<std::string>(&c.inputFile), "test beam data file containing the led run")
    ("uplink,u", po::value<unsigned int>(&c.uplinkNumber), "number of the uplink")
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

  lhcb::lhcbStyle();
  if(!c.debug){
    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
    RooMsgService::instance().setSilentMode(true);
  }
  TFile inputFile(c.inputFile.c_str(), "READ");
  TTree* inputTree = dynamic_cast<TTree*>(inputFile.Get("rawData"));

  std::vector<unsigned int> uplinks;
  std::vector<unsigned int> adcIDs;
  unsigned int uplinkMin;
  unsigned int uplinkMax;
  unsigned int adcIDmin;
  unsigned int adcIDmax;
  searchBranches(inputTree, uplinkMin, uplinkMax, adcIDmin, adcIDmax);

  GetGain(inputTree, uplinks, adcIDs, "gainsFile");

  return 0;
}


