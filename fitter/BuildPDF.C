// Author: Batoul Diab
#include <vector>
#include <stdio.h>
#include <iostream>
#include <algorithm>
#include "TF1.h"

#include "RooCategory.h"
#include "RooNumIntConfig.h"
#include "RooPlotable.h"
#include "RooRealVar.h"
#include "RooAbsReal.h"
#include "RooGaussian.h"
#include "RooUniform.h"
#include "RooPlot.h"
#include "RooDataHist.h"
#include "RooHist.h"
#include "RooDataSet.h"
#include "RooArgList.h"
#include "RooAddPdf.h"
#include "RooPoisson.h"
#include "RooExponential.h"
#include "RooFormulaVar.h"
#include "RooBifurGauss.h"
#include "RooExtendPdf.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"
#include "RooPolynomial.h"
#include "RooFormulaVar.h"

#include "SetRangesAndLabels.C"

using namespace  RooFit;

void BuildPDF(RooWorkspace* ws, map<string, string> parIni, TH1* dist) {

  for(auto it = parIni.cbegin(); it != parIni.cend(); ++it) {
    if(it->first.find("fit")!=std::string::npos) continue;
    
    RooRealVar* var (0x0);
    // read the initial parameter string and get the initial values
    cout <<"[INFO] Let's treat variable "<<it->first<<" with values "<<it->second<<endl;
    
    char delimiter = ',';
    std::string processedInput = removeSurroundingBrackets(it->second);
    std::vector<std::string> substrings = splitString(processedInput, delimiter);
    if (substrings.size() == 3)
      var = new RooRealVar(it->first.c_str(),it->first.c_str(), std::stof(substrings[0]), std::stof(substrings[1]), std::stof(substrings[2]));
    else if (substrings.size() == 1)
      var = new RooRealVar(it->first.c_str(),it->first.c_str(), std::stof(substrings[0]));
    else { cout <<"[ERROR] This initial parameter does not have the correct number of values."<<endl; continue;}

    ws->import(*var);
  }

  RooRealVar* phi_ws = ws->var("phi"); //just to make it easier to read
  
  RooRealVar* f_ns = new RooRealVar("f_ns","f_ns", 0.2*dist->GetEntries(), 0.1, dist->GetEntries());
  RooRealVar* f_fs = new RooRealVar("f_fs","f_fs", 0.1*dist->GetEntries(), 0.1, dist->GetEntries());
  RooRealVar* f_b = new RooRealVar("f_b","f_b", dist->GetMinimum()*dist->GetNbinsX(), dist->GetMinimum()*dist->GetNbinsX()-3, dist->GetMinimum()*dist->GetNbinsX()+3);

  phi_ws->setRange("nearWindow", -0.5*TMath::Pi(), 0.5*TMath::Pi());
  phi_ws->setRange("nearPeakWindow", -0.2*TMath::Pi(), 0.2*TMath::Pi());
  phi_ws->setRange("farWindow", 0.5*TMath::Pi(), 1.5*TMath::Pi());

  // Create extended PDFs
  RooExtendPdf* nearSidePDF(0x0);
  RooExtendPdf* farSidePDF(0x0);
  RooExtendPdf* bkgPDF(0x0);
  
  if (parIni["fitBkg"]=="Uniform") {
    //RooPolynomial* bkgUni = new RooPolynomial("bkgUni", "Uniform background", *phi_ws, RooArgList(*ws->var("norm_bkg"))); //ws->import(*bkgUni);
    RooUniform* bkgUni = new RooUniform("bkgUni", "Uniform background", RooArgSet(*phi_ws, *ws->var("norm_bkg")));
    bkgPDF = new RooExtendPdf("bkgPDF", "bkgPDF", *bkgUni, *f_b);
  }

  if (parIni["fitSig"]=="Gauss") { 
    RooFormulaVar* sigma_ns2 = new RooFormulaVar("sigma_ns2", "sigmaRatio_ns * sigma_ns", RooArgList(*ws->var("sigmaRatio_ns"), *ws->var("sigma_ns"))); 
    RooGaussian* nearGaus1 = new RooGaussian("nearGaus1", "near side Gaussian 1", *phi_ws, *ws->var("mean_ns"), *ws->var("sigma_ns")); 
    RooGaussian* nearGaus2 = new RooGaussian("nearGaus2", "near side Gaussian 2", *phi_ws, *ws->var("mean_ns"), *sigma_ns2);
    // Define the ranges for nearGaus1 and nearGaus2

    // Create extended PDFs for nearGaus1 and nearGaus2 with specified ranges
    RooExtendPdf* nearGaus1Ext = new RooExtendPdf("nearGaus1Ext", "near side Gaussian 1 extended", *nearGaus1, *ws->var("f_ns_1"), "nearWindow");
    RooExtendPdf* nearGaus2Ext = new RooExtendPdf("nearGaus2Ext", "near side Gaussian 2 extended", *nearGaus2, *ws->var("f_ns_1"), "nearPeakWindow");

    
    RooAddPdf* nearGaus = new RooAddPdf("nearGaus", "fit of near side", RooArgList(*nearGaus1Ext, *nearGaus2Ext), *ws->var("f_ns_1")); 
    RooGaussian* farGaus = new RooGaussian("farGaus", "far side Gaussian", *phi_ws, *ws->var("mean_fs"), *ws->var("sigma_fs")); 
    nearSidePDF = new RooExtendPdf("nearSidePDF", "nearSidePDF", *nearGaus, *f_ns, "nearWindow");
    farSidePDF = new RooExtendPdf("farSidePDF", "farSidePDF", *farGaus, *f_fs, "farWindow");
  }
    
  // Combine the components into a composite model
  RooAddPdf* model =  new RooAddPdf ("model", "model", RooArgList(*nearSidePDF, *farSidePDF, *bkgPDF));
   ws->import(*model);

}
