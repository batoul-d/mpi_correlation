// Author: Batoul Diab
#include <vector>
#include <iostream>
#include <algorithm>
#include <stdio.h>
#include <TStyle.h>
#include "TFile.h"
#include "TDirectoryFile.h"
#include "TList.h"
#include "THnSparse.h"
#include "TH1.h"
#include "TLine.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TString.h"
#include "TMath.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMatrixD.h"
#include "TFitResult.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include <fstream>
#include <string>
#include <sstream>

#include "BuildPDF.C"


TH1 *GetCorrHisto(const char *filePattern, string rangeLabel);
//___________________________________________________________________________________________________________
//___________________________________________________________________________________________________________
map<string, double> SignalExtraction1Fit(const char *file, int minMult, int maxMult, float minPTtrig, float maxPTtrig, float minPTassoc, float maxPTassoc, map<string, string> parIni, RooWorkspace* ws, string rangeLabel) {
  
  map<string, double> results;
  TH1 *dPhiDist = GetCorrHisto(file, rangeLabel);
  
  RooRealVar* phi = new RooRealVar("phi","#Delta#phi", -0.5*TMath::Pi(), 1.5*TMath::Pi());
  RooDataHist* phiDist = new RooDataHist("phiDist", "#Delta#phi distribution", *phi, dPhiDist);

  ws->import(*phi); ws->import(*phiDist);

  // bulding the model from input
  BuildPDF(ws, parIni, dPhiDist);

  // the actual fit
  RooFitResult* fitResult = ws->pdf("model")->fitTo(*phiDist, Extended(kTRUE), SumW2Error(true), RooFit::Save());
  
  std::string canvTitle = Form("dphi_%s_%s_%s", parIni["fitBkg"].c_str(), parIni["fitSig"].c_str(), rangeLabel.c_str());
  
  TCanvas *cPhiDist = new TCanvas("cPhiDist", "c", 800, 800);
  TPad *padHist = new TPad("padHist","",0,.23,1,1);
      padHist->SetBottomMargin(0.015);
      TPad *padPull = new TPad("padPull","",0,0,1,.228);
      padPull->SetTopMargin(0.016);
      padPull->SetBottomMargin(0.3);
      padHist->cd();

      
      RooPlot* phiframe = ws->var("phi")->frame();
      phiDist->plotOn(phiframe, DataError(RooAbsData::SumW2));
      ws->pdf("model")->plotOn(phiframe, Name("bkg"), Components(*ws->pdf("bkgPDF")), DrawOption("F"), LineColor(kBlack), FillColor(kGray));
      
      phiDist->plotOn(phiframe);
      ws->pdf("model")->plotOn(phiframe, Name("model"), LineColor(kRed));

      //the pull and chi2 need to be calculated right after plotting the model
      RooHist *hpull = phiframe->pullHist();
      RooPlot* pullframe = ws->var("phi")->frame(Title("Pull Distribution"));
      pullframe->addPlotable(hpull,"P");

      int nPar = ws->pdf("model")->getParameters(*phiDist)->selectByAttrib("Constant",kFALSE)->getSize();
      double chi2 = phiframe->chiSquare(nPar);

      ws->pdf("model")->plotOn(phiframe, Name("nearSide"), Components(*ws->pdf("nearSidePDF")), LineStyle(kDashed), LineColor(kBlue));
      ws->pdf("model")->plotOn(phiframe, Name("farSide"), Components(*ws->pdf("farSidePDF")), LineStyle(kDashed), LineColor(kMagenta+2));
      
      //cout <<"[INFO] done with plotting the pdfs"<<endl;

      
      phiframe->SetTitle(Form("%g < p_{T,trig} < %g, %g < p_{T, assoc} < %g, mult. %d-%d", minPTtrig, maxPTtrig, minPTassoc, maxPTassoc, minMult, maxMult));
      phiframe->GetYaxis()->SetLabelSize(0.04);
      phiframe->GetYaxis()->SetTitleSize(0.04);
      phiframe->GetYaxis()->SetTitleOffset(1.);
      phiframe->GetYaxis()->SetTitleFont(42);
      phiframe->GetYaxis()->CenterTitle(kTRUE);

      phiframe->GetXaxis()->SetLabelSize(0);
      phiframe->GetXaxis()->SetTitle("");

      phiframe->Draw();

      TLatex *textAlice = new TLatex();
      textAlice->SetNDC();
      textAlice->SetTextAlign(12);
      textAlice->SetTextFont(43);
      textAlice->SetTextSize(17); // Size in pixel height
      textAlice->DrawLatex(0.5, 0.86, "ALICE pp @ #sqrt{s} = 13.6 TeV");

      TLatex *textVar = new TLatex();
      textVar->SetNDC();
      textVar->SetTextAlign(12);
      textVar->SetTextFont(43);
      textVar->SetTextSize(17); // Size in pixel height
      float yText = 0.55;
      textVar->DrawLatex(0.4, yText, Form("%s = %.3f #pm  %.3f", varFancyLabel("f_ns").c_str(), ws->var("f_ns")->getValV(), ws->var("f_ns")->getError())); yText = yText-0.03;
      textVar->DrawLatex(0.4, yText, Form("%s = %.3f #pm  %.3f", varFancyLabel("f_fs").c_str(), ws->var("f_fs")->getValV(), ws->var("f_fs")->getError())); yText = yText-0.03;
      textVar->DrawLatex(0.4, yText, Form("%s = %.3f #pm  %.3f", varFancyLabel("f_b").c_str(), ws->var("f_b")->getValV(), ws->var("f_b")->getError())); yText = yText-0.03;
      
      for (auto it = parIni.cbegin(); it != parIni.cend(); ++it) {
	if(it->first.find("fit")!=std::string::npos) continue;
	if(it->first.find("mean")!=std::string::npos) continue;
	textVar->DrawLatex(0.4, yText, Form("%s = %.3f #pm  %.3f", varFancyLabel(it->first.c_str()).c_str(), ws->var(it->first.c_str())->getValV(), ws->var(it->first.c_str())->getError()));
	yText = yText-0.03;
      }
      textVar->DrawLatex(0.4, yText, Form("#chi^{2}/ndof =  %.3f",chi2 ));
      /////////////////////////////////
      
      TLegend *legend_lines = new TLegend(0.15, 0.4, 0.3, 0.55); // Define the legend position and size
      legend_lines->SetBorderSize(0);
      legend_lines->SetFillStyle(0);
      legend_lines->SetTextAlign(12);
      legend_lines->SetTextFont(43);
      legend_lines->SetTextSize(17);
      legend_lines->AddEntry(phiframe->RooPlot::findObject("phiDist"), "Data", "EP");
      legend_lines->AddEntry(phiframe->RooPlot::findObject("model"), "Total fit", "L");
      legend_lines->AddEntry(phiframe->RooPlot::findObject("nearSide"), "Near side", "L");
      legend_lines->AddEntry(phiframe->RooPlot::findObject("farSide"), "Far side", "L");
      legend_lines->AddEntry(phiframe->RooPlot::findObject("bkg"), "Bakground", "f");
      legend_lines->Draw("same");
      
      // make the histograms look better
      padPull->cd();
      pullframe->SetTitle("");
      pullframe->GetYaxis()->SetTitle("Pull");
      pullframe->GetYaxis()->SetRangeUser(-5,5);
      pullframe->GetYaxis()->SetLabelSize(0.07);
      pullframe->GetYaxis()->SetTitleSize(0.1);
      pullframe->GetYaxis()->SetTitleOffset(0.2);
      pullframe->GetYaxis()->SetTitleFont(42);
      pullframe->GetYaxis()->CenterTitle(kTRUE);
      
      pullframe->GetXaxis()->SetLabelSize(0.1);
      pullframe->GetXaxis()->CenterTitle(kTRUE);

      pullframe->GetXaxis()->SetTitleSize(0.15);
      pullframe->GetXaxis()->SetTitleFont(42);
      pullframe->GetXaxis()->SetTitleOffset(0.7);
     
      pullframe->Draw();
      //hpull->Draw();
      TLine* pL = new TLine(-0.5*TMath::Pi(),0, 1.5*TMath::Pi(),0);
      pL->SetLineColor(kRed);
      pL->SetLineStyle(kDashed);
      pL->SetLineWidth(2);
      pL->Draw();
      
      cPhiDist->cd();
      padHist->Draw();
      padPull->Draw();

      std::string pdfPath = "";
      cPhiDist->SaveAs(("output/" + pdfPath + canvTitle + ".pdf").c_str());
 
  for (auto it = parIni.cbegin(); it != parIni.cend(); ++it) {
    if(it->first.find("fit")!=std::string::npos) continue;
    results[it->first.c_str()] = ws->var(it->first.c_str())->getValV();
    results[Form("%s_err",it->first.c_str())] = ws->var(it->first.c_str())->getError();
  }
  
  results["f_ns"] = ws->var("f_ns")->getValV();
  results["f_ns_err"] = ws->var("f_ns")->getError();
  results["f_fs"] = ws->var("f_fs")->getValV();
  results["f_fs_err"] = ws->var("f_fs")->getError();
  results["f_b"] = ws->var("f_b")->getValV();
  results["f_b_err"] = ws->var("f_b")->getError();
    
  results["chi2ndf"] = chi2;
  results["ndf"] = nPar;//fitResult->ndf();

  return results;
  
} //end of SignalExtraction Function

//___________________________________________________________________________________________________________
//__________________________________________________________________________________________________________

// Function to get histogram
TH1 *GetCorrHisto(const char *filePattern, string rangeLabel) {
  TFile *f = new TFile(filePattern);
  string histName = Form("dphi_%s", rangeLabel.c_str());
  cout <<"[INFO] getting the histogram "<<histName<<endl;
  TH2D *hpx_2d = f->Get<TH2D>(histName.c_str());
  TH1D *hpx = (TH1D*) hpx_2d->ProjectionX(Form("%s_projX", histName.c_str()));
  return (TH1*) hpx;
}

//////////////////

