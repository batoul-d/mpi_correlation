// Author: Batoul Diab
#include <vector>
#include <iostream>
#include <algorithm>
#include <stdio.h>
#include <TSystem.h>
#include <TFile.h>
#include <TH1.h>
#include <TH3.h>
#include <TCanvas.h>
#include <TString.h>
#include <TMath.h>
#include <TLatex.h>
#include <TDirectory.h>
#include <TROOT.h>
#include <TKey.h>

#include "SignalExtraction.C"

void signalExtraction(const char *inputFile = "/home/bdiab/alice/mpi_correlation/dphi_corr.root", const char *caseName = "gaussian_uniform_all");
void plotResult(const char *caseName = "gaussian_uniform_all", string axisName = "ptTrig", int incMinMult=0, int incMaxMult=99, float incMinPTtrig=0.4, float incMaxPTtrig=150., float incMinPTassoc=0.4, float incMaxPTassoc=50.);

void InputToResults(bool justPlot = false, const char *inputFile = "/home/bdiab/alice/mpi_correlation/dphi_corr.root", const char *caseName = "gaussian_uniform_all") {
  
  if (!justPlot) signalExtraction(inputFile, caseName);
  
  int minMult = 0;
  int maxMult = 99;
  float minPTtrig = 0.7;
  float maxPTtrig = 150.;
  float minPTassoc = 0.4;
  float maxPTassoc = 0.7;
  string axisName = "multiplicity";
  plotResult(caseName, axisName.c_str(), minMult, maxMult, minPTtrig, maxPTtrig, minPTassoc, maxPTassoc);
}


void signalExtraction(const char *inputFile, const char *caseName) {

  string outputName = Form("output/output_%s.root",caseName);
  
  //do the fits or at least some of them
  vector< struct KinCuts >       cutVector;
  vector< map<string, string> >  parIniVector;
  vector< map<string, double> >  allResults;
  
  if (!addParameters(Form("input/input_%s.txt",caseName), cutVector, parIniVector)) { return; }
  std::ifstream fileCheck(outputName.c_str());
  bool fileExists = fileCheck.good();
  fileCheck.close();
  TFile* fSave;
  
  if (fileExists) 
    fSave = TFile::Open(outputName.c_str(),"UPDATE");
  else
    fSave = TFile::Open(outputName.c_str(),"RECREATE");
  if (!fSave || fSave->IsZombie()) {
    cout<<"[ERROR] Problem with the output file "<<outputName<<endl;
    return;
  }
  for (uint j = 0; j < cutVector.size(); j++) {
    if (parIniVector[j]["fitStat"] != "todo") {
      cout<<"[INFO] This fit does not have the 'todo' status, so it will be skipped!"<<endl;
      continue;
    }
    
    std::string rangeLabel = Form("ptTrig_%dp%d_%dp%d_ptAssoc_%dp%d_%dp%d_mult_%d_%d",
				  (int) cutVector[j].ptTrig.Min,
				  (int) (cutVector[j].ptTrig.Min*10)%10,
				  (int) cutVector[j].ptTrig.Max,
				  (int) (cutVector[j].ptTrig.Max*10)%10,
				  (int) cutVector[j].ptAssoc.Min,
				  (int) (cutVector[j].ptAssoc.Min*10)%10,
				  (int) cutVector[j].ptAssoc.Max,
				  (int) (cutVector[j].ptAssoc.Max*10)%10,
				  (int) cutVector[j].mult.Start,
				  (int) cutVector[j].mult.End);

    RooWorkspace* ws = new RooWorkspace(Form("ws_%s", rangeLabel.c_str()));
    map<string, double> resultsFit;
      
    resultsFit.clear();
    resultsFit = SignalExtraction1Fit(inputFile, cutVector[j].mult.Start, cutVector[j].mult.End, cutVector[j].ptTrig.Min, cutVector[j].ptTrig.Max, cutVector[j].ptAssoc.Min, cutVector[j].ptAssoc.Max, parIniVector[j], ws, rangeLabel);
    resultsFit["multMin"] = cutVector[j].mult.Start;
    resultsFit["multMax"] = cutVector[j].mult.End;
    resultsFit["ptTrigMin"] = cutVector[j].ptTrig.Min;
    resultsFit["ptTrigMax"] = cutVector[j].ptTrig.Max;
    resultsFit["ptAssocMin"] = cutVector[j].ptAssoc.Min;
    resultsFit["ptAssocMax"] = cutVector[j].ptAssoc.Max;
    
    allResults.push_back(resultsFit);

    //transform the results into trees
    fSave->cd();
    TTree* resTree = new TTree(Form("tree_%s_new", rangeLabel.c_str()), "Tree of results");
    
    std::map<std::string, double*> branches;
    
    // Create branches for each map entry
    for (const auto& entry : resultsFit) {
        branches[entry.first] = new double;
        resTree->Branch(entry.first.c_str(), branches[entry.first]);
    }

    // Fill the resTree with data from the map
    for (const auto& entry : resultsFit) {
        *(branches[entry.first]) = entry.second;
	cout <<"filling tree with "<<entry.first<<" = "<<entry.second<<endl;
    }
    resTree->Fill(); //just fill one entry in the tree (make a different tree for each fit to allow updating)
    
    //check if the ws exists, if yes delete it and save it again
    TObject* obj = fSave->Get(Form("ws_%s", rangeLabel.c_str()));
    if (obj) {
      fSave->Delete(Form("ws_%s;*",rangeLabel.c_str()));
      //fSave->Delete(Form("hist_%s;*",rangeLabel.c_str()));
      fSave->Delete(Form("tree_%s;*",rangeLabel.c_str()));
    }
      
    fSave->cd();
    ws->Write();
    resTree->Write(Form("tree_%s", rangeLabel.c_str()));
  }
  fSave->Close();
}


void plotResult(const char *caseName, string axisName, int incMinMult, int incMaxMult, float incMinPTtrig, float incMaxPTtrig, float incMinPTassoc, float incMaxPTassoc){
  cout <<"[INFO] The plotting function is yet to be done"<<endl;
  gStyle->SetOptStat(0);
  
  map<string, double> resultsFit;
  double *binEdges = new double[100];
  int nbins;
  double *resVal = new double[100];
  double *resVal_ns = new double[100];
  double *resVal_fs = new double[100];
  double *resErr = new double[100];
  double *resErr_ns = new double[100];
  double *resErr_fs = new double[100];
  
  string fileName = Form("output/output_%s.root",caseName);
  TFile* fHist = TFile::Open(fileName.c_str(),"READ");
  if (!fHist || fHist->IsZombie()) {
    cout<<"[ERROR] Problem with the result file that contains all the histograms"<<fileName<<endl;
    return;
  }
    
  cout<<"[INFO] Opening the result input file"<<endl;
  TIter next(fHist->GetListOfKeys());
  TKey *key; int iBin = 0; double lastEdge = 0.;
  while ((key = (TKey*)next())) {
    // Get the class name of the object
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (!cl) continue;
    if (cl->InheritsFrom("TTree")) {
      TTree *resTree = (TTree*) key->ReadObj();//fHist->Get(key->GetName());
      
      cout<<"[INFO] Reading the input tree"<<resTree->GetName()<<endl;

      double multMin = 0.; resTree->SetBranchAddress("multMin", &multMin);
      double multMax = 0.; resTree->SetBranchAddress("multMax", &multMax);
      double ptTrigMin = 0.; resTree->SetBranchAddress("ptTrigMin", &ptTrigMin);
      double ptTrigMax = 0.; resTree->SetBranchAddress("ptTrigMax", &ptTrigMax);
      double ptAssocMin = 0.; resTree->SetBranchAddress("ptAssocMin", &ptAssocMin);
      double ptAssocMax = 0.; resTree->SetBranchAddress("ptAssocMax", &ptAssocMax);
      
      double f_ns = 0.; resTree->SetBranchAddress("f_ns", &f_ns);
      double f_ns_err = 0.; resTree->SetBranchAddress("f_ns_err", &f_ns_err);
      double f_fs = 0.; resTree->SetBranchAddress("f_fs", &f_fs);
      double f_fs_err = 0.; resTree->SetBranchAddress("f_fs_err", &f_fs_err);

      resTree->GetEntry(0);
      
      //cout<<"[INFO] just before the if multiplicity f_ns = "<<f_ns<<endl;
      if (axisName.find("multiplicity")!=std::string::npos) { //if plotting as function of multiplicity, the multiplicity of the result needs to be in the range but the pt need to be the exact edges
	
        if (multMin < incMinMult || multMax > incMaxMult)
          continue;
        if (fabs(ptTrigMin - incMinPTtrig)>0.00001 || fabs(ptTrigMax - incMaxPTtrig)>0.00001)
          continue;
        if (fabs(ptAssocMin - incMinPTassoc)>0.00001 || fabs(ptAssocMax - incMaxPTassoc)>0.00001)
          continue;
      
	cout<<"[INFO] This tree passed the bin selection, getting the results for ptTrig ["<<ptTrigMin<<"-"<<ptTrigMax<<"], ptAssoc ["<<ptAssocMin<<"-"<<ptAssocMax<<"], mult ["<<multMin<<"-"<<multMax<<"]"<<endl;
	binEdges[iBin] = multMin;
	lastEdge =multMax+1;
	//cout<<"[INFO] binEdge = "<<binEdges[iBin]<<endl;
      }//end of if multiplicity
      else if (axisName.find("pTtrig")!=std::string::npos) {
	binEdges[iBin] = ptTrigMin;
	lastEdge = ptTrigMax;
      }//end of if pTtrig
      else if (axisName.find("pTassoc")!=std::string::npos) {
	binEdges[iBin] = ptAssocMin;
	lastEdge = ptAssocMax;
      }//end of pTassoc

      resVal_ns[iBin] = f_ns;
      resErr_ns[iBin] = f_ns_err;
      resVal_fs[iBin] = f_fs;
      resErr_fs[iBin] = f_fs_err;
      resVal[iBin] = f_ns+f_fs;
      resErr[iBin] = sqrt(f_ns_err*f_ns_err+f_fs_err*f_fs_err);

      //cout <<"[INFO] for iBin = "<<iBin<<", resVal_ns[iBin] = "<<resVal_ns[iBin]<<", resErr_ns[iBin] = "<<resErr_ns[iBin]<<endl;
      iBin++;
    }//end of if TTree
  }//end of while loop on file components

  //cout<<"iBin"
  binEdges[iBin] = lastEdge;
  
  TCanvas *c1 = new TCanvas("c1", "", 800, 800);
  c1->cd();
    
  nbins = iBin;
  TH1D* results = new TH1D(Form("results_%s", axisName.c_str()), "", nbins, binEdges);

  results->SetLineColor(kBlue+2);
  results->SetLineWidth(2);
  results->SetMarkerStyle(kFullCircle);
  results->SetMarkerColor(kBlue+1);
  
  results->GetYaxis()->SetLabelSize(0.04);
  results->GetYaxis()->SetTitleSize(0.04);
  results->GetYaxis()->SetTitleOffset(1.2);
  results->GetYaxis()->SetTitleFont(42);
  results->GetYaxis()->CenterTitle(kTRUE);
  
  results->GetXaxis()->SetLabelSize(0.04);
  results->GetXaxis()->CenterTitle(kTRUE);
  results->GetXaxis()->SetTitleSize(0.04);
  results->GetXaxis()->SetTitleFont(42);
  results->GetXaxis()->SetTitleOffset(1.);
  if (axisName.find("multiplicity")!=std::string::npos) {
    results->GetXaxis()->SetTitle("Multiplicity"); 
  }
  else if (axisName.find("ptTrig")!=std::string::npos) {
    results->GetXaxis()->SetTitle("p_{T, trig}"); 
  }
  else if (axisName.find("ptAssoc")!=std::string::npos) {
    results->GetXaxis()->SetTitle("p_{T, assoc}"); 
  }

  TH1D* results_ns = (TH1D*) results->Clone(Form("results_ns_%s", axisName.c_str()));
  TH1D* results_fs = (TH1D*) results->Clone(Form("results_fs_%s", axisName.c_str()));

  results->GetYaxis()->SetTitle("<N_{assoc}>");
  results->GetYaxis()->SetRangeUser(0,*std::max_element(resVal, resVal+nbins)*1.2);
  results_ns->GetYaxis()->SetTitle("<N_{assoc, near side}>");
  results_ns->GetYaxis()->SetRangeUser(0,*std::max_element(resVal_ns, resVal_ns+nbins)*1.2);
  results_fs->GetYaxis()->SetTitle("<N_{assoc, away side}>");
  results_fs->GetYaxis()->SetRangeUser(0,*std::max_element(resVal_fs, resVal_fs+nbins)*1.2);
  //results->Draw();
  
  for (int j=0; j<nbins; j++) {
    //cout<<"[INFO] filling the histogram at x = "<<binEdges[j]<<", with res = "<<resVal_ns[j]<<", err = "<<resErr_ns[j]<<endl;  
    int jBin = results->FindFixBin((binEdges[j]+binEdges[j+1])/2);
    results->SetBinContent(jBin, resVal[j]);
    results->SetBinError(jBin, resErr[j]);
    results_ns->SetBinContent(jBin, resVal_ns[j]);
    results_ns->SetBinError(jBin, resErr_ns[j]);
    results_fs->SetBinContent(jBin, resVal_fs[j]);
    results_fs->SetBinError(jBin, resErr_fs[j]);
    //cout<<"[INFO] the histogram at x = "<<results_ns->GetBinLowEdge(jBin)<<" is filled with "<<results_ns->GetBinContent(jBin)<<" with error "<<results_ns->GetBinError(jBin)<<endl;
  }

  string rangeName = Form("ptTrig_%dp%d_%dp%d_ptAssoc_%dp%d_%dp%d_mult_%d_%d",
			  (int) incMinPTtrig,
			  (int) (incMinPTtrig*10)%10,
			  (int) incMaxPTtrig,
			  (int) (incMaxPTtrig*10)%10,
			  (int) incMinPTassoc,
			  (int) (incMinPTassoc*10)%10,
			  (int) incMaxPTassoc,
			  (int) (incMaxPTassoc*10)%10,
			  (int) incMinMult,
			  (int) incMaxMult);
  
  results->Draw("E1");
  c1->SaveAs(Form("output/results_assoc_%s_%s.pdf",caseName,rangeName.c_str()));
  results_ns->Draw("E1");
  c1->SaveAs(Form("output/results_assoc_NS_%s_%s.pdf",caseName,rangeName.c_str()));
  results_fs->Draw("E1");
  c1->SaveAs(Form("output/results_assoc_FS_%s_%s.pdf",caseName,rangeName.c_str()));
  fHist->Close();
}
