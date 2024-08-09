// Author: Batoul Diab
#include <iostream>
#include <stdio.h>
#include "TString.h"


bool parseFile(string FileName, vector< map<string, string> >& data);
bool parseString(string input, string delimiter, vector<double>& output);
bool readFile(string FileName, vector< vector<string> >& content, const int nCol=-1, int nRow=-1);
bool setParameters(map<string, string> row, struct KinCuts& cut, map<string, string>& parIni);
bool addParameters(string InputFile,  vector< struct KinCuts >& cutVector, vector< map<string, string> >&  parIniVector);

string varFancyLabel(string varLabel) {
  std::map<std::string, std::string> fancyLabels;

  fancyLabels["f_ns"] = "N_{ns}";
  fancyLabels["f_fs"] = "N_{fs}";
  fancyLabels["f_b"] = "N_{bkg}";
  fancyLabels["f_ns_1"] = "f_{ns}";
  fancyLabels["mean_ns"] = "mean_{ns}";
  fancyLabels["sigma_ns"] = "#sigma^{1}_{ns}";
  fancyLabels["sigmaRatio_ns"] = "(#sigma_{2}/#sigma_{1})_{ns}";
  fancyLabels["mean_fs"] = "mean_{fs}";
  fancyLabels["sigma_fs"] = "#sigma_{fs}";
  fancyLabels["norm_bkg"] = "c_{bkg}";
  
  if (fancyLabels.find(varLabel) != fancyLabels.end())
    return fancyLabels[varLabel]; 
  else return "";
}

//___________________________________________________________________________________________________________
//___________________________________________________________________________________________________________

typedef struct StartEnd {
  int Start, End;
} StartEnd;

typedef struct MinMax {
  double Min, Max;
} MinMax;

typedef struct KinCuts {
  StartEnd   mult;
  MinMax ptTrig ;
  MinMax ptAssoc;
} KinCuts;

bool isEqualKinCuts(struct KinCuts cutA, struct KinCuts cutB, bool isPbPb) 
{
  bool cond = true;
  cond = cond && (cutA.mult.Start    == cutB.mult.Start);
  cond = cond && (cutA.mult.End      == cutB.mult.End);

  cond = cond && (cutA.ptTrig.Min    == cutB.ptTrig.Min);
  cond = cond && (cutA.ptTrig.Max      == cutB.ptTrig.Max);
  
  cond = cond && (cutA.ptAssoc.Min    == cutB.ptAssoc.Min);
  cond = cond && (cutA.ptAssoc.Max      == cutB.ptAssoc.Max);

  return cond;
    
}

//___________________________________________________________________________________________________________
//___________________________________________________________________________________________________________
//functions to handle strings (mainly for input parameters)

std::vector<std::string> splitString(const std::string& str, char delimiter) {
    std::vector<std::string> result;
    std::stringstream ss(str);
    std::string token;
    while (std::getline(ss, token, delimiter)) {
        result.push_back(token);
    }
    return result;
}

std::string removeSurroundingBrackets(const std::string& str) {
    if (str.size() >= 2 && str.front() == '[' && str.back() == ']') {
        return str.substr(1, str.size() - 2);
    }
    return str;
}


bool addParameters(string InputFile,  vector< struct KinCuts >& cutVector, vector< map<string, string> >&  parIniVector)
{
  vector< map<string, string> >  data;
  if(!parseFile(InputFile, data)) { return false; }
  
    for(vector< map<string, string> >::iterator row=data.begin(); row!=data.end(); ++row) {
      struct KinCuts cut; map<string, string> parIni;
      if(!setParameters(*row, cut, parIni)) { return false; }
      cutVector.push_back(cut);  parIniVector.push_back(parIni);
    }
  return true;
};

bool setParameters(map<string, string> row, struct KinCuts& cut, map<string, string>& parIni)
{
  cut.ptTrig.Min = 0.;
  cut.ptTrig.Max = 150.;

  cut.ptAssoc.Min = 0.;
  cut.ptAssoc.Max = 150.;

  cut.mult.Start = 0;
  cut.mult.End = 100;
  
  // set parameters from file
  for(map<string, string>::iterator col=row.begin(); col!=row.end(); ++col) {

    string label = col->first;
    //cout <<"checking parameter "<<col->first<<", with value"<<col->second<<endl;
    if (label=="ptTrig") {
      if (col->second=="" || col->second.find("-")==std::string::npos) {
        cout << "[ERROR] Input column 'ptTrig' has invalid value: " << col->second << endl; return false;
      }  
      std::vector<double> v; 
      if(!parseString(col->second, "-", v)) { return false; }
      if (v.size()!=2) {
        cout << "[ERROR] Input column 'ptTrig' has incorrect number of values, it should have 2 values but has: " << v.size() << endl; return false;
      }  
      cut.ptTrig.Min = v.at(0); 
      cut.ptTrig.Max = v.at(1);
    } 
    else if (label=="ptAssoc"){
      if (col->second=="" || col->second.find("-")==std::string::npos) {
        cout << "[ERROR] Input column 'ptAssoc' has invalid value: " << col->second << endl; return false;
      }
      std::vector<double> v; 
      if(!parseString(col->second, "-", v)) { return false; }
      if (v.size()!=2) {
        cout << "[ERROR] Input column 'ptAssoc' has incorrect number of values, it should have 2 values but has: " << v.size() << endl; return false;
      }  
      cut.ptAssoc.Min = v.at(0); 
      cut.ptAssoc.Max = v.at(1);
    } 
    else if (label=="mult"){
        if (col->second=="" || col->second.find("-")==std::string::npos) {
          cout << "[ERROR] Input column 'mult' has invalid value: " << col->second << endl; return false;
        }
        std::vector<double> v;
        if(!parseString(col->second, "-", v)) { return false; }
        if (v.size()!=2) {
          cout << "[ERROR] Input column 'mult' has incorrect number of values, it should have 2 values but has: " << v.size() << endl; return false;
        }
        cut.mult.Start = (int) (v.at(0));
        cut.mult.End = (int) (v.at(1));
    }
    else if (label.find("fit")!=std::string::npos){ //fitStat, fitSig, fitBkg
      if (col->second=="") {
        cout << "[ERROR] Input column "<<label<<" has empty value" << endl; return false;
      }
      parIni[col->first] = col->second;
    } 
    else {
      if (col->second != "") {
	string value = col->second;
	// check that initial parameters format is correct: [ num, num, num ]
	if ((value.find("[")==std::string::npos)||(value.find("]")==std::string::npos)) {
	  // Special cases like parameter constrains could be set here but for now, let's keep it simple
	  cout << "[ERROR] Either ']' or '[' are missing in the initial parameter values for " <<label<<"!"<< endl; return false;
	} else {
	  value.erase(value.find("["), string("[").length());
	  value.erase(value.find("]"), string("]").length());
	}
	std::vector<double> v; 
	if(!parseString(value, ",", v)) { return false; }
        if (v.size()>3 || v.size()<1) {
          cout << "[ERROR] Initial parameter " << col->first << " has incorrect number of values, it has: " << v.size() << endl; return false;
        } 
	// everything seems alright, then proceed to save the values
	if (v.size()==1){
	  // if only one value is given i.e. [ num ], consider it a constant value
	  parIni[col->first] = Form("[%.6f]", v.at(0));
	}
	else
	  {
	    parIni[col->first] = col->second;
	  }
      } else {
        parIni[col->first] = "";
      }
    }
  }
  return true;
};

bool parseString(string input, string delimiter, vector<double>& output)
{
  // remove spaces from input string 
  input.erase(std::remove(input.begin(), input.end(), ' '), input.end());
  // proceed to parse input string
  char *end;
  while(input!="") {
    double d = strtod(input.c_str(), &end);
    if (end != input) {
      output.push_back(d);
    } else {
      cout << "[ERROR] The conversion from string to double failed!"; return false;
    }
    input = end; 
    if(input.find(delimiter.c_str())!= std::string::npos){ input.erase(input.find(delimiter.c_str()), delimiter.length()); }
  }
  return true;
};


bool parseFile(string FileName, vector< map<string, string> >& data)
{
  vector< vector<string> > content, tmp; 
  if(!readFile(FileName, tmp, -1, 1)){ return false; }
  vector<string> header = tmp.at(0);
  if (header.size()==0) { cout << "[ERROR] The header is null!" << endl; return false; }
  if(!readFile(FileName, content, header.size())){ return false; }
  for(vector<string>::iterator rHeader=header.begin(); rHeader!=header.end(); ++rHeader) {
    if (*rHeader=="") { cout << "[ERROR] A column has no label!" << endl; return false; }
  }

  for(vector< vector<string> >::iterator row=content.begin()+1; row!=content.end(); ++row) {
    map<string, string> col;
    for (unsigned int i=0; i<header.size(); i++) {
      if (i<row->size()) {
	col[header.at(i)] = row->at(i);
      } else {
	col[header.at(i)] = "";
      }
    }
    data.push_back(col);
  }

  return true;
};	

bool readFile(string FileName, vector< vector<string> >& content, const int nCol, int nRow)
{
  if (nCol==0 || nRow==0) { 
    cout << "[WARNING] Ignoring content of File: " << FileName << endl; return true; 
  }
  if (nRow!=1) { cout << "[INFO] Reading file: " << FileName << endl; }
  ifstream myfile(FileName.c_str());
  if (myfile.is_open()){ 
    string line;
    while ( getline(myfile, line) ){
      if (nRow==0) break; else {nRow=nRow-1;}
      stringstream row(line);
      vector<string> cols; int i=0;
      while (true){
	string col; getline(row, col, ';');
	if ( (nCol>=0) ? (i>=nCol) : (col=="") ){ break; }
	cols.push_back(col);
	i++;
      }
      content.push_back(cols);
    }
  } else {
    cout << "[ERROR] File: " << FileName << " was not found!" << endl; return false;
  }
  return true;
};
