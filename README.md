# MPI correlations
Repository for code regarding Multiple Particle Interaction analysis in Run 3 with ALICE

## fitter
* InputToResults.C: master file
  - This file can be run with the following command: root -l -b -q InputToResults.C+'(justPlot, "inputFileName", "caseName")'
  - justPlot is a flag that can be turned on to skip the fitting procedure and just read the result file
  - inputFileName is the path to the root file with the phi distributions
  - caseName is the text input file indicator where all the initial parameters for each fit are saved
* SignalExtraction.C: extracts the results from one bin
* BuildPDF.C: builds the pdf model and imports the parameters to the RooWorkSpace used for the fits
* SetRangesAndLabels.C: includes many functions that help to extract/set ranges and labels
* input/input_caseName.txt
  - This file gives the fitter all the initial parameters for the fits. conventions to be respected:
    - parameters need to be seperated by ;
    - ptTrig;ptAssoc;mult give the bin edges with min-max format
    - fitStat is a flag to indicate which fits need to be done. Set to _todo_ to fit. if a fit results already exist for this bin, the results will be updated without changing, nor loosing, the other results
    - fitSig;fitBkg indicate the fitting models for the signal and background. for new models, please add them in BuildPDF.C
    - other paramters are initial values for the fits if one value is provided with [x] format, it is fixed. If 3 values are provided with [x,xmin,xmax] x will be the initial value and xmin-xmax will be the range
    - as many particles as needed can be added
