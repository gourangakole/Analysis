#ifndef _histogramplotter_h_ 
#define _histogramplotter_h_ 
 
#if !defined(__CINT__) || defined(__MAKECINT__) 
 
#include <iomanip> 
#include <iostream> 
#include <fstream> 
#include <map>
 
#include "TRandom2.h" 
#include "TMatrixD.h" 
#include "TF1.h" 
#include "TProfile.h" 
#include "TObjArray.h" 
#include "TMatrixD.h" 
#include "TH1.h" 
#include "TH2.h"
#include "TTimeStamp.h" 
#include "TFile.h"
#include <exception> 
 
#endif 
 
class HistogramPlotter{ 

public :
  HistogramPlotter(){}
  virtual ~HistogramPlotter(){}
  ///~HistogramPlotter(){}

  void InitHist(TString dirname, TString parentDir, TFile *file, bool inNvtx);
  void addHisto(TString name, TString dirname, int range, double min, double max);
  void add2DHisto(TString name, TString dirname, int range1, double min1, double max1, int range2, double min2, double max2);
  void fillHisto(TString name, TString dirname, double value, double weight=1.0);
  TH1* getHisto(TString name, TString dirname);
  void CreateAnalHistos(TString flowType, TFile* outFile_);

private :
  std::map<TString, TH1*> histos1_;
  std::map<TString, TH2*> histos2_;
  //static std::map<TString, histos1_> histDir1_;
  //static std::map<TString, histos2_> histDir2_;
  //std::map<TString, TDirectory*> hDir_;

  ClassDef(HistogramPlotter, 1)
};
#endif



