#ifndef ANALYZETPROXYTBSM_H
#define ANALYZETPROXYTBSM_H
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "NtupleVarsTProxy.h"
#include "TH1F.h"
#include "TTree.h"
#include "TH2.h"
#include "TProfile.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TDirectory.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"

//#pragma link C++ class std::vector< std::vector >+; 
//#pragma link C++ class std::vector< TLorentzVector >+;
//#pragma link C++ class NtupleVarsTProxy+;

//Define Root LorentzVector
typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<float> > myLV;

class AnalyzeTProxytBSM : public NtupleVarsTProxy{

 public:
  AnalyzeTProxytBSM(const TString &inputFileList="foo.txt", const char *outFileName="histo.root",const char *dataset="data", const char *sample="sample", const char* LostlepFlag ="Flag", const char* phoID="phoID");
  ~AnalyzeTProxytBSM();
  //void     EventLoop(const char *,const char *,const char *,const char *, const char*, const char*);
  //  void   EventLoop(std::string buffer);
  void   EventLoop(std::string buffer,const char *,const char *);

  void   BookHistogram(const char *);//, const char *);
  Bool_t Process(Long64_t entry);
  myLV   getBestPhoton(int);
  double getGenLep(myLV);
  int    bestPhotonIndxAmongPhotons=-100;
  void   CrossSection_Map_Init();
  void   DoCutFlow(int k, int decade, bool printsummary=false, bool debug=false);
  bool   RemoveSampleOverlap(TString s_sample, myLV bestPhoton);
  
  enum evtSel {selNone=0, selPho=1, selMET=2, selNJet=3, selST=4, selTrgEff=5, selEvtCln=6,selJetMetPhi=7, selEMu=8, selIsoTrk=9, selRmOverlap=10};
  std::vector<double> NEvtSel;//(10,0.0);
  double EvtWeight;
  //trigger efficiency weights
  const double p0=1.787e+02,p1=6.657e+01,p2=9.47e-01;

  bool genphocheck=false;
  int  genphomatch_before=0, genphomatch_after=0;

  
  TFile *oFile;
  TH1D *h_MET_test;
  TH1D *h_MET;
  TH1F *h_ele_pT, *h_ele_eta, *h_ele_phi;
  TH1F *h_pho_pT, *h_pho_eta, *h_pho_phi;
  TH1F *h_gen_pT, *h_gen_eta, *h_gen_phi;
};
#endif


#ifdef ANALYZETPROXYTBSM_cxx

//void AnalyzeLightBSM::BookHistogram(const char *outFileName, const char *N2_mass) {
void AnalyzeTProxytBSM::BookHistogram(const char *outFileName) {
  std::cout << "AnalyzeLightBSM::BookHistogram " << std::endl;

  oFile = new TFile(outFileName, "recreate");
  oFile->cd();

  // Book your histograms & summary counters here 
  h_MET_test= new TH1D("h_MET_test", "h_MET_test", 100, 0.0, 4000.0);
  h_MET     = new TH1D("h_MET", "h_MET", 100, 0.0, 4000.0);
  h_ele_pT  = new TH1F("h_ele_pT",  "h_ele_pT", 100, 0.0, 1000.0);
  h_ele_eta = new TH1F("h_ele_eta", "h_ele_eta",100, -3.0, 3.0);
  h_ele_phi = new TH1F("h_ele_phi", "h_ele_phi",100, -3.2, 3.2);
  h_pho_pT  = new TH1F("h_pho_pT",  "h_pho_pT", 100, 0.0, 1000.0);
  h_pho_eta = new TH1F("h_pho_eta", "h_pho_eta",100, -3.0, 3.0);
  h_pho_phi = new TH1F("h_pho_phi", "h_pho_phi",100, -3.2, 3.2);
  h_gen_pT  = new TH1F("h_gen_pT",  "h_gen_pT", 100, 0.0, 1000.0);
  h_gen_eta = new TH1F("h_gen_eta", "h_gen_eta",100, -6.0, 6.0);
  h_gen_phi = new TH1F("h_gen_phi", "h_gen_phi",100, -3.2, 3.2);

}

AnalyzeTProxytBSM::AnalyzeTProxytBSM(const TString &inputFileList, const char *outFileName,const char *dataset, const char *sample, const char* LostlepFlag, const char* phoID) {

  NEvtSel.assign(11, 0.0);
  EvtWeight = 1.0;
  
  std::cout << outFileName << std::endl;

  string nameData=dataset;//vvv

  if(nameData!="signalH") nameData="BG";
  if(nameData=="signalH") nameData="signal";
  cout<<"Treating the input files as "<<nameData<<" for setting tree branches"<<endl;
  BookHistogram(outFileName); //, N2_mass);
  CrossSection_Map_Init();
}

void AnalyzeTProxytBSM::CrossSection_Map_Init() {
  char *f_name_EH = new char[2000];
  sprintf(f_name_EH,"./map_crosssection_SMprocess_v1.txt");//,chi2_method);
  std::ifstream in_EH(f_name_EH);
  if(!in_EH) {
    cout<<"ERROR => "<<f_name_EH<<" Not found"<<endl;
    //return;
    exit(0);
  }
  string process_name;
  float value, entries;
  cout<<"File name = "<<f_name_EH<<endl;
  while(in_EH>>process_name>>value>>entries){
    std::pair<std::string, float> temp_pair;    
    float weight =value/entries;
    cout << "values: " << value << " entries: " << entries << endl;
    temp_pair = std::make_pair(process_name,weight);
    cross_sectionValues.insert(temp_pair);
  }
}


AnalyzeTProxytBSM::~AnalyzeTProxytBSM() { 
  if (!fChain) return;
  delete fChain->GetCurrentFile();
  oFile->cd();
  oFile->Write();
  oFile->Close();
}

#endif // AnalyzeTProxytBSM_cxx
