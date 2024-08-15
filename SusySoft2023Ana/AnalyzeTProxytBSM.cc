#define ANALYZETPROXYTBSM_cxx

#include "AnalyzeTProxytBSM.h"

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <vector>
#include <cstring>
#include <string>
#include <fstream>
#include"TGraphErrors.h"
#include"TGraphAsymmErrors.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
//#pragma link C++ class std::vector< std::vector >+; 
//#pragma link C++ class std::vector< TLorentzVector >+;
//#ifdef __MAKECINT__
//#pragma link C++ class NtupleVarsTProxy+;
//#endif

using namespace TMVA;
int main(int argc, char* argv[])
{

  if (argc < 6) {
    cerr << "Please give 5 arguments " << "runList " << " " << "outputFileName" << " " << "which year dataset" <<" "<<"which Process"<< " "<<"which Lostlep bkg"<< " "<<"Which pho_ID"<<endl;
    return -1;
  }
  const char *inputFileList = argv[1];
  const char *outFileName   = argv[2];
  const char *data          = argv[3];
  const char *sample=argv[4];
  const char *elec = argv[5];
  const char *phoID = argv[6];
  //TString pho_ID = phoID;

  AnalyzeTProxytBSM ana(inputFileList, outFileName, data,sample, elec,phoID);

  //=== === Loop over input files === ====
  int iFile = 0; 
  ifstream infile(inputFileList, ifstream::in);
  std::string buffer;
  if(!infile.is_open()) {
    std::cerr << "** ERROR: Can't open '" << inputFileList << "' for input" << std::endl;
    return kFALSE;
  }
  while(1) {
    infile >> buffer;
    if(!infile.good()) break;
    iFile++;
    std::cout << "===========================================================================" << std::endl;
    std::cout << "iFile " << iFile << "  Analyzing tree from " << buffer.c_str() << std::endl;
    std::cout << "===========================================================================" << std::endl;

    // for skimmed tree
    TFile *fin = new TFile(buffer.c_str());
    TTree *chain = (TTree*) fin->FindObjectAny("PreSelection");
    //TChain *chain = new TChain("PreSelection");
    //chain->Add(buffer.c_str());
    //TTree *chain; fin->GetObject("PreSelection",chain);
    std::cout << "main(): chain->GetEntries() "<<  chain->GetEntries() <<std::endl;    
    //ana.fChain = chain;
    //ana.fDirector.SetTree(chain);
    ana.Init(chain);
    ana.EventLoop(buffer.c_str());
    delete chain; 
    delete fin;
  }

  // === === some random summary === ===
  cout << "dataset " << data << " " << endl;
  cout<<"If analyzing the lost electron estimation ? "<<"  "<<elec<<endl;
  cout<<"Which pho_ID: "<<"\t"<<phoID<<endl;
  return 0;
}

//void AnalyzeLightBSM::EventLoop(const char *data,const char *inputFileList, const char *sample , const char *outFileName, const char *elec, const char* phoID) {
void AnalyzeTProxytBSM::EventLoop(std::string buffer) {
  
  std::cout << "AnalyzeTProxytBSM::EventLoop() " << std::endl;

  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  std::cout << "Analyzing " << buffer.c_str() << " nentries " << nentries << std::endl;  

  Long64_t nbytes = 0, nb = 0;
  int decade = 0;

  for (Long64_t jentry=0; jentry<fChain->GetEntries();jentry++) {

    fDirector.SetReadEntry(jentry);

    // == == print number of events done == == == == == == == =
    double progress = 10.0 * jentry / (1.0 * nentries);
    int k = int (progress);
    if (k > decade)
      cout << 10 * k << " %" << endl;
    //decade = k;

    //std::cout << jentry << " " << MET << std::endl;
    double wt = 0.12;
    h_MET->Fill(MET, wt);
    h_MET_test->Fill(1.0, wt);

    if (k > decade){
      std::cout << "h_MET->Integral()      " << h_MET->Integral()      << std::endl;
      std::cout << "h_MET_test->Integral() " << h_MET_test->Integral() << std::endl;
    }
    decade = k;
    
    //    if(jentry<10 ) {
      //std::cout<< "jentry " << jentry << " RunNum " << RunNum << std::endl;
      //std::cout << "GenParticles->size() "<< GenParticles->size() << std::endl;
      for(Long64_t ii=0; ii<GenParticles->size(); ii++){
	ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<float> > mygen = GenParticles[(int)ii];
	if(GenParticles[(int)ii].Pt()>1.0){
	h_gen_pT  ->Fill(GenParticles[(int)ii].Pt());
	h_gen_eta ->Fill(GenParticles[(int)ii].Eta());
	h_gen_phi ->Fill(GenParticles[(int)ii].Phi());
	}
	//std::cout <<" ii, Pt, Eta, Phi, E " << ii << " " << GenParticles[(int)ii].Pt() << " " << GenParticles[(int)ii].Eta() << " " << GenParticles[(int)ii].Phi() << " " << GenParticles[(int)ii].E() << " pdgid, parentid, status " << GenParticles_PdgId[(int)ii] << " " << GenParticles_ParentId[(int)ii] << " " << GenParticles_Status[(int)ii]  << std::endl;
      }
      
      //std::cout << std::endl; 
      //std::cout << "Electrons->size() "<< Electrons->size() << std::endl;
      for(Long64_t ii=0; ii<Electrons->size(); ii++){
	ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<float> > myele = Electrons[(int)ii];
	//std::cout <<" ii, Pt, Eta, Phi, E " << ii << " " << Electrons[(int)ii].Pt() 		  << " " << Electrons[(int)ii].Eta() << " " << Electrons[(int)ii].Phi() 		  << " " << Electrons[(int)ii].E()  		  << " iso, mediumID " << Electrons_iso[(int)ii] << " " << Electrons_mediumID[(int)ii]		  << std::endl;
	h_ele_pT  ->Fill(Electrons[(int)ii].Pt());
	h_ele_eta ->Fill(Electrons[(int)ii].Eta());
	h_ele_phi ->Fill(Electrons[(int)ii].Phi());
      }
      
      //std::cout << std::endl; 
      //std::cout << "Photons->size() "<< Photons->size() << std::endl;
      for(Long64_t ii=0; ii<Photons->size(); ii++){
	ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<float> > mypho = Photons[(int)ii];
	//std::cout <<" ii, Pt, Eta, Phi, E " << ii << " " << Photons[(int)ii].Pt() 		  << " " << Photons[(int)ii].Eta() << " " << Photons[(int)ii].Phi() 		  << " " << Photons[(int)ii].E()  		  << " mvavalueID, pfGammaIso " << Photons_mvaValuesID[(int)ii] << " " << Photons_pfGammaIso[(int)ii]		  << std::endl;
	h_pho_pT  ->Fill(Photons[(int)ii].Pt());
	h_pho_eta ->Fill(Photons[(int)ii].Eta());
	h_pho_phi ->Fill(Photons[(int)ii].Phi());
      }
      //      std::cout << "================================" << std::endl;
      //    } // if(jentry .. 

  }
}


myLV AnalyzeTProxytBSM::getBestPhoton(int pho_ID){
  //vector<TLorentzVector> goodPho;
  vector<myLV> goodPho;
  vector<int> goodPhoIndx;
  for(int iPho=0;iPho<Photons->size();iPho++){
    //if(((*Photons_hasPixelSeed)[iPho]<0.001) && ( (*Photons_fullID)[iPho]))
    if(((*Photons_hasPixelSeed)[iPho]<0.001) && ( (*Photons_fullID)[iPho] && ((*Photons_hasPixelSeed)[iPho]<0.001) &&( pho_ID==0 || (pho_ID==1 &&(((*Photons_cutBasedID)[iPho]==1 || (*Photons_cutBasedID)[iPho]==2))) || (pho_ID==2 && (*Photons_cutBasedID)[iPho]==2) || (pho_ID==3 && (*Photons_mvaValuesID)[iPho]>-0.02) || (pho_ID==4 && (*Photons_mvaValuesID)[iPho]>0.42))) ) 
      {
	goodPho.push_back(Photons[iPho] );
	goodPhoIndx.push_back(iPho);
      }
  }
  
  int highPtIndx=-100;
   for(int i=0;i<goodPho.size();i++){
     if(i==0) highPtIndx=0;
     else if( (goodPho[highPtIndx].Pt()) < (goodPho[i].Pt()) ){highPtIndx=i;}
   }
   
   if(highPtIndx>=0){
     bestPhotonIndxAmongPhotons = goodPhoIndx[highPtIndx];
   }
   else bestPhotonIndxAmongPhotons = -100;
   if(highPtIndx==-100){myLV v0;return v0;}
   else return goodPho[highPtIndx];
   
}
  
//SS//== not using this functions == 
Bool_t AnalyzeTProxytBSM::Process(Long64_t entry) {

  std::cout << entry << std::endl;
   fDirector.SetReadEntry(entry);
   std::cout<< "entry " << entry << " RunNum " << RunNum << std::endl;
   std::cout << "GenParticles->size() "<< GenParticles->size() << std::endl;
  return 0;
}

