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
#include "AnalyzeCutflow.cc"

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
    std::cout << "main(): chain->GetEntries() "<<  chain->GetEntries() <<std::endl;    
    ana.Init(chain);
    ana.EventLoop(buffer.c_str(),data,sample);
    //ana.EventLoop(buffer.c_str());
    delete chain; 
    delete fin;

    // method below worked for older prelegacy ra2b ntuples. Need to call event loop for each file now.
    //TChain *chain = new TChain("PreSelection");
    //chain->Add(buffer.c_str());
    //TTree *chain; fin->GetObject("PreSelection",chain);
    //ana.fChain = chain;
    //ana.fDirector.SetTree(chain);
  }

  // === === some random summary === ===
  cout << "dataset " << data << " " << endl;
  cout<<"If analyzing the lost electron estimation ? "<<"  "<<elec<<endl;
  cout<<"Which pho_ID: "<<"\t"<<phoID<<endl;
  return 0;
}

//void AnalyzeLightBSM::EventLoop(const char *data,const char *inputFileList, const char *sample , const char *outFileName, const char *elec, const char* phoID) {
//void AnalyzeTProxytBSM::EventLoop(std::string buffer) {
void AnalyzeTProxytBSM::EventLoop(std::string buffer, const char *data, const char *sample) {
  
  std::cout << "AnalyzeTProxytBSM::EventLoop() " << std::endl;

  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  std::cout << "Analyzing " << buffer.c_str() << " nentries " << nentries << std::endl;  

  Long64_t nbytes = 0, nb = 0;
  int decade = 0;

  //<-- get sample specific cross-section -->
  char* s_cross = new char[100];
  sprintf(s_cross,"%s.%s",data,sample);
  std::string s_process = s_cross;
  TString s_Process = s_process;
  TString s_sample= sample;
  std::cout << "s_sample " << s_sample << std::endl;


  double cross_section = getCrossSection(s_process);
  /*
  double cross_section = getCrossSection(s_process);
  double wt = -999.0;
  std::cout << "wt " << wt << std::endl;
  std::cout << cross_section << "\t" <<"analyzed process"<<"\t"<<s_cross<<endl;
  */
  std::cout << cross_section << "\t" <<"analyzed process"<<"\t"<<s_cross<<endl;

  //std::cout << "Pass_Pho_pT  "<< Pass_Pho_pT << " " <<"Pass_MET "<< Pass_MET<<std::endl;

  for (Long64_t jentry=0; jentry<fChain->GetEntries();jentry++) {
   
    fDirector.SetReadEntry(jentry);

    // == == print number of events done == == == == == == == =
    double progress = 10.0 * jentry / (1.0 * nentries);
    int k = int (progress);
    if (k > decade)
      cout << 10 * k << " %" << endl;
    //decade = k;

    //get event weight (xsec*lumi)
    EvtWeight = getEventWeight(s_process, cross_section);
    double wt = EvtWeight;
    if(jentry<10) std::cout << "wt " << wt << std::endl;
    

    // define analysis objects
    //get best photon
    int    pho_ID=0;  //for simplicity taking only soft one
    myLV bestPhoton=getBestPhoton(pho_ID);

    //remove jet-photon overlap
    vector<myLV> hadJets, bjets;
    vector<int>  jetMatchindx;
    int    hadJetID = -1;
    int    NHadJets = 0;
    double minDR=99999;
    int    minDRindx=-1, bJet1Idx=-1, phoMatchingJetIndx=-1; 
    double deepCSVvalue = 0.4148; 
    double dPhi_METjet1, dPhi_METjet2;
    float ST=0;
    float Jets_pT_Sum=0;
    int NEMu = NElectrons + NMuons;

    for(int i=0;i<Jets->size();i++){
      if( (Jets[i].Pt() > 30.0) && (abs(Jets[i].Eta()) <= 2.4) ){
	if (Photons->size()!=0) {
	  double dR=DeltaR(bestPhoton.Eta(),bestPhoton.Phi(),Jets[i].Eta(),Jets[i].Phi());
	  if(dR<minDR) {minDR=dR;minDRindx=i;}
	}
      } 
    }
    for(int i=0;i<Jets->size();i++){
      if( (Jets[i].Pt() > 30.0) && (abs(Jets[i].Eta()) <= 2.4) ){	  
	if( !(minDR < 0.3 && i==minDRindx) ){		
	  hadJetID= (*Jets_ID)[i];
	  if(hadJetID) {
	    hadJets.push_back(Jets[i]);
	    if((*Jets_bJetTagDeepCSVBvsAll)[i] > deepCSVvalue){
	      bjets.push_back(Jets[i]); bJet1Idx = i;}		  
	    jetMatchindx.push_back(i);
	  }
	}
      }
    }

    // define ST - add jets
    for(int i=0;i<hadJets.size();i++){
      if( (abs(hadJets[i].Eta()) < 2.4) ){
	NHadJets++; 
	ST=ST+(hadJets[i].Pt());
      }
    }
    // define ST - add photon
    if( minDR<0.3) {
      ST=ST+bestPhoton.Pt();
      phoMatchingJetIndx = minDRindx;
    }
    
    // calculate dPhi jets & met
    if(NHadJets>=1)
      dPhi_METjet1 = abs(DeltaPhi(METPhi,hadJets[0].Phi()));
    if(NHadJets>=2)
      dPhi_METjet2 = abs(DeltaPhi(METPhi,hadJets[1].Phi()));
    

    Pass_Pho_pT  =false, Pass_MET100=false;
    Pass_NHadJets=false, Pass_ST=false;
    Pass_EvtCln  =false, Pass_JetMetPhi=false;
    Pass_EMu_veto=false, Pass_Iso_trk_veto=false;
    rmOverlap    =false;
    // Set event selection weights
    if(bestPhoton.Pt() >40) Pass_Pho_pT = true;
    if(MET > 100)           Pass_MET100 = true;
    if (NHadJets >=2)       Pass_NHadJets = true;		
    if(ST > 300)            Pass_ST     = true;
    ApplyTrgEff = true;
    if(PrimaryVertexFilter==1 && globalSuperTightHalo2016Filter==1 &&
       HBHENoiseFilter==1 &&HBHEIsoNoiseFilter==1 &&
       EcalDeadCellTriggerPrimitiveFilter == 1 && BadPFMuonFilter==1 &&
       BadPFMuonDzFilter==1 && eeBadScFilter==1 && ecalBadCalibFilter==1 &&
       NVtx>0 && PFCaloMETRatio < 5 &&
       (!(phoMatchingJetIndx>=0 && (Jets[phoMatchingJetIndx].Pt()/bestPhoton.Pt() < 1.0)) &&
	phoMatchingJetIndx >= 0) ) Pass_EvtCln = true;
    if(dPhi_METjet1 > 0.3 && dPhi_METjet2 > 0.3) Pass_JetMetPhi = true;
    if(RemoveSampleOverlap(s_sample, bestPhoton)) rmOverlap = true;
    if (NEMu == 0) Pass_EMu_veto = true;
    if (!(isoElectronTracks || isoMuonTracks || isoPionTracks)) Pass_Iso_trk_veto = true;

    bool selBaseline =  Pass_Pho_pT   && Pass_MET100 && 
      Pass_NHadJets && Pass_ST && 
      Pass_EvtCln   && Pass_JetMetPhi && 
      Pass_EMu_veto && Pass_Iso_trk_veto ;
    
    // call function for cutflow
    DoCutFlow(k, decade);
    
    //std::cout << jentry << " " << MET << std::endl;
    //double wt = 1.0; //0.12;
    h_MET->Fill(MET, wt);
    h_MET_test->Fill(100.0, wt);

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

  std::cout << "Final counts " << std::endl;
  std::cout << "h_MET->Integral()      " << h_MET->Integral()+h_MET->GetBinContent(101)      << std::endl;
  std::cout << "h_MET_test->Integral() " << h_MET_test->Integral() << std::endl;

  DoCutFlow(0, 0, true);
} // end EventLoop function



//===== get best photon ====
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

double AnalyzeTProxytBSM::getGenLep(myLV bestPhoton){
  //vector<TLorentzVector> v_genLep2;
  vector<myLV> v_genLep2;
  //TLorentzVector genMu1, genEle1;
  myLV genMu1, genEle1;
  // if(flag)
  //   {
  for(int i=0 ; i < GenElectrons->size(); i++)
    {
      if(GenElectrons[i].Pt()!=0)
	{
	  genEle1 = (GenElectrons[i]);
	  v_genLep2.push_back(genEle1);
	}
      
    }
  //   }
  // else
  //   {
  for(int i=0 ; i < GenMuons->size(); i++)
        {
          if(GenMuons[i].Pt()!=0)
            {
              genMu1 = (GenMuons[i]);
              v_genLep2.push_back(genMu1);
            }
        }
      //  }
      return MinDr_myLV(bestPhoton,v_genLep2);
}

// ====================
bool AnalyzeTProxytBSM::RemoveSampleOverlap(TString s_sample, myLV bestPhoton) {

  //script to define conditions to remove ovrlap
  
  bool cont1=true, cont2=true, cont3=true, cont4=true, cont5=true,
       cont6=true, cont7=true, cont8=true, cont9=true, cont10=true; 

  if((s_sample.Contains("TTJets_HT") || s_sample.Contains("TTJets-HT")) && madHT<600)
    cont1=false;

  if((s_sample.Contains("TTJets_inc")|| s_sample.Contains("TTJets_SingleLept") ||
      s_sample.Contains("TTJets_DiLept") || s_sample.Contains("TTJets_Leptons") ||
      s_sample.Contains("TTJets_Leptons")) && madHT>600)
    cont2=false;

  if(!genphocheck) {
    genphomatch_before++;
    double mindr_Pho_genlep=getGenLep(bestPhoton);
    
    if( s_sample.Contains("TTG") ) {
      if(!hasGenPromptPhoton) {
	//h_selectBaselineYields_v1->Fill("No gen prompt #gamma",wt);
	//if(jentry==0)cout<<"**********processing "<<s_sample<<" with non-prompt Gen photon"<<endl;
      } else if(hasGenPromptPhoton) {
	//h_selectBaselineYields_v1->Fill("Gen prompt #gamma",wt);
	if(!(madMinPhotonDeltaR >= 0.5 && mindr_Pho_genlep >=0.5 )) {
	 //h_phoPt_promptPho_rejected->Fill(bestPhoton.Pt(),wt);
	  //if(madMinPhotonDeltaR<0.5)h_selectBaselineYields_v1->Fill("madMinPhotonDeltaR <0.5",wt);
	  //if(mindr_Pho_genlep<0.5) h_selectBaselineYields_v1->Fill("mindr_Pho_genlep<0.5",wt);
	  cont3=false;
	} else {
	  //if(madMinPhotonDeltaR >= 0.5) h_selectBaselineYields_v1->Fill("mindR(q/g, #gamma)",wt);
	  //if(mindr_Pho_genlep >=0.5)    h_selectBaselineYields_v1->Fill("mindR(l, #gamma)",wt);
	}
      }
    }

    if(s_sample.Contains("WGJets_MonoPhoton_PtG-40to130UL") ||
       s_sample.Contains("WGJets_MonoPhoton_PtG-130UL")) {
      //if(s_sample.Contains("WGJets_MonoPhoton_PtG-40to130UL"||"WGJets_MonoPhoton_PtG-130UL"))
      if(!hasGenPromptPhoton) {
	//h_selectBaselineYields_v1->Fill("No gen prompt #gamma",wt);
	//if(jentry==0)cout<<"**********processing "<<s_sample<<" with non-prompt Gen photon"<<endl;
      } else if(hasGenPromptPhoton) {
	//h_selectBaselineYields_v1->Fill("Gen prompt #gamma",wt);
	if(!(madMinPhotonDeltaR >= 0.5 && mindr_Pho_genlep >=0.5 ))
	  {//h_phoPt_promptPho_rejected->Fill(bestPhoton.Pt(),wt);
	    //if(madMinPhotonDeltaR<0.5) h_selectBaselineYields_v1->Fill("madMinPhotonDeltaR <0.5",wt);
	    //if(mindr_Pho_genlep<0.5)   h_selectBaselineYields_v1->Fill("mindr_Pho_genlep<0.5",wt);
	    cont4=false;
	  } else {
	  //if(madMinPhotonDeltaR >= 0.5) h_selectBaselineYields_v1->Fill("mindR(q/g, #gamma)",wt);
	  //if(mindr_Pho_genlep >=0.5)    h_selectBaselineYields_v1->Fill("mindR(l, #gamma)",wt);
	}
      }
    }
    
    if(s_sample.Contains("WJets")) {
      if(!hasGenPromptPhoton) {
	//h_selectBaselineYields_v1->Fill("No gen prompt #gamma",wt);
	//if(jentry==0)cout<<"**********processing "<<s_sample<<" with non-prompt Gen photon"<<endl; 
      } else if(hasGenPromptPhoton) {
	//h_selectBaselineYields_v1->Fill("Gen prompt #gamma",wt);
	if(!(madMinPhotonDeltaR < 0.5 || mindr_Pho_genlep < 0.5)) {
	  //h_phoPt_promptPho_rejected->Fill(bestPhoton.Pt(),wt);
	  cont5=false;
	} else {
	  //if(madMinPhotonDeltaR >= 0.5) h_selectBaselineYields_v1->Fill("pass_mindR(q/g, #gamma)",wt); 
	  //if(mindr_Pho_genlep >=0.5)    h_selectBaselineYields_v1->Fill("pass_mindR(l, #gamma)",wt); 
	}
      }
    }
    
    if(s_sample.Contains("TTJets_HT") || s_sample.Contains("TTJets-HT")||
       s_sample.Contains("TTJets-inc")|| s_sample.Contains("TTJets_inc") ||
       s_sample.Contains("TTJets2_v17")||s_sample.Contains("TTJets")  ||
       s_sample.Contains("TTJets_Leptons")) {
      if(hasGenPromptPhoton) {	
	if(!(madMinPhotonDeltaR < 0.5 || mindr_Pho_genlep < 0.5)) {
	  cont6=false;
	}
      }
    }
	
    if(hasGenPromptPhoton && (s_sample.Contains("GJets"))) {
      if(!(madMinPhotonDeltaR>0.4)) cont7=false;
    }
	
    if(hasGenPromptPhoton && (s_sample.Contains("QCD"))) {
      if((madMinPhotonDeltaR>0.4 && hasGenPromptPhoton))
	cont8=false;
    }
	
    if(hasGenPromptPhoton && ((s_sample.Contains("ZG"))|| (s_sample.Contains("ZNuNuG"))
			      || s_sample.Contains("ZNuNuGJets"))) {
	    if(!(madMinPhotonDeltaR>0.5))
	      cont9=false;
    }
	
    if(hasGenPromptPhoton && ((s_sample.Contains("ZJets"))|| (s_sample.Contains("ZNuNuJets")))) {
      if(!(madMinPhotonDeltaR<=0.5))
	cont10=false;
    }
    genphomatch_after++;
  }

  bool rmOvrlp = cont1 && cont2 && cont3 && cont4 && cont5 && cont6 && cont7 && cont8 && cont9 && cont10;

  return rmOvrlp;
}

//SS//== not using this functions == 
Bool_t AnalyzeTProxytBSM::Process(Long64_t entry) {

  std::cout << entry << std::endl;
   fDirector.SetReadEntry(entry);
   std::cout<< "entry " << entry << " RunNum " << RunNum << std::endl;
   std::cout << "GenParticles->size() "<< GenParticles->size() << std::endl;
  return 0;
}

