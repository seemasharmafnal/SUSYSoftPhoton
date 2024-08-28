#ifdef ANALYZETPROXYTBSM_cxx

void AnalyzeTProxytBSM::DoCutFlow(int k, int decade, bool printsummary, bool debug) {

    if (k > decade && debug){
      std::cout << "In function DoCutFlow " << std::endl;
      //std::cout << "h_MET->Integral()      " << h_MET->Integral()      << std::endl;
      //std::cout << "h_MET_test->Integral() " << h_MET_test->Integral() << std::endl;
      std::cout << "NEvtSel[selNone]       " <<  NEvtSel[selNone]      << std::endl;
      std::cout << "NEvtSel[selPho]        " <<  NEvtSel[selPho]      << std::endl;
      std::cout << "NEvtSel[selMET]        " <<  NEvtSel[selMET]      << std::endl;
      std::cout << "NEvtSel[selNJet]       " <<  NEvtSel[selNJet]      << std::endl;
      std::cout << "NEvtSel[selST]         " <<  NEvtSel[selST]      << std::endl;
      std::cout << "NEvtSel[selTrgEff]     " <<  NEvtSel[selTrgEff]      << std::endl;
      std::cout << "NEvtSel[selEvtCln]     " <<  NEvtSel[selEvtCln]      << std::endl;
      std::cout << "NEvtSel[selJetMetPhi]  " <<  NEvtSel[selJetMetPhi]      << std::endl;
      std::cout << "NEvtSel[selRmOverlap]  " <<  NEvtSel[selRmOverlap]      << std::endl;
      std::cout << "NEvtSel[selEMu]        " <<  NEvtSel[selEMu]      << std::endl;
      std::cout << "NEvtSel[selIsoTrk]     " <<  NEvtSel[selIsoTrk]      << std::endl;
    }
  
  double wt = EvtWeight;
  // enum evtSel {selNone=0, selPho=1, selMET=2, selNJet=3, selST=4};
    NEvtSel[selNone] = NEvtSel[selNone] + wt;
    if(Pass_Pho_pT) {
      NEvtSel[selPho] = NEvtSel[selPho] + wt;
      if(Pass_MET100) {
	NEvtSel[selMET] = NEvtSel[selMET] + wt;
	if(Pass_NHadJets) {
	  NEvtSel[selNJet] = NEvtSel[selNJet] + wt;
	  if(Pass_ST) {
	    NEvtSel[selST] = NEvtSel[selST] + wt;
	    if (ApplyTrgEff) wt = wt * (((TMath::Erf((MET - p0)/p1)+1)/2.0)*p2);
	    NEvtSel[selTrgEff] = NEvtSel[selTrgEff] + wt;
	    if(Pass_EvtCln){
	      NEvtSel[selEvtCln] = NEvtSel[selEvtCln] + wt;
	      if(Pass_JetMetPhi){
		NEvtSel[selJetMetPhi] = NEvtSel[selJetMetPhi] + wt;
		if(rmOverlap) NEvtSel[selRmOverlap] = NEvtSel[selRmOverlap] + wt;
		if(Pass_EMu_veto){
		  NEvtSel[selEMu] = NEvtSel[selEMu] + wt;
		  if(Pass_Iso_trk_veto)
		    NEvtSel[selIsoTrk] = NEvtSel[selIsoTrk] + wt;
		}
	      }
	    }
	  }
	}
      }
    }

    if(printsummary) {
      std::cout << "NEvtSel[selNone]       " <<  NEvtSel[selNone]      << std::endl;
      std::cout << "NEvtSel[selPho]        " <<  NEvtSel[selPho]      << std::endl;
      std::cout << "NEvtSel[selMET]        " <<  NEvtSel[selMET]      << std::endl;
      std::cout << "NEvtSel[selNJet]       " <<  NEvtSel[selNJet]      << std::endl;
      std::cout << "NEvtSel[selST]         " <<  NEvtSel[selST]      << std::endl;
      std::cout << "NEvtSel[selTrgEff]     " <<  NEvtSel[selTrgEff]      << std::endl;
      std::cout << "NEvtSel[selEvtCln]     " <<  NEvtSel[selEvtCln]      << std::endl;
      std::cout << "NEvtSel[selJetMetPhi]  " <<  NEvtSel[selJetMetPhi]      << std::endl;
      std::cout << "NEvtSel[selRmOverlap]  " <<  NEvtSel[selRmOverlap]      << std::endl;
      std::cout << "NEvtSel[selEMu]        " <<  NEvtSel[selEMu]      << std::endl;
      std::cout << "NEvtSel[selIsoTrk]     " <<  NEvtSel[selIsoTrk]      << std::endl;
    }
}

#endif // AnalyzeTProxytBSM_cxx
