#include "ObjectSelection.h"

namespace objectSelection
{

  bool selectElectrons(std::vector<TLorentzVectorWithIndex> & selElectrons, std::vector<TLorentzVectorWithIndex> & extraElectrons, std::vector<float> *ElPt, std::vector<float> *ElEta, std::vector<float> *ElPhi, std::vector<float> *ElE, std::vector<unsigned int> *ElId, std::vector<float> *ElEtaSc)
  {
    for(unsigned int i = 0 ; i<ElPt->size() ; i++){
      bool passEta = false, passId = false, passPt = false, passLoosePt = false, passLooseId = false;
      TLorentzVectorWithIndex currentLepton = TLorentzVectorWithIndex::PtEtaPhiEIndex(ElPt->at(i),ElEta->at(i),ElPhi->at(i),ElE->at(i), i);
      passId = ElId->at(i) & (1<<3);
      passLooseId = ElId->at(i) & (1<<1);
      double eta = fabs(ElEtaSc->at(i));//I took the supercluster eta since it's really the geometry which is taken here.
      passEta = (eta<=2.5 && (eta>=1.5660 || eta<=1.4442));
      passPt = (currentLepton.Pt() >=25);
      passLoosePt = (currentLepton.Pt() >=10);
      bool isLooseElectron = passEta && passLooseId && passLoosePt;
      bool isGoodElectron = passEta && passId && passPt; //Isolation is embedded in the ID
      if(isLooseElectron && !isGoodElectron) extraElectrons.push_back(currentLepton);
      if(isGoodElectron && selElectrons.size()==2) extraElectrons.push_back(currentLepton);
      if(isGoodElectron && selElectrons.size()<2) selElectrons.push_back(currentLepton);
    }
    return true;
  }

  bool selectMuons(std::vector<TLorentzVectorWithIndex> & selMuons, std::vector<TLorentzVectorWithIndex> & extraMuons, std::vector<float> *MuPt, std::vector<float> *MuEta, std::vector<float> *MuPhi, std::vector<float> *MuE, std::vector<unsigned int> *MuId, std::vector<unsigned int> *MuIdTight, std::vector<unsigned int> *MuIdSoft, std::vector<float> *MuPfIso)
  {
    //ID and ISO from this TWiki https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2
    for(unsigned int i = 0 ; i<MuPt->size() ; i++){
      bool passEta = false, passIso = false, passId = false, passPt = false, passLoosePt = false, passLooseIso = false, passLooseId = false, passSoftId = false, passSoftPt = false;
      TLorentzVectorWithIndex currentLepton = TLorentzVectorWithIndex::PtEtaPhiMIndex(MuPt->at(i),MuEta->at(i),MuPhi->at(i),0.1056, i);
      passId = MuIdTight->at(i) & (1<<0); //Look at the first vertex, hence the bit 0.
      passLooseId = MuId->at(i) & (1<<0);
      passSoftId = MuIdSoft->at(i) & (1<<0);
      double eta = fabs(MuEta->at(i));
      passEta = (eta<=2.4);
      //Iso //We use MuPfIso for now, we'll see after if it's mandatory to refine it.
      passIso = (MuPfIso->at(i)<0.15);
      passLooseIso = (MuPfIso->at(i)<0.25);
      passPt = (currentLepton.Pt() >=25);
      passLoosePt = (currentLepton.Pt() >=10);
      passSoftPt = (currentLepton.Pt() >=3);
      bool isLooseMuon = passEta && ( (passLooseId && passLoosePt && passLooseIso) || (passSoftId && passSoftPt) ); //Accounts for both loose or soft muons.
      bool isGoodMuon = passEta && passIso && passId && passPt;
      if(isLooseMuon && !isGoodMuon) extraMuons.push_back(currentLepton);
      if(isGoodMuon && selMuons.size()==2) extraMuons.push_back(currentLepton);
      if(isGoodMuon && selMuons.size()<2) selMuons.push_back(currentLepton);
    }
    return true;
  }

  bool selectPhotons(std::vector<TLorentzVectorWithIndex> & selPhotons, std::vector<float> *PhotPt, std::vector<float> *PhotEta, std::vector<float> *PhotPhi, std::vector<unsigned int> *PhotId, std::vector<float> *PhotScEta, std::vector<bool> *PhotHasPixelSeed, std::vector<float> *PhotSigmaIetaIeta, std::vector<TLorentzVectorWithIndex> & selMuons, std::vector<TLorentzVectorWithIndex> & selElectrons)
  {
    for(unsigned int i = 0 ; i<PhotPt->size() ; i++){
      bool passId = false, passPt = false, passEta = false, passLeptonCleaning = false, passSpikes = false;
      TLorentzVectorWithIndex currentPhoton = TLorentzVectorWithIndex::PtEtaPhiEIndex(PhotPt->at(i),PhotEta->at(i),PhotPhi->at(i),utils::getPhotonEnergy(PhotPt->at(i),PhotEta->at(i)), i); //photon energy is completely given by Pt and Eta.
      passId = PhotId->at(i) & (1<<2); //tight, according to llvv_fwk the code. FIXME: check that it's not better to redefine everything ourselves.
      passPt = (currentPhoton.Pt() >= 55);
      passEta = (fabs(PhotScEta->at(i))<=1.4442);
      passSpikes = PhotSigmaIetaIeta->at(i) > 0.001; //added to the ID atfer PhotonCR study. For now it is just on ieta but it should also be on iphi ==>FIXME
      double minDRlg(9999.); for(unsigned int ilep=0; ilep<selMuons.size(); ilep++) minDRlg = TMath::Min( minDRlg, utils::deltaR(currentPhoton,selMuons[ilep]) );
      for(unsigned int ilep=0; ilep<selElectrons.size(); ilep++) minDRlg = TMath::Min( minDRlg, utils::deltaR(currentPhoton,selElectrons[ilep]) );
      passLeptonCleaning = (minDRlg>=0.1); //according to the llvv_fwk code.
      if(passId && passPt && passEta && passLeptonCleaning && !PhotHasPixelSeed->at(i) && passSpikes) selPhotons.push_back(currentPhoton); //We ask for no pixel seed for the photons.
    }
    return true;
  }

  bool selectPartiallyPhotons(std::vector<TLorentzVectorWithIndex> & selPhotons, std::vector<float> *PhotPt, std::vector<float> *PhotEta, std::vector<float> *PhotPhi, std::vector<float> *PhotScEta, std::vector<bool> *PhotHasPixelSeed, std::vector<TLorentzVectorWithIndex> & selMuons, std::vector<TLorentzVectorWithIndex> & selElectrons, std::vector<float> *PhotHoE, std::vector<float> *PhotSigmaIetaIeta, std::vector<float> *PhotPfIsoChHad, std::vector<float> *PhotPfIsoNeutralHad, std::vector<float> *PhotPfIsoPhot, Float_t EvtFastJetRho, TString selectionLevel)
  {
    for(unsigned int i = 0 ; i<PhotPt->size() ; i++){
      bool passId = false, passPt = false, passEta = false, passLeptonCleaning = false;
      TLorentzVectorWithIndex currentPhoton = TLorentzVectorWithIndex::PtEtaPhiEIndex(PhotPt->at(i),PhotEta->at(i),PhotPhi->at(i),utils::getPhotonEnergy(PhotPt->at(i),PhotEta->at(i)), i); //photon energy is completely given by Pt and Eta.
      //Tight ID by hand with available variables for barrel photon
      double pt = currentPhoton.Pt();
      double sceta = PhotScEta->at(i);
      passId = true;
      if(selectionLevel.Contains("1") && PhotHoE->at(i) > 0.0269) passId = false; //HoE
      if(selectionLevel.Contains("2") && PhotSigmaIetaIeta->at(i) > 0.00994) passId = false;
      if(selectionLevel.Contains("3") && utils::photon_rhoCorrectedIso(PhotPfIsoChHad->at(i), EvtFastJetRho, sceta, "chIso") > 0.202) passId = false;
      if(selectionLevel.Contains("4") && utils::photon_rhoCorrectedIso(PhotPfIsoNeutralHad->at(i), EvtFastJetRho, sceta, "nhIso") > 0.264+0.0148*pt+0.000017*pt*pt) passId = false;
      if(selectionLevel.Contains("5") && utils::photon_rhoCorrectedIso(PhotPfIsoPhot->at(i), EvtFastJetRho, sceta, "gIso") > 2.362+0.0047*pt) passId = false;
      passPt = (pt >= 55);
      passEta = (fabs(sceta)<=1.4442);
      double minDRlg(9999.); for(unsigned int ilep=0; ilep<selMuons.size(); ilep++) minDRlg = TMath::Min( minDRlg, utils::deltaR(currentPhoton,selMuons[ilep]) );
      for(unsigned int ilep=0; ilep<selElectrons.size(); ilep++) minDRlg = TMath::Min( minDRlg, utils::deltaR(currentPhoton,selElectrons[ilep]) );
      passLeptonCleaning = (minDRlg>=0.1); //according to the llvv_fwk code.
      if(passId && passPt && passEta && passLeptonCleaning && !PhotHasPixelSeed->at(i)) selPhotons.push_back(currentPhoton); //We ask for no pixel seed for the photons.
    }
    return true;
  }

  bool selectJets(std::vector<TLorentzVectorWithIndex> & selJets, std::vector<double> & btags, std::vector<float> *JetAk04Pt, std::vector<float> *JetAk04Eta, std::vector<float> *JetAk04Phi, std::vector<float> *JetAk04E, std::vector<float> *JetAk04Id, std::vector<float> *JetAk04NeutralEmFrac, std::vector<float> *JetAk04NeutralHadAndHfFrac, std::vector<float> *JetAk04NeutMult, std::vector<float> *JetAk04BDiscCisvV2, const std::vector<TLorentzVectorWithIndex> & selMuons, const std::vector<TLorentzVectorWithIndex> & selElectrons, const std::vector<TLorentzVectorWithIndex> & selPhotons)
  {
    for(unsigned int i =0 ; i<JetAk04Pt->size() ; i++){
      bool passSelPt = false, passEta = false, passTightEta = false, passId = false, passLeptonCleaning = false, passPhotonCleaning = false;
      TLorentzVectorWithIndex currentJet = TLorentzVectorWithIndex::PtEtaPhiEIndex(JetAk04Pt->at(i),JetAk04Eta->at(i),JetAk04Phi->at(i),JetAk04E->at(i), i);
      passSelPt = (currentJet.Pt() >=30);
      double eta = fabs(currentJet.Eta());
      passEta = (eta <= 4.7);
      passTightEta = (eta <= 2.5);
      if(eta<2.7){
        passId = (JetAk04Id->at(i) >= 1); //This case is simple, it corresponds exactly to what we apply for now
        }
      float nef = JetAk04NeutralEmFrac->at(i);
      float nhf = JetAk04NeutralHadAndHfFrac->at(i);
      float nnp = JetAk04NeutMult->at(i);
      if(eta<3.0 && eta >=2.7){
        passId = (nef > 0.01 && nhf < 0.98 && nnp > 2);
        //passId = (JetAk04Id->at(i) >= 2); //Simpler criterium, but not equivalent to what is mentionned in the AN. Was applied before having access to nnp.
        }
      if(eta>=3.0){
        passId = (nef<0.90 && nnp > 10);
        //passId = (JetAk04Id->at(i) >= 3); //Simpler criterium, but not equivalent to what is mentionned in the AN. Was applied before having access to nnp.
        }
      double minDRlj(9999.); for(unsigned int ilep=0; ilep<selMuons.size(); ilep++) minDRlj = TMath::Min( minDRlj, utils::deltaR(currentJet,selMuons[ilep]) );
      for(unsigned int ilep=0; ilep<selElectrons.size(); ilep++) minDRlj = TMath::Min( minDRlj, utils::deltaR(currentJet,selElectrons[ilep]) );
      passLeptonCleaning = (minDRlj>=0.4);
      double minDRgj(9999.); for(unsigned int ilep=0; ilep<selPhotons.size(); ilep++) minDRgj = TMath::Min( minDRgj, utils::deltaR(currentJet,selPhotons[ilep]) );
      passPhotonCleaning = (minDRgj>=0.4);
      if(passSelPt && passEta && passId && passLeptonCleaning && passPhotonCleaning) selJets.push_back(currentJet);
      if(passSelPt && passTightEta && passId && passLeptonCleaning && passPhotonCleaning) btags.push_back(JetAk04BDiscCisvV2->at(i));
      }
    return true;
  }

  bool cleanPathologicEventsInPhotons(TString datasetName, float EvtRunNum, float EvtLumiNum, float EvtNum){
    bool eventShouldBeCleaned = false;
    //We remove events in the spike for those samples 
    if(datasetName.Contains("QCD_HT100to200")){
      //Spike at 190~200GeV (cut if pt > 190GeV)
      if( EvtRunNum == 1 && EvtLumiNum == 21997  && EvtNum == 32986438 ) eventShouldBeCleaned = true;
      if( EvtRunNum == 1 && EvtLumiNum == 123682 && EvtNum == 185472705) eventShouldBeCleaned = true;
      if( EvtRunNum == 1 && EvtLumiNum == 133696 && EvtNum == 200489234) eventShouldBeCleaned = true;
      if( EvtRunNum == 1 && EvtLumiNum == 301998 && EvtNum == 452875030) eventShouldBeCleaned = true;
      if( EvtRunNum == 1 && EvtLumiNum == 237717 && EvtNum == 356480214) eventShouldBeCleaned = true;
      
      //Spike at ~330-340GeV (cut if pt > 330GeV)
      if( EvtRunNum == 1 && EvtLumiNum == 405615 && EvtNum == 608258936) eventShouldBeCleaned = true;
      if( EvtRunNum == 1 && EvtLumiNum == 238040 && EvtNum == 356963627) eventShouldBeCleaned = true;
      if( EvtRunNum == 1 && EvtLumiNum == 192575 && EvtNum == 288784917) eventShouldBeCleaned = true;
      if( EvtRunNum == 1 && EvtLumiNum == 405110 && EvtNum == 607502440) eventShouldBeCleaned = true;
      if( EvtRunNum == 1 && EvtLumiNum == 398170 && EvtNum == 597094584) eventShouldBeCleaned = true;
      if( EvtRunNum == 1 && EvtLumiNum == 242217 && EvtNum == 363227739) eventShouldBeCleaned = true;
      if( EvtRunNum == 1 && EvtLumiNum == 175934 && EvtNum == 263829468) eventShouldBeCleaned = true;
      if( EvtRunNum == 1 && EvtLumiNum == 336765 && EvtNum == 505011533) eventShouldBeCleaned = true;
      if( EvtRunNum == 1 && EvtLumiNum == 252239 && EvtNum == 378257731) eventShouldBeCleaned = true;
      if( EvtRunNum == 1 && EvtLumiNum == 239877 && EvtNum == 359718643) eventShouldBeCleaned = true;
    }
    if(datasetName.Contains("QCD_Pt-30to50_EMEnriched")){
      //Spike at ~330-340GeV (cut if pt > 330GeV)
      if( EvtRunNum == 1 && EvtLumiNum == 5960 && EvtNum == 16066353) eventShouldBeCleaned = true;
    }
    if(datasetName.Contains("QCD_Pt-20toInf_MuEnrichedPt15")){
      //Spike at ~330-340GeV (cut if pt > 330GeV)
      if( EvtRunNum == 1 && EvtLumiNum == 181694 && EvtNum == 2097414398) eventShouldBeCleaned = true;
      if( EvtRunNum == 1 && EvtLumiNum == 40564  && EvtNum == 3384544677) eventShouldBeCleaned = true;
      if( EvtRunNum == 1 && EvtLumiNum == 49551  && EvtNum == 1742764771) eventShouldBeCleaned = true;
      if( EvtRunNum == 1 && EvtLumiNum == 94563  && EvtNum == 2145642832) eventShouldBeCleaned = true;
    }
    return eventShouldBeCleaned;
  }
    bool selectJpsi (std::vector<TLorentzVectorWithIndex> & selJpsi,
                   std::vector<float>  *J_mass,
                   std::vector<float>  *J_px,
                   std::vector<float>  *J_py,
                   std::vector<float>  *J_pz,
                   std::vector<float>  *J_px1,
                   std::vector<float>  *J_py1,
                   std::vector<float>  *J_pz1,
                   std::vector<int>    *J_charge1,
                   std::vector<float>  *J_px2,
                   std::vector<float>  *J_py2,
                   std::vector<float>  *J_pz2,
                   std::vector<int>    *J_charge2,
                   std::vector<float>  *mumC2,
                   std::vector<int>    *mumNHits,
                   std::vector<int>    *mumNPHits,
                   std::vector<float>  *mupC2,
                   std::vector<int>    *mupNHits,
                   std::vector<int>    *mupNPHits,
                   std::vector<float>  *mumdxy,
                   std::vector<float>  *mupdxy,
                   std::vector<float>  *mumdz,
                   std::vector<float>  *mupdz,
                   std::vector<float>  *muon_dca,
                   std::vector<bool>   *mu1soft,
                   std::vector<bool>   *mu2soft,
                   std::vector<bool>   *mu1tight,
                   std::vector<bool>   *mu2tight,
                   std::vector<bool>   *mu1PF,
                   std::vector<bool>   *mu2PF,
                   std::vector<bool>   *mu1loose,
                   std::vector<bool>   *mu2loose,
                   std::vector<float>   *J_lxy,
                   std::vector<float>   *J_lxyErr,
                   std::vector<float>   *J_vertexchi2,
                   std::vector<float>   *J_Prob){
      for(unsigned int i = 0 ; i<J_mass->size() ; i++){
        bool passMuonLooseID   = false;
        bool passJpsiPt        = false;
        bool passMuonHits      = false;
        bool passMuonPixelHits = false;
        TLorentzVector Jpsi ;
        Jpsi.SetXYZM(J_px->at(i),J_py->at(i),J_pz->at(i),J_mass->at(i));
        passMuonLooseID = mu1loose->at(i) && mu2loose->at(i);
        passJpsiPt = Jpsi.Pt() >10 ;
        passMuonHits = mupNHits->at(i) >6 && mumNHits->at(i) >6;
        passMuonPixelHits = mupNPHits->at(i) >0 &&  mumNPHits->at(i) >0;
        TLorentzVectorWithIndex JpsiWithIndex = TLorentzVectorWithIndex(Jpsi,i);
        if (passMuonLooseID && passJpsiPt && passMuonHits && passMuonPixelHits) {selJpsi.push_back(JpsiWithIndex);}
        if (selJpsi.size() >1 ) {
          auto maxprob = std::max_element(J_Prob->begin(), J_Prob->begin()+i);
          int jpsiN =std::distance(J_Prob->begin(), maxprob) ;
          selJpsi.clear();
          TLorentzVector Jpsi_tmp;
          Jpsi_tmp.SetXYZM(J_px->at(jpsiN),J_py->at(jpsiN),J_pz->at(jpsiN),J_mass->at(jpsiN));
          TLorentzVectorWithIndex Jpsi_tmpWithIndex = TLorentzVectorWithIndex(Jpsi_tmp,jpsiN);
          selJpsi.push_back(Jpsi_tmpWithIndex);
          //std::cout <<"push_back extrajpsi"<<std::endl;
        }
      }
      return selJpsi.size() > 0;
    }

  bool selectPions(std::vector<TLorentzVectorWithIndex> & selPion1,
                   std::vector<TLorentzVectorWithIndex> & selPion2,
                   std::vector<TLorentzVectorWithIndex> & selJpsi,
                   std::vector<float>   *JPiPi_lxy,
                   std::vector<float>   *JPiPi_lxyErr,
                   std::vector<float>   *JPiPi_x,
                   std::vector<float>   *JPiPi_y,
                   std::vector<float>   *JPiPi_z,
                   std::vector<int>     *Pi_nhits1,
                   std::vector<int>     *Pi_npixelhits1,
                   std::vector<int>     *Pi_nhits2,
                   std::vector<int>     *Pi_npixelhits2,
                   std::vector<float>   *Pi_eta1,
                   std::vector<float>   *Pi_eta2,
                   std::vector<float>   *Pi_phi1,
                   std::vector<float>   *Pi_phi2,
                   std::vector<float>   *Pi_pt1,
                   std::vector<float>   *Pi_pt2,
                   std::vector<float>   *Pi_e1,
                   std::vector<float>   *Pi_e2,
                   std::vector<int>     *Pi_charge1,
                   std::vector<int>     *Pi_charge2,
                   std::vector<float>   *Pi_dxy1,
                   std::vector<float>   *Pi_dxy2,
                   std::vector<float>   *Pi_dxyerr1,
                   std::vector<float>   *Pi_dxyerr2,
                   std::vector<float>   *Pi_vertexchisq1,
                   std::vector<float>   *Pi_vertexchisq2,
                   std::vector<float>   *Pi1_hcalFraction,
                   std::vector<float>   *Pi2_hcalFraction,
                   std::vector<float>   *Pi1_vertexNdof,
                   std::vector<float>   *Pi2_vertexNdof,
                   std::vector<float>   *Pi1_vertexNchi2,
                   std::vector<float>   *Pi2_vertexNchi2,
                   std::vector<float>   *Pi1_lambda,
                   std::vector<float>   *Pi2_lambda,
                   std::vector<float>   *Pi1_lambdaError,
                   std::vector<float>   *Pi2_lambdaError,
                   std::vector<float>   *Pi1_qoverp,
                   std::vector<float>   *Pi2_qoverp,
                   std::vector<float>   *Pi1_qoverpError,
                   std::vector<float>   *Pi2_qoverpError,
                   std::vector<float>   *Pi1_validTkFraction,
                   std::vector<float>   *Pi2_validTkFraction,
                   std::vector<int>     *Pi1_numberOfMothers,
                   std::vector<int>     *Pi2_numberOfMothers,
                   std::vector<int>     *Pi1_pdgId,
                   std::vector<int>     *Pi2_pdgId,
                   std::vector<int>     *Pi1_numberOfValidHitsOnTrack,
                   std::vector<int>     *Pi2_numberOfValidHitsOnTrack,
                   std::vector<bool>    *Pi1_isCaloMuon,
                   std::vector<bool>    *Pi2_isCaloMuon,
                   std::vector<bool>    *Pi1_isConvertedPhoton,
                   std::vector<bool>    *Pi2_isConvertedPhoton,
                   std::vector<bool>    *Pi1_isElectron,
                   std::vector<bool>    *Pi2_isElectron,
                   std::vector<bool>    *Pi1_isMuon,
                   std::vector<bool>    *Pi2_isMuon,
                   std::vector<bool>    *Pi1_isPhoton,
                   std::vector<bool>    *Pi2_isPhoton ,
                   std::vector<bool>    *Pi1_isGlobalMuon,
                   std::vector<bool>    *Pi2_isGlobalMuon,
                   std::vector<bool>    *Pi1_isJet,
                   std::vector<bool>    *Pi2_isJet,
                   std::vector<bool>    *Pi1_isLonglived,
                   std::vector<bool>    *Pi2_isLonglived,
                   std::vector<bool>    *Pi1_massConstraint,
                   std::vector<bool>    *Pi2_massConstraint,
                   std::vector<float>   *J_Prob,
                   std::vector<float>   *JPi_Prob,
                   std::vector<float>   *JPiPi_Prob){
    for(unsigned int i = 0 ; i<Pi_pt1->size() ; i++){
      bool passDeltaRJP1 = false;
      bool passDeltaRJP2 = false;
      bool passDeltaRPP  = false;
      bool pipiInPhiMass = false;
      bool pipiInXMass = false;
      bool pipiInLowMass =false;
      bool passPiPimassregion = false;
      bool passVertexNormalizedChi2 = false;
      bool passEta =false;
      bool passPt  = false;
      bool passIsNotOtherObject = false;
      TLorentzVector pion1,pion2;
      float Pion_mass = 0.13957061;
      pion1.SetPtEtaPhiE(Pi_pt1->at(i),Pi_eta1->at(i),Pi_phi1->at(i),Pi_e1->at(i));
      pion2.SetPtEtaPhiE(Pi_pt2->at(i),Pi_eta2->at(i),Pi_phi2->at(i),Pi_e2->at(i));
      pipiInPhiMass = (pion1+pion2).M() >1.01  && (pion1+pion2).M() <1.03 ;
      pipiInXMass =  (pion1+pion2).M() >0.81  && (pion1+pion2).M() <0.97 ;
      pipiInLowMass = (pion1+pion2).M() <0.35;
      passPiPimassregion = !pipiInXMass && !pipiInLowMass && !pipiInPhiMass;
      pion1.SetPtEtaPhiM(Pi_pt1->at(i),Pi_eta1->at(i),Pi_phi1->at(i),Pion_mass);
      pion2.SetPtEtaPhiM(Pi_pt2->at(i),Pi_eta2->at(i),Pi_phi2->at(i),Pion_mass);
      passDeltaRJP1 = selJpsi.at(0).DeltaR(pion1) < 1.2;
      passDeltaRJP2 = selJpsi.at(0).DeltaR(pion2) < 1.2;
      passDeltaRPP  = pion1.DeltaR(pion2) <1.8;
      passVertexNormalizedChi2 = Pi1_vertexNchi2->at(i) < 10.0 && Pi2_vertexNchi2->at(i) <10.0;
      passEta = Pi_eta1->at(i) <2.0 && Pi_eta1->at(i) >-2.0 && Pi_eta2->at(i) <2.0 &&Pi_eta2->at(i) >-2.0;
      passPt =  Pi_pt1->at(i) > 0.8 && Pi_pt2->at(i) >0.8;
      passIsNotOtherObject = true;//! (Pi1_isGlobalMuon || Pi2_isGlobalMuon);
      TLorentzVectorWithIndex pion1WithIndex = TLorentzVectorWithIndex(pion1,i);
      TLorentzVectorWithIndex pion2WithIndex = TLorentzVectorWithIndex(pion2,i);

      if(passDeltaRJP1 && passDeltaRJP2 && passDeltaRPP &&passPiPimassregion && passVertexNormalizedChi2 && passEta && passPt &&passIsNotOtherObject) {
        //std::cout <<"push_back pions"<<std::endl;
        selPion1.push_back(pion1WithIndex);
        selPion2.push_back(pion2WithIndex);
      }
      if(selPion1.size()>1){
          auto maxprob = std::max_element(JPiPi_Prob->begin(), JPiPi_Prob->begin() + i);
          int piN =std::distance(JPiPi_Prob->begin(), maxprob) ;
          //std::cout <<"a little test"<<std::endl;
          selPion1.clear();
          selPion2.clear();
          //TLorentzVector pion1_tmp;
          //TLorentzVector pion2_tmp;
          //pion1_tmp.SetPtEtaPhiM(Pi_pt1->at(piN),Pi_eta1->at(piN),Pi_phi1->at(piN),Pion_mass);
          //pion2_tmp.SetPtEtaPhiM(Pi_pt2->at(piN),Pi_eta2->at(piN),Pi_phi2->at(piN),Pion_mass);
          //TLorentzVectorWithIndex pion1_tmpWithIndex = ;
         // TLorentzVectorWithIndex pion2_tmpWithIndex = ;
          //std::cout <<"push_back extra pions"<<std::endl;
          selPion1.push_back(TLorentzVectorWithIndex::PtEtaPhiMIndex(Pi_pt1->at(piN),Pi_eta1->at(piN),Pi_phi1->at(piN),Pion_mass, piN));
          selPion2.push_back(TLorentzVectorWithIndex::PtEtaPhiMIndex(Pi_pt2->at(piN),Pi_eta2->at(piN),Pi_phi2->at(piN),Pion_mass, piN));
      }
    }
    return selPion1.size()>0 ;
  }

}
