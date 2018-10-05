#include "ObjectSelection.h"

namespace objectSelection
{
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
        bool passMuonPt        = false;

      ///////////////////////////////////
      /* define variables for selection*/
      ///////////////////////////////////
        TLorentzVector Jpsi ;
        TLorentzVector muon1,muon2;
        float muon_mass = 0.10565837;
        muon1.SetXYZM(J_px1->at(i),J_py1->at(i),J_pz1->at(i),muon_mass);
        muon2.SetXYZM(J_px2->at(i),J_py2->at(i),J_pz2->at(i),muon_mass);
        Jpsi.SetXYZM(J_px->at(i),J_py->at(i),J_pz->at(i),J_mass->at(i));

      ///////////////////////////////////
      /* initiating the bool variables */
      ///////////////////////////////////
        passMuonLooseID = mu1loose->at(i) && mu2loose->at(i);
        passMuonPt = muon1.Pt() >4 && muon2.Pt() >4;
        passJpsiPt = Jpsi.Pt() >5 ;
        passMuonHits = mupNHits->at(i) >5 && mumNHits->at(i) >5;
        passMuonPixelHits = mupNPHits->at(i) >0 &&  mumNPHits->at(i) >0;
        TLorentzVectorWithIndex JpsiWithIndex = TLorentzVectorWithIndex(Jpsi,i);
        if (passMuonPt && passMuonLooseID && passJpsiPt && passMuonHits && passMuonPixelHits) {selJpsi.push_back(JpsiWithIndex);}
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
                   std::vector<float>   *JPiPi_Prob,
                   std::vector<float>   *JPiPi_px,
                   std::vector<float>   *JPiPi_py,
                   std::vector<float>   *JPiPi_pz){
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
      bool passLambda = false;
      bool passLxy =false;
      bool passProb = false;
      bool passPdgId = false;
      bool passCosine = false;
      bool passPiPidxy = false;
      bool passNhits = false;
      ///////////////////////////////////
      /* define variables for selection*/
      ///////////////////////////////////
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
      float px = JPiPi_px->at(i);
      float py = JPiPi_py->at(i);
      float rx = JPiPi_x->at(i);
      float ry = JPiPi_y->at(i);  
      float cosine = (px*rx+py*ry)/((px*px+py*py)*(rx*rx+ry*ry));

      ///////////////////////////////////
      /* initiating the bool variables */
      ///////////////////////////////////
      passCosine = cosine >0.95 ;
      passDeltaRJP1 = selJpsi.at(0).DeltaR(pion1) < 1.2;
      passDeltaRJP2 = selJpsi.at(0).DeltaR(pion2) < 1.2;
      passDeltaRPP  = pion1.DeltaR(pion2) <1.8;
      passVertexNormalizedChi2 = true; // Pi1_vertexNchi2->at(i) < 10.0 && Pi2_vertexNchi2->at(i) <10.0;
      passEta = Pi_eta1->at(i) <2.0 && Pi_eta1->at(i) >-2.0 && Pi_eta2->at(i) <2.0 &&Pi_eta2->at(i) >-2.0;
      passPt =  Pi_pt1->at(i) > 0.8 && Pi_pt2->at(i) >0.8;
      passIsNotOtherObject =  !(Pi1_isGlobalMuon->at(i) ||  Pi2_isGlobalMuon->at(i) );
      passLambda = Pi1_lambda->at(i) <1.0 && Pi1_lambda->at(i) >-1.0 &&Pi2_lambda->at(i) <1.0 &&Pi2_lambda->at(i) >-1.0;
      passLxy = JPiPi_lxy->at(i) >0.025;
      passProb = JPiPi_Prob->at(i) >0.1 ;
      passPdgId = (Pi1_pdgId->at(i) == 211 || Pi1_pdgId->at(i) == -211) && (Pi2_pdgId->at(i) == 211 || Pi2_pdgId->at(i) == -211);
      passPiPidxy = Pi_dxy1->at(i)/Pi_dxyerr1->at(i) >3.0 && Pi_dxy2->at(i)/pi_dxy2->at(i) >2.0;
      passNhits = Pi_nhits1->at(i) >5 && Pi_nhits2->at(i) >5 ;
      /////////////////////////////////////
      /*push back if pions pass all cuts */
      /////////////////////////////////////
      TLorentzVectorWithIndex pion1WithIndex = TLorentzVectorWithIndex(pion1,i);
      TLorentzVectorWithIndex pion2WithIndex = TLorentzVectorWithIndex(pion2,i);
      if( passPiPidxy&&
          passNhits &&
          passCosine && 
          passPdgId && 
          passProb && 
          passLxy && 
          passDeltaRJP1 && 
          passDeltaRJP2 && 
          passDeltaRPP &&
          passPiPimassregion && 
          passVertexNormalizedChi2 && 
          passEta && 
          passPt &&
          passIsNotOtherObject && 
          passLambda) {

        selPion1.push_back(pion1WithIndex);
        selPion2.push_back(pion2WithIndex);
      }
      if(selPion1.size()>1){
          auto maxprob = std::max_element(JPiPi_Prob->begin(), JPiPi_Prob->begin() + i);
          int piN =std::distance(JPiPi_Prob->begin(), maxprob) ;
          selPion1.clear();
          selPion2.clear();
          selPion1.push_back(TLorentzVectorWithIndex::PtEtaPhiMIndex(Pi_pt1->at(piN),Pi_eta1->at(piN),Pi_phi1->at(piN),Pion_mass, piN));
          selPion2.push_back(TLorentzVectorWithIndex::PtEtaPhiMIndex(Pi_pt2->at(piN),Pi_eta2->at(piN),Pi_phi2->at(piN),Pion_mass, piN));
      }
    }
    return selPion1.size()>0 ;
  }

}
