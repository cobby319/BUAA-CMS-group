#ifndef objectselection_h
#define objectselection_h

#include <iostream>
#include <string>
#include <TMath.h>
#include <vector>
#include "TLorentzVectorWithIndex.h"
#include "Utils.h"
#include "RoccoR.h"

namespace objectSelection
{

  bool selectElectrons(std::vector<TLorentzVectorWithIndex> & selElectrons, std::vector<TLorentzVectorWithIndex> & extraElectrons, std::vector<float> *ElPt, std::vector<float> *ElEta, std::vector<float> *ElPhi, std::vector<float> *ElE, std::vector<unsigned int> *ElId, std::vector<float> *ElEtaSc);

  bool selectMuons(std::vector<TLorentzVectorWithIndex> & selMuons, std::vector<TLorentzVectorWithIndex> & extraMuons, std::vector<float> *MuPt, std::vector<float> *MuEta, std::vector<float> *MuPhi, std::vector<float> *MuE, std::vector<unsigned int> *MuId, std::vector<unsigned int> *MuIdTight, std::vector<unsigned int> *MuIdSoft, std::vector<float> *MuPfIso);

  bool selectPhotons(std::vector<TLorentzVectorWithIndex> & selPhotons, std::vector<float> *PhotPt, std::vector<float> *PhotEta, std::vector<float> *PhotPhi, std::vector<unsigned int> *PhotId, std::vector<float> *PhotScEta, std::vector<bool> *PhotHasPixelSeed, std::vector<float> *PhotSigmaIetaIeta, std::vector<TLorentzVectorWithIndex> & selMuons, std::vector<TLorentzVectorWithIndex> & selElectrons);

  bool selectPartiallyPhotons(std::vector<TLorentzVectorWithIndex> & selPhotons, std::vector<float> *PhotPt, std::vector<float> *PhotEta, std::vector<float> *PhotPhi, std::vector<float> *PhotScEta, std::vector<bool> *PhotHasPixelSeed, std::vector<TLorentzVectorWithIndex> & selMuons, std::vector<TLorentzVectorWithIndex> & selElectrons, std::vector<float> *PhotHoE, std::vector<float> *PhotSigmaIetaIeta, std::vector<float> *PhotPfIsoChHad, std::vector<float> *PhotPfIsoNeutralHad, std::vector<float> *PhotPfIsoPhot, Float_t EvtFastJetRho, TString selectionLevel);

  bool selectJets(std::vector<TLorentzVectorWithIndex> & selJets, std::vector<double> & btags, std::vector<float> *JetAk04Pt, std::vector<float> *JetAk04Eta, std::vector<float> *JetAk04Phi, std::vector<float> *JetAk04E, std::vector<float> *JetAk04Id, std::vector<float> *JetAk04NeutralEmFrac, std::vector<float> *JetAk04NeutralHadAndHfFrac, std::vector<float> *JetAk04NeutMult, std::vector<float> *JetAk04BDiscCisvV2, const std::vector<TLorentzVectorWithIndex> & selMuons, const std::vector<TLorentzVectorWithIndex> & selElectrons, const std::vector<TLorentzVectorWithIndex> & selPhotons);

  bool cleanPathologicEventsInPhotons(TString datasetName, float EvtRunNum, float EvtLumiNum, float EvtNum);


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
                   std::vector<float>   *J_vertexchi2);

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
                   std::vector<float>   *JPiPi_Prob);
}

#endif
