#ifndef objectselection_h
#define objectselection_h

#include <iostream>
#include <string>
#include <TMath.h>
#include <vector>
#include "TLorentzVectorWithIndex.h"
#include "Utils.h"
#include "../Common/RoccoR.h"

namespace objectSelection
{

  bool selectElectrons(std::vector<TLorentzVectorWithIndex> & selElectrons, std::vector<TLorentzVectorWithIndex> & extraElectrons, std::vector<float> *ElPt, std::vector<float> *ElEta, std::vector<float> *ElPhi, std::vector<float> *ElE, std::vector<unsigned int> *ElId, std::vector<float> *ElEtaSc);

  bool selectMuons(std::vector<TLorentzVectorWithIndex> & selMuons, std::vector<TLorentzVectorWithIndex> & extraMuons, std::vector<float> *MuPt, std::vector<float> *MuEta, std::vector<float> *MuPhi, std::vector<float> *MuE, std::vector<unsigned int> *MuId, std::vector<unsigned int> *MuIdTight, std::vector<unsigned int> *MuIdSoft, std::vector<float> *MuPfIso);

  bool selectPhotons(std::vector<TLorentzVectorWithIndex> & selPhotons, std::vector<float> *PhotPt, std::vector<float> *PhotEta, std::vector<float> *PhotPhi, std::vector<unsigned int> *PhotId, std::vector<float> *PhotScEta, std::vector<bool> *PhotHasPixelSeed, std::vector<float> *PhotSigmaIetaIeta, std::vector<TLorentzVectorWithIndex> & selMuons, std::vector<TLorentzVectorWithIndex> & selElectrons);

  bool selectPartiallyPhotons(std::vector<TLorentzVectorWithIndex> & selPhotons, std::vector<float> *PhotPt, std::vector<float> *PhotEta, std::vector<float> *PhotPhi, std::vector<float> *PhotScEta, std::vector<bool> *PhotHasPixelSeed, std::vector<TLorentzVectorWithIndex> & selMuons, std::vector<TLorentzVectorWithIndex> & selElectrons, std::vector<float> *PhotHoE, std::vector<float> *PhotSigmaIetaIeta, std::vector<float> *PhotPfIsoChHad, std::vector<float> *PhotPfIsoNeutralHad, std::vector<float> *PhotPfIsoPhot, Float_t EvtFastJetRho, TString selectionLevel);

  bool selectJets(std::vector<TLorentzVectorWithIndex> & selJets, std::vector<double> & btags, std::vector<float> *JetAk04Pt, std::vector<float> *JetAk04Eta, std::vector<float> *JetAk04Phi, std::vector<float> *JetAk04E, std::vector<float> *JetAk04Id, std::vector<float> *JetAk04NeutralEmFrac, std::vector<float> *JetAk04NeutralHadAndHfFrac, std::vector<float> *JetAk04NeutMult, std::vector<float> *JetAk04BDiscCisvV2, const std::vector<TLorentzVectorWithIndex> & selMuons, const std::vector<TLorentzVectorWithIndex> & selElectrons, const std::vector<TLorentzVectorWithIndex> & selPhotons);

  bool cleanPathologicEventsInPhotons(TString datasetName, float EvtRunNum, float EvtLumiNum, float EvtNum);
}

#endif
