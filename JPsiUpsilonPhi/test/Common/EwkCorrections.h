#ifndef EwkCorrections_h
#define EwkCorrections_h

#include <vector>
#include <map>
#include <string>
#include <utility>
#include <iostream>
#include <fstream>
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "../Includes/SmartSelectionMonitor_hzz.h"

namespace EwkCorrections
{
  std::vector<std::vector<float>> readFile_and_loadEwkTable(TString catalogInputFile);
  std::vector<float> findCorrection(const std::vector<std::vector<float>> & Table_EWK, float sqrt_s_hat, float t_hat);
  std::map<std::string,std::pair<TLorentzVector,TLorentzVector>> reconstructGenLevelBosons(std::vector<float> *GLepBarePt, std::vector<float> *GLepBareEta, std::vector<float> *GLepBarePhi, std::vector<float> *GLepBareE, std::vector<int> *GLepBareId, std::vector<int> *GLepBareSt, std::vector<int> *GLepBareMomId);
  double getEwkCorrections(TString catalogInputFile, std::map<std::string,std::pair<TLorentzVector,TLorentzVector>> genLevelLeptons, const std::vector<std::vector<float>> & Table, double & ewkCorrections_error, std::vector<float> *GPdfx1, std::vector<float> *GPdfx2, std::vector<int> *GPdfId1, std::vector<int> *GPdfId2);
}

#endif
