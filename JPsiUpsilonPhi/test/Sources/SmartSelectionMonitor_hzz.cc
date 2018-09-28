#include "../Includes/SmartSelectionMonitor_hzz.h"

bool SmartSelectionMonitor_hzz::declareHistos(){ //FIXME: Later, will take an array as input for the binnings.
  return true;
}

bool SmartSelectionMonitor_hzz::declareHistos_InstrMET(){ 
  return true;  
} 

bool SmartSelectionMonitor_hzz::declareHistos_NRB(){ 
  return true;  
}
bool SmartSelectionMonitor_hzz::declareHistos_jpsipipi(){
 
 addHistogram( new TH1F("Num_JPsi",";N_{J/#psi};Events",5,0,5));
 addHistogram( new TH1F("M_JPsi",";M_{J/#psi};Events",100,2.8,3.3));
 addHistogram( new TH1F("pT_JPsi",";p_{T,J/#Psi};Events",200,0,100)); 
 addHistogram( new TH1F("M_JPsiPiPi3-8",";M_{J/#psi#pi+#pi-};Events",100,3,8));
 addHistogram( new TH1F("M_JPsiPiPi4-5",";M_{J/#psi#pi+#pi-};Events",20,4,5));
 addHistogram( new TH1F("cosine of P&r",";cos#alpha;Events",20,0,1)); 
 addHistogram( new TH1F("Deta_JpsiPi1",";#delta #eta;Events",80,-4,4)); 
 addHistogram( new TH1F("Deta_JpsiPi2",";#delta #eta;Events",80,-4,4)); 
 addHistogram( new TH1F("Dphi_JpsiPi1",";#delta #eta;Events",80,-4,4)); 
 addHistogram( new TH1F("Dphi_JpsiPi2",";#delta #eta;Events",80,-4,4)); 
 addHistogram( new TH1F("DeltaR_JpsiPi1",";#delta R;Events",40,0,4)); 
 addHistogram( new TH1F("DeltaR_JpsiPi2",";#delta R;Events",40,0,4)); 
 addHistogram( new TH1F("lxy_jpsipipi",";L_{xy};Events",100,0,0.1)); 
 addHistogram( new TH1F("M_JPsiPicut4.1-4.2",";m_{J/psi,Pi};Events",40,3.5,4.3));
 addHistogram( new TH1F("M_JPsiPicut4.2-4.25",";m_{J/psi,Pi};Events",40,3.5,4.3)); 
 addHistogram( new TH1F("M_JPsiPicut4.25-4.3",";m_{J/psi,Pi};Events",40,3.5,4.3)); 
 addHistogram( new TH1F("M_JPsiPicut4.3-4.4",";m_{J/psi,Pi};Events",40,3.5,4.3)); 
 addHistogram( new TH1F("M_JPsiPicut4.4-4.7",";m_{J/psi,Pi};Events",40,3.5,4.3)); 
 addHistogram( new TH1F("M_JPsiPicut4.7-5.0",";m_{J/psi,Pi};Events",40,3.5,4.3)); 
 addHistogram( new TH2F("M_JpsiPi1&M_JpsiPi2",";m^{2}_{J/psi,Pi1};m^{2}_{J/psi,Pi2}",100,0,50,100,0,50));
 addHistogram( new TH1F("J_mass",";;Events",80,2.9,3.3)); 
 addHistogram( new TH1F("J_px",  ";;Events",200,-100,100)); 
 addHistogram( new TH1F("J_py",  ";;Events",200,-100,100)); 
 addHistogram( new TH1F("J_pz",  ";;Events",200,-100,100)); 
 addHistogram( new TH1F("J_px1", ";;Events",200,-100,100)); 
 addHistogram( new TH1F("J_py1", ";;Events",200,-100,100)); 
 addHistogram( new TH1F("J_pz1", ";;Events",200,-100,100)); 
 addHistogram( new TH1F("J_charge1",";;Events",2,-1,1)); 
 addHistogram( new TH1F("J_px2",";;Events",200,-100,100)); 
 addHistogram( new TH1F("J_py2",";;Events",200,-100,100)); 
 addHistogram( new TH1F("J_pz2",";;Events",200,-100,100)); 
 addHistogram( new TH1F("J_charge2",";;Events",2,-1,1)); 
 addHistogram( new TH1F("mumC2",";;Events",20,0,10)); 
 addHistogram( new TH1F("mumNHits",";;Events",35,0,35)); 
 addHistogram( new TH1F("mumNPHits",";;Events",10,0,10)); 
 addHistogram( new TH1F("mupC2",";;Events",30,0,10)); 
 addHistogram( new TH1F("mupNHits",";;Events",35,0,35)); 
 addHistogram( new TH1F("mupNPHits",";;Events",10,0,10)); 
 addHistogram( new TH1F("mumdxy",";;Events",100,-0.05,0.05)); 
 addHistogram( new TH1F("mupdxy",";;Events",100,-0.05,0.05)); 
 addHistogram( new TH1F("mumdz",";;Events",100,-0.05,0.05)); 
 addHistogram( new TH1F("mupdz",";;Events",100,-0.05,0.05)); 
 addHistogram( new TH1F("muon_dca",";;Events",10,0,0.1)); 
 addHistogram( new TH1F("mu1soft",";;Events",2,0,2)); 
 addHistogram( new TH1F("mu2soft",";;Events",2,0,2)); 
 addHistogram( new TH1F("mu1tight",";;Events",2,0,2)); 
 addHistogram( new TH1F("mu2tight",";;Events",2,0,2)); 
 addHistogram( new TH1F("mu1PF",";;Events",2,0,2)); 
 addHistogram( new TH1F("mu2PF",";;Events",2,0,2)); 
 addHistogram( new TH1F("mu1loose",";;Events",2,0,2)); 
 addHistogram( new TH1F("mu2loose",";;Events",2,0,2)); 
 addHistogram( new TH1F("J_lxy",";;Events",100,0,0.1)); 
 addHistogram( new TH1F("J_lxyErr",";;Events",100,0,0.005)); 
 addHistogram( new TH1F("J_vertexchi2",";;Events",20,0,4)); 
 addHistogram( new TH1F("Pi_dJP",";;Events",10,-0.1,0.1)); 
 addHistogram( new TH1F("JPi_lxy",";;Events",100,0,0.1)); 
 addHistogram( new TH1F("JPi_lxyErr",";;Events",100,0,0.0002)); 
 addHistogram( new TH1F("JPiPi_lxy",";;Events",100,0,0.1)); 
 addHistogram( new TH1F("JPiPi_lxyErr",";;Events",100,0,0.0002)); 
 addHistogram( new TH1F("JPiPi_x",";;Events",100,-0.05,0.05)); 
 addHistogram( new TH1F("JPiPi_y",";;Events",100,-0.05,0.05)); 
 addHistogram( new TH1F("JPiPi_z",";;Events",100,-0.05,0.05)); 
 addHistogram( new TH1F("Pi_nhits1",";;Events",40,0,40)); 
 addHistogram( new TH1F("Pi_npixelhits1",";;Events",10,0,10)); 
 addHistogram( new TH1F("Pi_nhits2",";;Events",40,0,40)); 
 addHistogram( new TH1F("Pi_npixelhits2",";;Events",10,0,10)); 
 addHistogram( new TH1F("Pi_eta1",";;Events",60,-3,3)); 
 addHistogram( new TH1F("Pi_eta2",";;Events",60,-3,3)); 
 addHistogram( new TH1F("Pi_phi1",";;Events",80,-4,4)); 
 addHistogram( new TH1F("Pi_phi2",";;Events",40,-4,4)); 
 addHistogram( new TH1F("Pi_pt1",";;Events",40,0,20)); 
 addHistogram( new TH1F("Pi_pt2",";;Events",40,0,20)); 
 addHistogram( new TH1F("Pi_e1",";;Events",40,0,20)); 
 addHistogram( new TH1F("Pi_e2",";;Events",40,0,20)); 
 addHistogram( new TH1F("Pi_charge1",";;Events",2,-1,1)); 
 addHistogram( new TH1F("Pi_charge2",";;Events",2,-1,1)); 
 addHistogram( new TH1F("Pi_dxy1",";;Events",100,0,0.1)); 
 addHistogram( new TH1F("Pi_dxy2",";;Events",100,0,0.1)); 
 addHistogram( new TH1F("Pi_dxyerr1",";;Events",100,0,0.1)); 
 addHistogram( new TH1F("Pi_dxyerr2",";;Events",100,0,0.1)); 
 addHistogram( new TH1F("Pi_vertexchisq1",";;Events",50,0,20)); 
 addHistogram( new TH1F("Pi_vertexchisq2",";;Events",50,0,20)); 
 return true;           
}

template<class T>
bool SmartSelectionMonitor_hzz::fillHistoForAllCategories(TString name, double variable, T currentEvt, TString tag, double weight){
  fillHisto(name, tag+currentEvt.s_jetCat+currentEvt.s_lepCat, variable, weight, true);
  fillHisto(name, tag+currentEvt.s_lepCat, variable, weight, true); //all jet cats. No need for all lep cats since the s_lepCat tag already contains "_ll".
  return true;
}

template<class T>
bool SmartSelectionMonitor_hzz::fillProfileForAllCategories(TString name, double variableY, double variableX, T currentEvt, TString tag, double weight){
  fillProfile(name, tag+currentEvt.s_jetCat+currentEvt.s_lepCat, variableY, variableX, weight);
  fillProfile(name, tag+currentEvt.s_lepCat, variableY, variableX, weight); //all jet cats. No need for all lep cats since the s_lepCat tag already contains "_ll".
  return true;
}

template<class T>
bool SmartSelectionMonitor_hzz::fillAnalysisHistos_common(T currentEvt, TString tag, double weight){
  std::map<std::string, double> data;
  data["mT"] = currentEvt.MT;
  data["M_Boson"] = currentEvt.M_Boson;
  data["pT_Boson"] = currentEvt.pT_Boson;
  data["MET"] = currentEvt.MET;
  data["METphi"] = currentEvt.METphi;
  data["nJets"] = currentEvt.nJets;
  data["eta_Boson"] = currentEvt.eta_Boson;
  data["DeltaPhi_MET_Jet"] = currentEvt.deltaPhi_MET_Jet;
  data["DeltaPhi_MET_Boson"] = currentEvt.deltaPhi_MET_Boson;
  data["reco-vtx"] = currentEvt.nVtx;
  for(std::map<std::string,double>::iterator it = data.begin() ; it != data.end() ; it++) fillHistoForAllCategories(it->first, it->second, currentEvt, tag, weight);
  return true;
}

bool SmartSelectionMonitor_hzz::fillAnalysisHistos(photon_evt currentEvt, TString tag, double weight){
  fillAnalysisHistos_common(currentEvt, tag, weight);
  return true;
}

bool SmartSelectionMonitor_hzz::fillAnalysisHistos(evt currentEvt, TString tag, double weight){
  fillAnalysisHistos_common(currentEvt, tag, weight);
  std::map<std::string, double> data;
  if(tag.Contains("inZpeak")){
    data["pT_l1"] = ((evt) currentEvt).lep1pT;
    data["eta_l1"] = currentEvt.lep1eta;
    data["pT_l2"] = currentEvt.lep2pT;
    data["eta_l2"] = currentEvt.lep2eta;
    data["runNumber"] = currentEvt.runNumber;
  }
  for(std::map<std::string,double>::iterator it = data.begin() ; it != data.end() ; it++) fillHistoForAllCategories(it->first, it->second, currentEvt, tag, weight);
  return true;
}

bool SmartSelectionMonitor_hzz::fillPhotonIDHistos_InstrMET(photon_evt currentEvt, TString tag, double weight){
  fillAnalysisHistos(currentEvt, tag, weight);

  std::map<std::string, double> data;
  data["METoverPt"] = currentEvt.METoPT;
  data["METoverPt_zoom"] = currentEvt.METoPT;
  data["METpar"] = currentEvt.METpar;
  data["METperp"] = currentEvt.METperp;
  data["METsig"] = currentEvt.METsig;
  data["nvtx"] = currentEvt.nVtx;
  data["HoE"] = currentEvt.HoE;
  data["HoE_zoom"] = currentEvt.HoE;
  data["SigmaIetaIeta"] = currentEvt.sigmaIEtaIEta;
  data["SigmaIetaIeta_zoom"] = currentEvt.sigmaIEtaIEta;
  data["RhoCorrPfIsoChHad"] = currentEvt.chIsoRhoCorr;
  data["RhoCorrPfIsoChHad_zoom"] = currentEvt.chIsoRhoCorr;
  data["RhoCorrPfIsoNeutralHad"] = currentEvt.neuIsoRhoCorr;
  data["RhoCorrPfIsoNeutralHad_zoom"] = currentEvt.neuIsoRhoCorr;
  data["RhoCorrPfIsoPhot"] = currentEvt.phoIsoRhoCorr;
  data["RhoCorrPfIsoPhot_zoom"] = currentEvt.phoIsoRhoCorr;
  data["R9"] = currentEvt.R9;
  data["rho"] = currentEvt.rho;
  data["pT_jet0"] = currentEvt.jet0_pT;
  data["pT_jet1"] = currentEvt.jet1_pT;
  data["pT_jet2"] = currentEvt.jet2_pT;
  data["pT_jet3"] = currentEvt.jet3_pT;
  data["selJetsHT"] = currentEvt.HT_selJets;
  for(std::map<std::string,double>::iterator it = data.begin() ; it != data.end() ; it++) fillHistoForAllCategories(it->first, it->second, currentEvt, tag, weight);
  return true;
}

void SmartSelectionMonitor_hzz::WriteForSysts(TString systName, bool keepEverything){
  TString systNameToAppend = "";
  if(systName != "") systNameToAppend = "_"+systName;
  for(SmartSelectionMonitor::Monitor_t::iterator it = SmartSelectionMonitor::allMonitors_.begin(); it!= SmartSelectionMonitor::allMonitors_.end(); it++){
    std::map<TString, TH1*>* map = it->second;
    for(std::map<TString, TH1*>::iterator h =map->begin(); h!= map->end(); h++){
      if(h->first!="all") h->second->SetName(h->second->GetName() + systNameToAppend);
      if(h->first!="all") h->second->SetTitle(h->second->GetName());
      if(keepEverything && h->first!="all")h->second->Write();
      std::vector<TString> jetCat = {"_eq0jets","_geq1jets","_vbf"};  std::vector<TString> lepCat = {"_ee","_mumu"};
      if(!keepEverything){
        if(h->second->GetName()==TString("totEventInBaobab_tot"))h->second->Write();
        for(unsigned int i = 0; i < jetCat.size(); i++){
          for(unsigned int j = 0; j < lepCat.size(); j++){
            if(h->second->GetName()=="mT_final"+jetCat[i]+lepCat[j]+systNameToAppend)h->second->Write();
            if(h->second->GetName()=="MET_InstrMET_reweighting"+jetCat[i]+lepCat[j])h->second->Write();
          }
        }
      }
    }
  }
}

//Histo used for closure Test and check of Instr. MET
//template<class T>
bool SmartSelectionMonitor_hzz::fillInstrMETControlRegionHisto(base_evt currentEvt, TString tag, double weight){
  std::map<std::string, double> histo;
  histo["reco-vtx"] = currentEvt.nVtx;
  histo["pT_Boson"] = currentEvt.pT_Boson;
  histo["M_Boson"] = currentEvt.M_Boson;
  histo["MET"] = currentEvt.MET;
  histo["mT"] = currentEvt.MT;
  histo["DeltaPhi_MET_Boson"] = currentEvt.deltaPhi_MET_Boson;
  histo["DeltaPhi_MET_Jet"] = currentEvt.deltaPhi_MET_Jet;
  histo["METoverPt_zoom"] = currentEvt.METoPT;
  histo["eta_Boson"] = currentEvt.eta_Boson;
  histo["pT_jet0"] = currentEvt.jet0_pT;
  histo["nJets"] = currentEvt.nJets;
  histo["selJetsHT"] = currentEvt.HT_selJets;
  if(currentEvt.HT_selJets >300) histo["MET_HT300"] = currentEvt.MET;
  if(currentEvt.pT_Boson < 300) histo["MET_Pt0-300"] = currentEvt.MET;
  else if(currentEvt.pT_Boson < 400) histo["MET_Pt300-400"] = currentEvt.MET;
  else if(currentEvt.pT_Boson < 600) histo["MET_Pt400-600"] = currentEvt.MET;
  else histo["MET_Pt600-Inf"] = currentEvt.MET;
  
  std::map<std::string, std::pair<double,double> > profile;
  profile["METvsBosonPt"] = {currentEvt.pT_Boson, currentEvt.MET};
  profile["METvsMT"] = {currentEvt.MT, currentEvt.MET};
  profile["METvsDPhiMETBos"] = {currentEvt.deltaPhi_MET_Boson, currentEvt.MET};
  profile["METvsDPhiMETJet"] = {currentEvt.deltaPhi_MET_Jet, currentEvt.MET};
  profile["METvsJetPt"] = {currentEvt.jet0_pT, currentEvt.MET};
  profile["METvsNJets"] = {currentEvt.nJets, currentEvt.MET};
  profile["METvsBosonEta"] = {currentEvt.eta_Boson ,currentEvt.MET};
  profile["METvsHT"] = {currentEvt.HT_selJets ,currentEvt.MET};
  profile["HTvsBosonEta"] = {currentEvt.eta_Boson , currentEvt.HT_selJets};
  profile["HTvsBosonPt"] = {currentEvt.pT_Boson , currentEvt.HT_selJets};

  for(std::map<std::string, std::pair<double,double> >::iterator it = profile.begin() ; it != profile.end() ; it++) fillProfileForAllCategories(it->first, it->second.first, it->second.second, currentEvt, tag, weight);
  for(std::map<std::string,double>::iterator it = histo.begin() ; it != histo.end() ; it++) fillHistoForAllCategories(it->first, it->second, currentEvt, tag, weight);
  return true;
}
