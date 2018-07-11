#include "../Includes/SmartSelectionMonitor_hzz.h"

bool SmartSelectionMonitor_hzz::declareHistos(){ //FIXME: Later, will take an array as input for the binnings.
  addHistogram(new TH1F("totEventInBaobab",";Number of events in Baobab;Events",150,0,150));
  addHistogram(new TH1F("pile-up",";Number of PU events;Events",100,0,100));
  addHistogram(new TH1F("truth-pile-up",";Truth number of PU events;Events",100,0,100));
  addHistogram(new TH1F("reco-vtx",";Number of reco vtx;Events",100,0,100)); //Use for Photon reweighting method, don't change binning if you don't know what you're doing

  TH1F *h =(TH1F*) addHistogram(new TH1F("eventflow",";;Events",10,0,10));
  h->GetXaxis()->SetBinLabel(1,"skimmed");
  h->GetXaxis()->SetBinLabel(2,"#geq 2 iso leptons");
  h->GetXaxis()->SetBinLabel(3,"|M-91|<15");
  h->GetXaxis()->SetBinLabel(4,"p_{T}>55");
  h->GetXaxis()->SetBinLabel(5,"#Delta #phi(Z,E_{T}^{miss})>0.5");
  h->GetXaxis()->SetBinLabel(6,"3^{rd}-lepton veto");
  h->GetXaxis()->SetBinLabel(7,"b-veto");
  h->GetXaxis()->SetBinLabel(8,"#Delta #phi(jet,E_{T}^{miss})>0.5");
  h->GetXaxis()->SetBinLabel(9,"E_{T}^{miss}>80");
  h->GetXaxis()->SetBinLabel(10,"E_{T}^{miss}>125");
  addHistogram(new TH1F("pT_mu",";p_{T} of muon;Events",200,0,800));
  addHistogram(new TH1F("pT_e",";p_{T} of electron;Events",200,0,800));
  addHistogram(new TH1F("nb_mu",";number of muons;Events",10,0,10));
  addHistogram(new TH1F("nb_e",";number of electrons;Events",10,0,10));

  Double_t METaxis[]={0,5,10,15,20,25,30,35,40,45,50,55,60,70,80,90,100,125,150,175,200,250,300,400,500,600,700,800,900,1000};
  Double_t zptaxis[]= {0,15,30,45,60,75,90,105,120,135,150,165,180,195,210,225,240,255,270,285,300,315,330,345,360,375,390,405,435,465,495,525,555,585,615,675,735,795,855,975,1500};
  Double_t mTaxis[]={100,120,140,160,180,200,220,240,260,280,300,325,350,375,400,450,500,600,700,800,900,1000,1500,2000};
  Int_t nMETAxis=sizeof(METaxis)/sizeof(Double_t);
  Int_t nzptAxis=sizeof(zptaxis)/sizeof(Double_t);
  Int_t nmTAxis=sizeof(mTaxis)/sizeof(Double_t);
  addHistogram(new TH1F("M_Boson",";M_{Z};Events",100,76,106));
  addHistogram(new TH1F("pT_Boson",";p_{T,Z};Events",nzptAxis-1,zptaxis));  //Use for Photon reweighting method, don't change binning if you don't know what you're doing
  addHistogram(new TH1F("MET",";Missing transverse energy (GeV);Events",nMETAxis-1,METaxis));
  addHistogram(new TH1F("METphi",";#phi of missing transverse energy;Events",80,-4.,4.));
  addHistogram(new TH1F("mT",";m_{T} (GeV);Events",nmTAxis-1,mTaxis));

  TH1F *hc = (TH1F*) addHistogram(new TH1F("jetCategory",";Jet Category;Events",3,0,3));
  hc->GetXaxis()->SetBinLabel(1,"= 0 jets");
  hc->GetXaxis()->SetBinLabel(2,"#geq 1 jets");
  hc->GetXaxis()->SetBinLabel(3,"vbf");
  addHistogram(new TH1F("nJets",";N_{jets #geq 30 GeV};Events",20,0,20)); //Jets using the same criteria as the ones for jet bin category

  addHistogram(new TH1F("pT_l1",";p_{T} of lepton 1;Events",16,10,250));
  addHistogram(new TH1F("pT_l2",";p_{T} of lepton 2;Events",16,10,250));
  addHistogram(new TH1F("eta_l1",";#eta of lepton 1;Events",20,-2.5,2.5));
  addHistogram(new TH1F("eta_l2",";#eta of lepton 2;Events",20,-2.5,2.5));

  addHistogram(new TH1F("runNumber",";run number",100,273158, 284044));

  TH1F *h_metFilter =(TH1F*) addHistogram(new TH1F("metFilters",";;Events passing MET filters",26+1,0,26+1)); //We add +1 everywhere so we can start with bin 0
  h_metFilter->GetXaxis()->SetBinLabel(26+1,"all"); //1 is already taken so we take the last bin of the MET filters +1, i.e 25+1=26
  h_metFilter->GetXaxis()->SetBinLabel(0+1,"duplicateMuons");
  h_metFilter->GetXaxis()->SetBinLabel(1+1,"badMuons");
  h_metFilter->GetXaxis()->SetBinLabel(2+1,"noBadMuons");
  h_metFilter->GetXaxis()->SetBinLabel(14+1,"primary vertex filter");
  h_metFilter->GetXaxis()->SetBinLabel(8+1,"beam halo filter");
  h_metFilter->GetXaxis()->SetBinLabel(3+1,"HBHE noise filter");
  h_metFilter->GetXaxis()->SetBinLabel(4+1,"HBHEiso noise filter");
  h_metFilter->GetXaxis()->SetBinLabel(12+1,"ECAL TP filter");
  h_metFilter->GetXaxis()->SetBinLabel(15+1,"ee badSC noise filter");
  h_metFilter->GetXaxis()->SetBinLabel(24+1,"badPFMuon");
  h_metFilter->GetXaxis()->SetBinLabel(25+1,"badCharged hadron");

  addHistogram(new TH1F("custom_HT",";#Sigma p_{T}^{gen jets};Events",200,0,2000));
  addHistogram(new TH1F("selJetsHT",";#Sigma p_{T}^{sel jets};Events",200,0,2000));

  //Control histo for closure test
  addHistogram(new TH1F("DeltaPhi_MET_Boson",";#Delta #phi(Z,E_{T}^{miss});Events", 50, 0, 4));
  addHistogram(new TH1F("DeltaPhi_MET_Jet",";min(#Delta#phi(jet,E_{T}^{miss}));Events", 50, 0, 4));
  addHistogram(new TH1F("METoverPt_zoom",";MET/p^{Z}_{T};Events",40,0, 2));
  addHistogram(new TH1F("eta_Boson",";#eta_{Z};Events",80, -4, 4));
  addHistogram(new TH1F("pT_jet0",";p_{T} of jet 0;Events",200,0,800));
  addHistogram(new TH1F("MET_HT300",";Missing transverse energy (GeV);Events",nMETAxis-1,METaxis));
  addHistogram(new TH1F("MET_Pt0-300",";Missing transverse energy (GeV);Events (Pt0-300)",nMETAxis-1,METaxis));
  addHistogram(new TH1F("MET_Pt300-400",";Missing transverse energy (GeV);Events (Pt300-400)",nMETAxis-1,METaxis));
  addHistogram(new TH1F("MET_Pt400-600",";Missing transverse energy (GeV);Events (Pt400-600)",nMETAxis-1,METaxis));
  addHistogram(new TH1F("MET_Pt600-Inf",";Missing transverse energy (GeV);Events (Pt600-Inf)",nMETAxis-1,METaxis));

  

  //TProfile for closure
  addHistogram(new TProfile("METvsBosonPt",";p^{Z}_{T} (GeV);MET profile (GeV)", nzptAxis-1, zptaxis, 0, 500));
  addHistogram(new TProfile("METvsMT",";m_{T} (GeV);MET profile (GeV)", nmTAxis-1, mTaxis, 0, 500));
  addHistogram(new TProfile("METvsDPhiMETBos",";#Delta #phi(Z,E_{T}^{miss});MET profile (GeV)", 20, 0, 4, 0, 500));
  addHistogram(new TProfile("METvsDPhiMETJet",";min(#Delta#phi(jet,E_{T}^{miss}));MET profile (GeV)",20, 0, 4, 0, 500));
  addHistogram(new TProfile("METvsJetPt",";leading jet p_{T} (GeV);MET profile (GeV)", 40, 0, 800, 0, 500));
  addHistogram(new TProfile("METvsNJets",";# jets;MET profile (GeV)", 10, 0, 10, 0, 500));
  addHistogram(new TProfile("METvsBosonEta",";Z #eta;MET profile (GeV)", 40, -4, 4, 0, 500));
  addHistogram(new TProfile("METvsHT",";HT;MET profile (GeV)", 15, 0, 1500, 0, 500));
  addHistogram(new TProfile("HTvsBosonEta",";Z #eta;HT profile (GeV)", 40, -4, 4, 0, 500));

  return true;
}

bool SmartSelectionMonitor_hzz::declareHistos_InstrMET(){ 
  addHistogram(new TH1F("totEventInBaobab",";Number of events in Baobab;Events",150,0,150));
  addHistogram(new TH1F("pile-up",";Number of PU events;Events",100,0,100));
  addHistogram(new TH1F("truth-pile-up",";Truth number of PU events;Events",100,0,100));
  addHistogram(new TH1F("reco-vtx",";Number of reco vtx;Events",100,0,100));   //Use for Photon reweighting method, don't change binning if you don't know what you're doing
  TH1F *h =(TH1F*) addHistogram(new TH1F("eventflow",";;Events",17,0,17));
  h->GetXaxis()->SetBinLabel(1,"skimmed (ID, p_{T}, trigger, ...)"); //what is coming from the bonzai
  h->GetXaxis()->SetBinLabel(2,"prescale weight");
  h->GetXaxis()->SetBinLabel(3,"photon efficiency SF");
  h->GetXaxis()->SetBinLabel(4,"PU reweighting");
  h->GetXaxis()->SetBinLabel(5,"MET filters");
  h->GetXaxis()->SetBinLabel(6,"QCD/GJets double counting");
  h->GetXaxis()->SetBinLabel(7,"ZnunuG k-factor");
  h->GetXaxis()->SetBinLabel(8,"WJets double counting");
  h->GetXaxis()->SetBinLabel(9,"p_{T, #gamma}>55");
  h->GetXaxis()->SetBinLabel(10,"#Delta #phi(#gamma,E_{T}^{miss})>0.5");
  h->GetXaxis()->SetBinLabel(11,"no lepton veto");
  h->GetXaxis()->SetBinLabel(12,"Dilepton NVtx Reweighting");
  h->GetXaxis()->SetBinLabel(13,"Dilepton Pt Reweighting");
  h->GetXaxis()->SetBinLabel(14,"b-veto");
  h->GetXaxis()->SetBinLabel(15,"#Delta #phi(jet,E_{T}^{miss})>0.5");
  h->GetXaxis()->SetBinLabel(16,"E_{T}^{miss}>80");
  h->GetXaxis()->SetBinLabel(17,"E_{T}^{miss}>125");
  addHistogram(new TH1F("custom_HT",";#Sigma p_{T}^{gen jets};Events",200,0,2000));
  addHistogram(new TH1F("selJetsHT",";#Sigma p_{T}^{sel jets};Events",200,0,2000));
  addHistogram(new TH1F("pT_jet0",";p_{T} of jet 0;Events",200,0,800));
  addHistogram(new TH1F("pT_jet1",";p_{T} of jet 1;Events",200,0,800));
  addHistogram(new TH1F("pT_jet2",";p_{T} of jet 2;Events",200,0,800));
  addHistogram(new TH1F("pT_jet3",";p_{T} of jet 3;Events",200,0,800));

  Double_t METaxis[]={0,5,10,15,20,25,30,35,40,45,50,55,60,70,80,90,100,125,150,175,200,250,300,400,500,600,700,800,900,1000};
  Double_t zptaxis[]= {0,15,30,45,60,75,90,105,120,135,150,165,180,195,210,225,240,255,270,285,300,315,330,345,360,375,390,405,435,465,495,525,555,585,615,675,735,795,855,975,1500};  //Use for Photon reweighting method, don't change binning if you don't know what you're doing
  Double_t mTaxis[]={100,120,140,160,180,200,220,240,260,280,300,325,350,375,400,450,500,600,700,800,900,1000,1500,2000};
  Int_t nMETAxis=sizeof(METaxis)/sizeof(Double_t);
  Int_t nzptAxis=sizeof(zptaxis)/sizeof(Double_t);
  Int_t nmTAxis=sizeof(mTaxis)/sizeof(Double_t);
  addHistogram(new TH1F("M_Boson",";M_{#gamma};Events",100,76,106));
  addHistogram(new TH1F("pT_Boson",";p_{T,#gamma};Events",nzptAxis-1,zptaxis));
  addHistogram(new TH1F("pT_Boson_unif",";p_{T,#gamma};Events",300, 0, 1500));
  addHistogram(new TH1F("eta_Boson",";#eta_{#gamma};Events",80, -4, 4));
  addHistogram(new TH1F("MET",";Missing transverse energy (GeV);Events",nMETAxis-1,METaxis));
  addHistogram(new TH1F("MET_unif",";Missing transverse energy (GeV);Events",200,0,1000));
  addHistogram(new TH1F("mT",";m_{T} (GeV);Events",nmTAxis-1,mTaxis));
  addHistogram(new TH1F("mT_unif",";m_{T} (GeV);Events",300,0,1500));

  TH1F *hc = (TH1F*) addHistogram(new TH1F("jetCategory",";Jet Category;Events",3,0,3));
  hc->GetXaxis()->SetBinLabel(1,"= 0 jets");
  hc->GetXaxis()->SetBinLabel(2,"#geq 1 jets");
  hc->GetXaxis()->SetBinLabel(3,"vbf");
  addHistogram(new TH1F("nJets",";N_{jets #geq 30 GeV};Events",20,0,20)); //Jets using the same criteria as the ones for jet bin category

  TH1F *h_metFilter =(TH1F*) addHistogram(new TH1F("metFilters",";;Events passing MET filters",26+1,0,26+1)); //We add +1 everywhere so we can start with bin 0
  h_metFilter->GetXaxis()->SetBinLabel(26+1,"all"); //1 is already taken so we take the last bin of the MET filters +1, i.e 25+1=26
  h_metFilter->GetXaxis()->SetBinLabel(0+1,"duplicateMuons");
  h_metFilter->GetXaxis()->SetBinLabel(1+1,"badMuons");
  h_metFilter->GetXaxis()->SetBinLabel(2+1,"noBadMuons");
  h_metFilter->GetXaxis()->SetBinLabel(14+1,"primary vertex filter");
  h_metFilter->GetXaxis()->SetBinLabel(8+1,"beam halo filter");
  h_metFilter->GetXaxis()->SetBinLabel(3+1,"HBHE noise filter");
  h_metFilter->GetXaxis()->SetBinLabel(4+1,"HBHEiso noise filter");
  h_metFilter->GetXaxis()->SetBinLabel(12+1,"ECAL TP filter");
  h_metFilter->GetXaxis()->SetBinLabel(15+1,"ee badSC noise filter");
  h_metFilter->GetXaxis()->SetBinLabel(24+1,"badPFMuon");
  h_metFilter->GetXaxis()->SetBinLabel(25+1,"badCharged hadron");

  //Compa old new
  addHistogram(new TH1F("qt_rebin",      ";Transverse momentum [GeV];Events / GeV",nzptAxis-1,zptaxis));
  addHistogram(new TH1F("nvtx",";Vertices;Events",100, 0, 100)); //50,0,50) );
  addHistogram(new TH1F("rho",";#rho;Events",100,0,50) );
  addHistogram(new TH1F("MET_phi",";<MET #phi>;Events", 25, -4, 4));
  addHistogram(new TH1F("DeltaPhi_MET_Boson",";#Delta #phi(#gamma,E_{T}^{miss});Events", 50, 0, 4));
  addHistogram(new TH1F("DeltaPhi_MET_Jet",";min(#Delta#phi(jet,E_{T}^{miss}));Events", 50, 0, 4));

  //Photon Study
  //R9
  addHistogram( new TH1F ("R9", ";R9;Events", 150,0,1.5));
  //Sigma ieta ieta
  addHistogram( new TH1F ("SigmaIetaIeta", ";#sigma_{i#eta i#eta};Events", 50, 0, 0.05));
  addHistogram( new TH1F ("SigmaIetaIeta_zoom", ";#sigma_{i#eta i#eta};Events", 40, 0, 0.02));
  //Isolation
  Double_t PfIsoChHadBins[] = {1, 2, 10, 50};
  Int_t NPfIsoChHadBins = sizeof(PfIsoChHadBins)/sizeof(PfIsoChHadBins[0]) - 1;
  addHistogram( new TH1F ("PfIsoChHad_showQCD", ";PfIsoChHad;Events", NPfIsoChHadBins, PfIsoChHadBins));
  addHistogram( new TH1F ("PfIsoChHad", ";PfIsoChHad;Events", 50, 0, 10));
  addHistogram( new TH1F ("PfIsoNeutralHad", ";PfIsoNeutralHad;Events", 100, 0, 20));
  addHistogram( new TH1F ("PfIsoPhot", ";PfIsoPhot;Events", 100, 0, 20));
  addHistogram( new TH1F ("IsoEcal", ";IsoEcal;Events", 150, 0, 60));
  addHistogram( new TH1F ("IsoHcal", ";IsoHcal;Events", 100, 0, 40));
  addHistogram( new TH1F ("IsoTk", ";IsoTk;Events", 100, 0, 40));
  addHistogram( new TH1F ("PfIsoPuChHad", ";PfIsoPuChHad;Events", 150, 0, 75));
  addHistogram( new TH1F ("PfIsoEcalClus", ";PfIsoEcalClus;Events", 100, 0, 20));
  addHistogram( new TH1F ("PfIsoHcalClus", ";PfIsoHcalClus;Events", 100, 0, 20));
  //HoverE
  addHistogram( new TH1F ("HoE", ";HoE;Events", 40, 0, 0.20));
  addHistogram( new TH1F ("HoE_zoom", ";HoE;Events", 30, 0, 0.03));
  //Rho corrected photon ID Iso
  addHistogram( new TH1F ("RhoCorrPfIsoChHad", ";Rho Corrected PfIsoChHad;Events", 50, 0, 10));
  addHistogram( new TH1F ("RhoCorrPfIsoChHad_zoom", ";Rho Corrected PfIsoChHad;Events", 50, 0, 0.5));
  addHistogram( new TH1F ("RhoCorrPfIsoNeutralHad", ";Rho Corrected PfIsoNeutralHad;Events", 100, 0, 20));
  addHistogram( new TH1F ("RhoCorrPfIsoNeutralHad_zoom", ";Rho Corrected PfIsoNeutralHad;Events", 100, 0, 20));
  addHistogram( new TH1F ("RhoCorrPfIsoPhot", ";Rho Corrected PfIsoPhot;Events", 100, 0, 20));
  addHistogram( new TH1F ("RhoCorrPfIsoPhot_zoom", ";Rho Corrected PfIsoPhot;Events", 100, 0, 20));

  //More MET variables
  addHistogram(new TH1F("METsigx2",";Significance x^2 (C(0,0)) of the Missing transverse energy (GeV);Events",500,0, 1000));
  addHistogram(new TH1F("METsigxy",";Significance x-y (C(0,1)) of the Missing transverse energy (GeV);Events",500,0, 1000));
  addHistogram(new TH1F("METsigy2",";Significance y^2 (C(1,1)) of the Missing transverse energy (GeV);Events",500,0, 1000));
  addHistogram(new TH1F("METsig",";Significance of the Missing transverse energy (GeV);Events",100,0, 200));
  addHistogram(new TH1F("METoverPt",";MET/p^{boson}_{T};Events",100,0, 100));
  addHistogram(new TH1F("METoverPt_zoom",";MET/p^{boson}_{T};Events",40,0, 2));
  addHistogram(new TH1F("METperp","MET_{#perp}", 100, -200, 200)); 
  addHistogram(new TH1F("METpar","MET_{#parallel}+q_{T}", 100, -200, 200)); 

  addHistogram(new TH1F("MET_HT300",";Missing transverse energy (GeV);Events (HT>300GeV)",nMETAxis-1,METaxis));
  addHistogram(new TH1F("MET_Pt0-300",";Missing transverse energy (GeV);Events (Pt0-300)",nMETAxis-1,METaxis));
  addHistogram(new TH1F("MET_Pt300-400",";Missing transverse energy (GeV);Events (Pt300-400)",nMETAxis-1,METaxis));
  addHistogram(new TH1F("MET_Pt400-600",";Missing transverse energy (GeV);Events (Pt400-600)",nMETAxis-1,METaxis));
  addHistogram(new TH1F("MET_Pt600-Inf",";Missing transverse energy (GeV);Events (Pt600-Inf)",nMETAxis-1,METaxis));

  //TProfile for closure
  addHistogram(new TProfile("METvsBosonPt","; p^{#gamma}_{T} (GeV);MET profile (GeV)", nzptAxis-1, zptaxis, 0, 500));
  addHistogram(new TProfile("METvsMT",";m_{T} (GeV);MET profile (GeV)", nmTAxis-1, mTaxis, 0, 500));
  addHistogram(new TProfile("METvsDPhiMETBos",";#Delta #phi(#gamma,E_{T}^{miss});MET profile (GeV)", 20, 0, 4, 0, 500));
  addHistogram(new TProfile("METvsDPhiMETJet",";min(#Delta#phi(jet,E_{T}^{miss}));MET profile (GeV)",20, 0, 4, 0, 500));
  addHistogram(new TProfile("METvsJetPt",";leading jet p_{T} (GeV);MET profile (GeV)", 40, 0, 800, 0, 500));
  addHistogram(new TProfile("METvsNJets",";# jets;MET profile (GeV)", 10, 0, 10, 0, 500));
  addHistogram(new TProfile("METvsBosonEta",";#gamma #eta;MET profile (GeV)", 40, -4, 4, 0, 500));
  addHistogram(new TProfile("METvsHT",";HT;MET profile (GeV)", 15, 0, 1500, 0, 500));
  addHistogram(new TProfile("HTvsBosonEta",";#gamma #eta;HT profile (GeV)", 40, -4, 4, 0, 500));

  return true;  
} 

bool SmartSelectionMonitor_hzz::declareHistos_NRB(){ 
  Double_t mtaxis[]={100,120,140,160,180,200,220,240,260,280,300,325,350,375,400,450,500,600,700,800,900,1000,1500,2000};
  Int_t nmtAxis=sizeof(mtaxis)/sizeof(Double_t);
  Double_t metaxis[]={0,5,10,15,20,25,30,35,40,45,50,55,60,70,80,90,100,125,150,175,200,250,300,400,500,600,700,800,900,1000};
  Int_t nmetAxis=sizeof(metaxis)/sizeof(Double_t);
  addHistogram( new TH1F( "zmass_btag50", ";Mass [GeV];Events / 2 GeV", 100,40,200) );
  addHistogram( new TH1F( "zmass_bveto50",";Mass [GeV];Events / 2 GeV", 100,40,200) );
  addHistogram( new TH1F( "zmass_btag80", ";Mass [GeV];Events / 2 GeV", 100,40,200) );
  addHistogram( new TH1F( "zmass_bveto80",";Mass [GeV];Events / 2 GeV", 100,40,200) );
  addHistogram( new TH1F( "zmass_btag125", ";Mass [GeV];Events / 2 GeV", 100,40,200) );
  addHistogram( new TH1F( "zmass_bveto125",";Mass [GeV];Events / 2 GeV", 100,40,200) );
  addHistogram( new TH1F( "met_Inbtag",          ";Missing transverse energy [GeV];Events / GeV",nmetAxis-1,metaxis) ); //50,0,1000) );
  addHistogram( new TH1F( "met_Inbveto",          ";Missing transverse energy [GeV];Events / GeV",nmetAxis-1,metaxis) ); //50,0,1000) );
  addHistogram( new TH1F( "met_Outbtag",          ";Missing transverse energy [GeV];Events / GeV",nmetAxis-1,metaxis) ); //50,0,1000) );
  addHistogram( new TH1F( "met_Outbveto",          ";Missing transverse energy [GeV];Events / GeV",nmetAxis-1,metaxis) ); //50,0,1000) );
  addHistogram( new TH1F( "mt_Inbtag50"  ,         ";Transverse mass [GeV];Events / GeV",nmtAxis-1,mtaxis) );
  addHistogram( new TH1F( "mt_Inbveto50"  ,         ";Transverse mass [GeV];Events / GeV",nmtAxis-1,mtaxis) );
  addHistogram( new TH1F( "mt_Inbtag80"  ,         ";Transverse mass [GeV];Events / GeV",nmtAxis-1,mtaxis) );
  addHistogram( new TH1F( "mt_Inbveto80"  ,         ";Transverse mass [GeV];Events / GeV",nmtAxis-1,mtaxis) );
  addHistogram( new TH1F( "mt_Inbtag125"  ,         ";Transverse mass [GeV];Events / GeV",nmtAxis-1,mtaxis) );
  addHistogram( new TH1F( "mt_Inbveto125"  ,         ";Transverse mass [GeV];Events / GeV",nmtAxis-1,mtaxis) );
  addHistogram( new TH1F( "mt_Outbtag50"  ,         ";Transverse mass [GeV];Events / GeV",nmtAxis-1,mtaxis) );
  addHistogram( new TH1F( "mt_Outbveto50"  ,         ";Transverse mass [GeV];Events / GeV",nmtAxis-1,mtaxis) );
  addHistogram( new TH1F( "mt_Outbtag80"  ,         ";Transverse mass [GeV];Events / GeV",nmtAxis-1,mtaxis) );
  addHistogram( new TH1F( "mt_Outbveto80"  ,         ";Transverse mass [GeV];Events / GeV",nmtAxis-1,mtaxis) );
  addHistogram( new TH1F( "mt_Outbtag125"  ,         ";Transverse mass [GeV];Events / GeV",nmtAxis-1,mtaxis) );
  addHistogram( new TH1F( "mt_Outbveto125"  ,         ";Transverse mass [GeV];Events / GeV",nmtAxis-1,mtaxis) );
  return true;  
}
bool SmartSelectionMonitor_hzz::declareHistos_jpsipipi(){
 
 addHistogram( new TH1F("Num_J/Psi",";N_{J/#psi};Events",5,0,5));
 addHistogram( new TH1F("M_J/Psi",";M_{J/#psi};Events",50,2.8,3.3));
 addHistogram( new TH1F("pT_J/Psi",";p_{T,J/Psi};Events",120,0,120)); 
 addHistogram( new TH1F("M_J/PsiPicut4.2",";m_{J/psi,Pi};Events",40,3.5,4.3)); 
 addHistogram( new TH1F("M_J/PsiPicut4.25",";m_{J/psi,Pi};Events",40,3.5,4.3));
 addHistogram( new TH1F("M_J/PsiPicut4.3",";m_{J/psi,Pi};Events",40,3.5,4.3)); 
 addHistogram( new TH1F("M_J/PsiPicut4.4",";m_{J/psi,Pi};Events",40,3.5,4.3)); 
 addHistogram( new TH1F("M_J/PsiPicut4.7",";m_{J/psi,Pi};Events",40,3.5,4.3)); 
 addHistogram( new TH1F("M_J/PsiPicut5.0",";m_{J/psi,Pi};Events",40,3.5,4.3)); 
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
