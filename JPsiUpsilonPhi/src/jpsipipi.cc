// -*- C++ -*-
//
// Package:    jpsipipi
// Class:      jpsipipi
// 

//=================================================
// original author:  Jhovanny Andres Mejia        |
//         created:  Monday Aug 28 (2017)         |
//         <jhovanny.andres.mejia.guisao@cern.ch> | 
// modified by    :   Hanwen Wang
//         <allenwang@buaa.edu.cn>
//=================================================

// system include files
#include <memory>


#include "jpsipipi.h"


#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "MagneticField/Engine/interface/MagneticField.h"            
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TH2F.h"
#include "TMath.h"
//
// constants, enums and typedefs
//

  typedef math::Error<3>::type CovarianceMatrix;

//
// static data member definitions
//

//
// constructors and destructor
//

jpsipipi::jpsipipi(const edm::ParameterSet& iConfig)
  :
  dimuon_Label(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("dimuons"))),
  trakCollection_label(consumes<edm::View<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("Trak"))),
  primaryVertices_Label(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertices"))),
   
  isMC_(iConfig.getParameter<bool>("isMC")),


  tree_(0), 

  mumC2(0), mumNHits(0), mumNPHits(0),
  mupC2(0), mupNHits(0), mupNPHits(0),
  mumdxy(0), mupdxy(0), mumdz(0), mupdz(0),
  muon_dca(0),

  mu1soft(0), mu2soft(0), mu1tight(0), mu2tight(0), 
  mu1PF(0), mu2PF(0), mu1loose(0), mu2loose(0),
 
  nJ(0),nPiPair(0),

  J_mass(0), J_px(0), J_py(0), J_pz(0),J_energy(0),

  J_px1(0), J_py1(0), J_pz1(0),
  J_px2(0), J_py2(0), J_pz2(0), 
  J_charge1(0), J_charge2(0),
  J_lxy(0),
  J_lxyErr(0),
  J_vertexchi2(0),
  Pi_dJP(0),
  JPi_lxy(0),
  JPi_lxyErr(0),
  JPiPi_lxy(0),
  JPiPi_lxyErr(0),
  Pi_nhits1(0),
  Pi_npixelhits1(0),
  Pi_nhits2(0),
  Pi_npixelhits2(0),
  Pi_eta1(0),
  Pi_eta2(0),
  Pi_phi1(0),
  Pi_phi2(0),
  Pi_pt1(0),
  Pi_pt2(0),
  Pi_e1(0),
  Pi_e2(0),
  Pi_charge1(0),
  Pi_charge2(0),
  Pi_dxy1(0),
  Pi_dxy2(0),
  Pi_dxyerr1(0),
  Pi_dxyerr2(0),
  Pi_vertexchisq1(0),
  Pi_vertexchisq2(0),
  JPiPi_x(0),
  JPiPi_y(0),
  JPiPi_z(0),
  Pi1_hcalFraction(0),
  Pi2_hcalFraction(0),
  Pi1_vertexNdof(0),
  Pi2_vertexNdof(0),
  Pi1_vertexNchi2(0),
  Pi2_vertexNchi2(0),
  Pi1_lambda(0),
  Pi2_lambda(0),
  Pi1_lambdaError(0),
  Pi2_lambdaError(0),
  Pi1_qoverp(0),
  Pi2_qoverp(0),
  Pi1_qoverpError(0),
  Pi2_qoverpError(0),
  Pi1_validTkFraction(0),
  Pi2_validTkFraction(0),
  Pi1_numberOfMothers(0),
  Pi2_numberOfMothers(0),
  Pi1_numberOfSourceCandidatePtrs(0),
  Pi2_numberOfSourceCandidatePtrs(0),
  Pi1_pdgId(0),
  Pi2_pdgId(0),
  Pi1_numberOfValidHitsOnTrack(0),
  Pi2_numberOfValidHitsOnTrack(0),
  Pi1_innerDetId(0),
  Pi2_innerDetId(0),
  Pi1_innerOk(0),
  Pi2_innerOk(0),
  Pi1_isCaloMuon(0),
  Pi2_isCaloMuon(0),
  Pi1_isConvertedPhoton(0),
  Pi2_isConvertedPhoton(0),
  Pi1_isElectron(0),
  Pi2_isElectron(0),
  Pi1_isMuon (0),
  Pi2_isMuon(0),
  Pi1_isPhoton(0),
  Pi2_isPhoton (0),
  Pi1_isGlobalMuon(0),
  Pi2_isGlobalMuon(0),
  Pi1_isJet(0),
  Pi2_isJet(0),
  Pi1_isLonglived(0),
  Pi2_isLonglived(0),
  Pi1_massConstraint(0),
  Pi2_massConstraint(0),
  J_Prob(0),
JPi_Prob(0),
JPiPi_Prob(0)


 {
   //now do what ever initialization is needed
 }


jpsipipi::~jpsipipi()
{

}


//
// member functions
//

// ------------ method called to for each event  ------------
void jpsipipi::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using std::vector;
  using namespace edm;
  using namespace reco;
  using namespace std;
  
  //*********************************
  // Get event content information
  //*********************************  

  // Kinematic fit
  edm::ESHandle<TransientTrackBuilder> theB; 
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB); 

  edm::Handle<View<pat::PackedCandidate>> thePATTrackHandle;
  iEvent.getByToken(trakCollection_label,thePATTrackHandle);


  edm::Handle< View<pat::Muon> > thePATMuonHandle;
  iEvent.getByToken(dimuon_Label,thePATMuonHandle);
  
 
  //*********************************
  //Now we get the primary vertex 
  //*********************************

  reco::Vertex bestVtx;
  edm::Handle<reco::VertexCollection> primaryVertices_handle;
  iEvent.getByToken(primaryVertices_Label, primaryVertices_handle);

  bestVtx = *(primaryVertices_handle->begin());
  const GlobalPoint pvertex(bestVtx.x(),bestVtx.y(),bestVtx.z());
  //nVtx = primaryVertices_handle->size(); 
  
  //*****************************************
  //Let's begin by looking for J/psi

  //unsigned int nMu_tmp = thePATMuonHandle->size();
 std::vector<FreeTrajectoryState> JpsiFTS;
 std::vector<reco::TransientTrack> muontt1;
 std::vector<reco::TransientTrack> muontt2;

 ParticleMass muon_mass = 0.10565837; //pdg mass
//ParticleMass psi_mass = 3.096916;
 float muon_sigma = muon_mass*1.e-6;
 for(View<pat::Muon>::const_iterator iMuon1 = thePATMuonHandle->begin(); iMuon1 != thePATMuonHandle->end(); ++iMuon1) 
 {  
    //if(!(iMuon1->isGlobalMuon())) continue;
    //if(iMuon1->pt()<1.5) continue;
    //if(!(iMuon1->track())) continue;

    for(View<pat::Muon>::const_iterator iMuon2 = iMuon1+1; iMuon2 != thePATMuonHandle->end(); ++iMuon2) 
	{
	  if(iMuon1==iMuon2) continue;
	  
	  //opposite charge 
	  if( (iMuon1->charge())*(iMuon2->charge()) == 1) continue;
      if (!(iMuon1->isGlobalMuon()) && !(iMuon2->isGlobalMuon()) ) continue;
      //if (!(iMuon1->isTrackerMuon()) || !(iMuon2->isTrackerMuon()) ) continue;
	  TrackRef glbTrackP;	  
	  TrackRef glbTrackM;	  
	  
	  if(iMuon1->charge() == 1){ glbTrackP = iMuon1->track();}
	  if(iMuon1->charge() == -1){ glbTrackM = iMuon1->track();}
	  
	  if(iMuon2->charge() == 1) { glbTrackP = iMuon2->track();}
	  if(iMuon2->charge() == -1){ glbTrackM = iMuon2->track();}
	  
	  if( glbTrackP.isNull() || glbTrackM.isNull() ) 
	    {
	      //std::cout << "continue due to no track ref" << endl;
	      continue;
	    }

	  if(iMuon1->track()->pt()<4.0) continue;
	  if(iMuon2->track()->pt()<4.0) continue;

	  if(!(glbTrackM->quality(reco::TrackBase::highPurity))) continue;
	  if(!(glbTrackP->quality(reco::TrackBase::highPurity))) continue;	 
	  
	  reco::TransientTrack muon1TT((*theB).build(glbTrackP));
	  reco::TransientTrack muon2TT((*theB).build(glbTrackM));

	 // *****  Trajectory states to calculate DCA for the 2 muons *********************
	  FreeTrajectoryState mu1State = muon1TT.impactPointTSCP().theState();
	  FreeTrajectoryState mu2State = muon2TT.impactPointTSCP().theState();

	  if( !muon1TT.impactPointTSCP().isValid() || !muon2TT.impactPointTSCP().isValid() ) continue;

	  // Measure distance between tracks at their closest approach
	  ClosestApproachInRPhi cApp;
	  cApp.calculate(mu1State, mu2State);
	  if( !cApp.status() ) continue;
	  float dca = fabs( cApp.distance() );	  
	  if (dca < 0. || dca > 0.5) continue;
	  //cout<<" closest approach  "<<dca<<endl;

	  // ******  Methods to check to which category of muon candidates a given pat::Muon object belongs ****

	  /*
	 
 //if (iMuon1->isTrackerMuon() || iMuon2->isTrackerMuon())
	  //if (muon::isHighPtMuon(*iMuon1,bestVtx) || muon::isHighPtMuon(*iMuon2,bestVtx))
	  if (muon::isGoodMuon(*iMuon1,muon::TMLastStationAngTight) || muon::isGoodMuon(*iMuon2,muon::TMLastStationAngTight))
	    {
	      cout<<" is category muon  "<<endl;
	    }
	  else
	    {
	      cout<<" it is not category muon  "<<endl;
	    }
	  */
	  
	  // ******   Let's check the vertex and mass ****

	  
	  //The mass of a muon and the insignificant mass sigma 
	  //to avoid singularities in the covariance matrix.
	  
	  //float psi_sigma = psi_mass*1.e-6;
	  
	  //Creating a KinematicParticleFactory
	  KinematicParticleFactoryFromTransientTrack pFactory;
	  
	  //initial chi2 and ndf before kinematic fits.
	  float chi = 0.;
	  float ndf = 0.;
	  vector<RefCountedKinematicParticle> muonParticles;
	  try {
	    muonParticles.push_back(pFactory.particle(muon1TT,muon_mass,chi,ndf,muon_sigma));
	    muonParticles.push_back(pFactory.particle(muon2TT,muon_mass,chi,ndf,muon_sigma));
	  }
	  catch(...) { 
	    std::cout<<" Exception caught ... continuing 1 "<<std::endl; 
	    continue;
	  }
	  
	  KinematicParticleVertexFitter fitter;   
	  
	  RefCountedKinematicTree psiVertexFitTree;
	  try {
	    psiVertexFitTree = fitter.fit(muonParticles); 
	  }
	  catch (...) { 
	    std::cout<<" Exception caught ... continuing 2 "<<std::endl; 
	    continue;
	  }
	  
	  if (!psiVertexFitTree->isValid()) 
	  {
	      //std::cout << "caught an exception in the psi vertex fit" << std::endl;
	      continue; 
	  }
	  
	  psiVertexFitTree->movePointerToTheTop();
	  
	  RefCountedKinematicParticle psi_vFit_noMC = psiVertexFitTree->currentParticle();
	  RefCountedKinematicVertex psi_vFit_vertex_noMC = psiVertexFitTree->currentDecayVertex();
	  
	  if( psi_vFit_vertex_noMC->chiSquared() < 0 )
	    {
	      //std::cout << "negative chisq from psi fit" << endl;
	      continue;
	    }
	  
	  //some loose cuts go here
	  
	  //if(psi_vFit_vertex_noMC->chiSquared()>20.) continue;
	  double J_Prob_tmp   = TMath::Prob(psi_vFit_vertex_noMC->chiSquared(),(int)psi_vFit_vertex_noMC->degreesOfFreedom());
	   
	   if(J_Prob_tmp<0.01)
	     {
	       continue;
	     }

	  if(psi_vFit_noMC->currentState().mass()<2.92 || psi_vFit_noMC->currentState().mass()>3.25) continue;
	  float dx = psi_vFit_vertex_noMC->position().x()-bestVtx.x();
	  float dy = psi_vFit_vertex_noMC->position().y()-bestVtx.y();
	  float J_dxy = TMath::Sqrt(dx*dx+dy*dy);
	  float J_dxyerr = psi_vFit_vertex_noMC->error().rerr(pvertex);
      
	  if (J_dxy/J_dxyerr<5.0) continue;
	  muontt1.push_back(muon1TT);
	  muontt2.push_back(muon2TT);
	  J_Prob->push_back(J_Prob_tmp);
	  J_vertexchi2->push_back(psi_vFit_vertex_noMC->chiSquared());
      J_lxy->push_back(J_dxy);
	  J_lxyErr->push_back(J_dxyerr);
	  //fill variables?iMuon1->track()->pt()
	  JpsiFTS.push_back(psi_vFit_noMC->initialState().freeTrajectoryState());

	  J_mass->push_back( psi_vFit_noMC->currentState().mass() );
	  J_px->push_back( psi_vFit_noMC->currentState().globalMomentum().x() );
	  J_py->push_back( psi_vFit_noMC->currentState().globalMomentum().y() );
	  J_pz->push_back( psi_vFit_noMC->currentState().globalMomentum().z() );
	  
	  J_px1->push_back(iMuon1->track()->px());
	  J_py1->push_back(iMuon1->track()->py());
	  J_pz1->push_back(iMuon1->track()->pz());
	  J_charge1->push_back(iMuon1->charge());
	  
	  J_px2->push_back(iMuon2->track()->px());
	  J_py2->push_back(iMuon2->track()->py());
	  J_pz2->push_back(iMuon2->track()->pz());
	  J_charge2->push_back(iMuon2->charge());
	  
	   // ************
	  
	  mu1soft->push_back(iMuon1->isSoftMuon(bestVtx) );
	  mu2soft->push_back(iMuon2->isSoftMuon(bestVtx) );
	  mu1tight->push_back(iMuon1->isTightMuon(bestVtx) );
	  mu2tight->push_back(iMuon2->isTightMuon(bestVtx) );
	  mu1PF->push_back(iMuon1->isPFMuon());
	  mu2PF->push_back(iMuon2->isPFMuon());
	  mu1loose->push_back(muon::isLooseMuon(*iMuon1));
	  mu2loose->push_back(muon::isLooseMuon(*iMuon2));
	  
	  mumC2->push_back( glbTrackP->normalizedChi2() );
	  //mumAngT->push_back( muon::isGoodMuon(*iMuon1,muon::TMLastStationAngTight) ); // 
	  mumNHits->push_back( glbTrackP->numberOfValidHits() );
	  mumNPHits->push_back( glbTrackP->hitPattern().numberOfValidPixelHits() );	       
	  mupC2->push_back( glbTrackM->normalizedChi2() );
	  //mupAngT->push_back( muon::isGoodMuon(*iMuon2,muon::TMLastStationAngTight) );  // 
	  mupNHits->push_back( glbTrackM->numberOfValidHits() );
	  mupNPHits->push_back( glbTrackM->hitPattern().numberOfValidPixelHits() );
	  mumdxy->push_back(glbTrackP->dxy(bestVtx.position()) );
	  mupdxy->push_back(glbTrackM->dxy(bestVtx.position()) );
	  mumdz->push_back(glbTrackP->dz(bestVtx.position()) );
	  mupdz->push_back(glbTrackM->dz(bestVtx.position()) );
	  muon_dca->push_back(dca);          	  
	  
	  nJ++;	       
	  muonParticles.clear();
	  
	}
}


for(unsigned int i=0; i<JpsiFTS.size(); i++)
{
	reco::TransientTrack JpsiTT((*theB).build(JpsiFTS.at(i)));
	for(View<pat::PackedCandidate>::const_iterator iTrack1= thePATTrackHandle->begin(); iTrack1 != thePATTrackHandle->end();++iTrack1)
	{
		
		
  	    //if(iTrack1->pt()<0.8)continue;
  	    if(iTrack1->eta()>2.4||iTrack1->eta()<-2.4)continue;
  	    if(!(iTrack1->bestTrack())) continue;
  	    if(iTrack1->charge() == 0) continue; //NO neutral objects
  	    //if(fabs(iTrack1->pdgId()!= 211)) continue; //Due to the lack of the particle ID all the tracks for cms are pions(ID == 211)
  	    
  	    if(iTrack1->vertexNormalizedChi2()>10) continue;
  	    if(!(iTrack1->trackHighPurity())) continue;
  	    if(iTrack1->numberOfHits()<5) continue;
  	    if(iTrack1->numberOfPixelHits()<1) continue;
  	    if(iTrack1->dxy(bestVtx.position())/iTrack1->dxyError() < 2.0) continue;
        if ( IsTheSame(*iTrack1, muontt1.at(i).track()) || IsTheSame(*iTrack1,muontt2.at(i).track()) ) continue;
  	    reco::TransientTrack track1TT((*theB).build(iTrack1->pseudoTrack()));
  	    //FreeTrajectoryState pi_trajectory = track1TT.impactPointTSCP().theState();
  	    //ClosestApproachInRPhi JpsiPi;

  	    
  	    //JpsiPi.calculate(JpsiFTS.at(i), pi_trajectory);
  	    //if( !JpsiPi.status() ) continue;
	    float djp = 0; //fabs( JpsiPi.distance() );	  
	    //if (djp < 0. || djp > 0.03) continue;
	   
	    //begin vertex fit of Jpsi and pi1
       //ParticleMass Jpsi_mass = 3.0969;
        ParticleMass Pion_mass = 0.13957061;
	    //float Jpsi_sigma = 0.000006 ;
	    float Pion_sigma = Pion_mass*1.e-6;
	    //float psi_sigma = psi_mass*1.e-6;
	    float chi = 0.;
	    float ndf = 0.;
	    //Creating a KinematicParticleFactory
	    KinematicParticleFactoryFromTransientTrack pFactory;
	    

	    vector<RefCountedKinematicParticle> JpsiPi_fit;
	    try {
	      //JpsiPi_fit.push_back(pFactory.particle(JpsiTT,Jpsi_mass,J_vertexFitChi2->at(i),J_vertexFitNdf->at(i),Jpsi_sigma));
	     // JpsiPi_fit.push_back(pFactory.particle(track1TT,Pion_mass,iTrack1->vertexChi2(),iTrack1->vertexNdof(),Pion_sigma));
	      JpsiPi_fit.push_back(pFactory.particle(muontt1.at(i),muon_mass,chi,ndf,muon_sigma));
	      JpsiPi_fit.push_back(pFactory.particle(muontt2.at(i),muon_mass,chi,ndf,muon_sigma));
	      JpsiPi_fit.push_back(pFactory.particle(track1TT,     Pion_mass,chi,ndf,Pion_sigma));
	    }
	    catch(...) { 
	      std::cout<<" Exception caught ... continuing 1 "<<std::endl; 
	      continue;
	    }
	    
	    KinematicParticleVertexFitter fitter;   
	    
	    RefCountedKinematicTree psiVertexFitTree;
	    try {
	      psiVertexFitTree = fitter.fit(JpsiPi_fit); 
	    }
	    catch (...) { 
	      std::cout<<" Exception caught ... continuing 2 "<<std::endl; 
	      continue;
	    }
	    
	    if (!psiVertexFitTree->isValid()) 
	    {
	        //std::cout << "caught an exception in the psi vertex fit" << std::endl;
	        continue; 
	    }
	    
	    psiVertexFitTree->movePointerToTheTop();
	    
	    RefCountedKinematicParticle psi_vFit_noMC = psiVertexFitTree->currentParticle();
	    RefCountedKinematicVertex psi_vFit_vertex_noMC = psiVertexFitTree->currentDecayVertex();
	    
	    if( psi_vFit_vertex_noMC->chiSquared() < 0 )
	      {
	        //std::cout << "negative chisq from psi fit" << endl;
	        continue;
	      }
	    //if(psi_vFit_vertex_noMC->chiSquared()>15.) continue;
	     double Prob_tmp  = TMath::Prob(psi_vFit_vertex_noMC->chiSquared(),(int)psi_vFit_vertex_noMC->degreesOfFreedom());
		   if(Prob_tmp<0.01)
		     {
		       continue;
		     }
	    float dx1 = psi_vFit_vertex_noMC->position().x()-bestVtx.x();
	    float dy1 = psi_vFit_vertex_noMC->position().y()-bestVtx.y();
	    float JpsiPi_dxy = TMath::Sqrt(dx1*dx1+dy1*dy1);
	    float JpsiPi_dxyerr = psi_vFit_vertex_noMC->error().rerr(pvertex);
	    JpsiPi_fit.clear();
	    //if (JpsiPi_dxy/JpsiPi_dxyerr<3) continue;
	    for(View<pat::PackedCandidate>::const_iterator iTrack2= iTrack1+1; iTrack2 != thePATTrackHandle->end();++iTrack2)
	    {
            if((iTrack1->charge())*(iTrack2->charge())==1) continue;
            
            //if(iTrack2->pt()<0.8)continue;
  	        //if(iTrack2->eta()>2.4||iTrack2->eta()<-2.4)continue;
  	        if(!(iTrack2->bestTrack())) continue;
  	        if(iTrack2->charge() == 0) continue; //NO neutral objects
  	        //if(fabs(iTrack2->pdgId()!= 211)) continue; //Due to the lack of the particle ID all the tracks for cms are pions(ID == 211)
  	        if(!(iTrack2->trackHighPurity())) continue;
  	        if(iTrack2->numberOfHits()<5) continue;
  	        if(iTrack2->numberOfPixelHits()<1) continue;
  	        if(iTrack2->vertexNormalizedChi2()>10) continue;
  	        if(iTrack2->dxy(bestVtx.position())/iTrack2->dxyError() < 1.0) continue;
  	        if ( IsTheSame(*iTrack2, muontt1.at(i).track()) || IsTheSame(*iTrack2,muontt2.at(i).track()) ) continue;
            reco::TransientTrack track2TT((*theB).build(iTrack2->pseudoTrack()));
  	        //begin vertex fit of Jpsi and pi1
            //ParticleMass Jpsi_mass = 3.0969;
            //ParticleMass Pion_mass = 0.13957061;
	        //float Jpsi_sigma = 0.000006 ;
	        //float Pion_sigma = 0.00000024;
	        //float psi_sigma = psi_mass*1.e-6;
	        
	        //Creating a KinematicParticleFactory
	        KinematicParticleFactoryFromTransientTrack pFactory2;
	      
    
	        vector<RefCountedKinematicParticle> JpsiPiPi_fit;
	        try {
	          //JpsiPi_fit.push_back(pFactory.particle(JpsiTT,Jpsi_mass,J_vertexFitChi2->at(i),J_vertexFitNdf->at(i),Jpsi_sigma));
	         // JpsiPi_fit.push_back(pFactory.particle(track1TT,Pion_mass,iTrack1->vertexChi2(),iTrack1->vertexNdof(),Pion_sigma));
	          JpsiPiPi_fit.push_back(pFactory.particle(muontt1.at(i),muon_mass,chi,ndf,muon_sigma));
	          JpsiPiPi_fit.push_back(pFactory.particle(muontt2.at(i),muon_mass,chi,ndf,muon_sigma));
	          JpsiPiPi_fit.push_back(pFactory.particle(track1TT,     Pion_mass,chi,ndf,Pion_sigma));
	          JpsiPiPi_fit.push_back(pFactory.particle(track2TT,     Pion_mass,chi,ndf,Pion_sigma));
	        }
	        catch(...) { 
	          std::cout<<" Exception caught ... continuing 1 "<<std::endl; 
	          continue;
	        }
	        
	        KinematicParticleVertexFitter fitter2;   
	        
	        RefCountedKinematicTree psiVertexFitTree2;
	        try {
	          psiVertexFitTree2 = fitter2.fit(JpsiPiPi_fit); 
	        }
	        catch (...) { 
	          std::cout<<" Exception caught ... continuing 2 "<<std::endl; 
	          continue;
	        }
	        
	        if (!psiVertexFitTree2->isValid()) 
	        {
	            //std::cout << "caught an exception in the psi vertex fit" << std::endl;
	            continue; 
	        }
	        
	        psiVertexFitTree2->movePointerToTheTop();
	        
	        RefCountedKinematicParticle psi_vFit_noMC2 = psiVertexFitTree2->currentParticle();
	        RefCountedKinematicVertex psi_vFit_vertex_noMC2 = psiVertexFitTree2->currentDecayVertex();
	        
	        if( psi_vFit_vertex_noMC2->chiSquared() < 0 )
	          {
	            //std::cout << "negative chisq from psi fit" << endl;
	            continue;
	          }
	        //if(psi_vFit_vertex_noMC2->chiSquared() > psi_vFit_vertex_noMC->chiSquared() + 6.0 ) continue;
	        double Prob_tmp2  = TMath::Prob(psi_vFit_vertex_noMC2->chiSquared(),(int)psi_vFit_vertex_noMC2->degreesOfFreedom());
		   if(Prob_tmp2<0.01)
		     {
		       continue;
		     }
	        float px,py,rx,ry;
	        px = psi_vFit_noMC2->currentState().globalMomentum().x();
	        py = psi_vFit_noMC2->currentState().globalMomentum().y();
	        //pz = psi_vFit_noMC2->currentState().globalMomentum().z();
	        rx = psi_vFit_vertex_noMC2->position().x()-bestVtx.x();
	        ry = psi_vFit_vertex_noMC2->position().y()-bestVtx.y();
	        //rz = psi_vFit_vertex_noMC2->position().z()-bestVtx.z();
	        //float cosine = (px*rx+py*ry+pz*rz)/((px*px+py*py+pz*pz)*(rx*rx+ry*ry+rz*rz));
	        float cosine = (px*rx+py*ry)/((px*px+py*py)*(rx*rx+ry*ry));
	        if (cosine < 0.9 ) continue;
            float JpsiPiPi_dxy = TMath::Sqrt(rx*rx+ry*ry);
	        float JpsiPiPi_dxyerr = psi_vFit_vertex_noMC2->error().rerr(pvertex);
	        //if (JpsiPiPi_dxy/JpsiPiPi_dxyerr<2) continue;
            JPi_Prob->push_back(Prob_tmp);
            JPiPi_Prob->push_back(Prob_tmp2);
	        Pi_nhits1->push_back(iTrack1->numberOfHits());
            Pi_npixelhits1->push_back(iTrack1->numberOfPixelHits());
            Pi_nhits2->push_back(iTrack2->numberOfHits());
            Pi_npixelhits2->push_back(iTrack2->numberOfPixelHits());
            Pi_eta1->push_back(iTrack1->eta());
            Pi_eta2->push_back(iTrack2->eta());
            Pi_phi1->push_back(iTrack1->phi());
            Pi_phi2->push_back(iTrack2->phi());
            Pi_pt1->push_back(iTrack1->pt());
            Pi_pt2->push_back(iTrack2->pt());
            Pi_e1->push_back(iTrack1->energy());
            Pi_e2->push_back(iTrack2->energy());
            Pi_charge1->push_back(iTrack1->charge());
            Pi_charge2->push_back(iTrack2->charge());
            Pi_dxy1->push_back(iTrack1->dxy(bestVtx.position()));
            Pi_dxy2->push_back(iTrack2->dxy(bestVtx.position()));
            Pi_dxyerr1->push_back(iTrack1->dxyError());
            Pi_dxyerr2->push_back(iTrack2->dxyError());
            Pi_vertexchisq1->push_back(psi_vFit_vertex_noMC->chiSquared());
            Pi_vertexchisq2->push_back(psi_vFit_vertex_noMC2->chiSquared());
            Pi_dJP->push_back(djp);
            JPi_lxy->push_back(JpsiPi_dxy);
            JPi_lxyErr->push_back(JpsiPi_dxyerr);
            JPiPi_lxy->push_back(JpsiPiPi_dxy);
            JPiPi_lxyErr->push_back(JpsiPiPi_dxyerr);
            JPiPi_x->push_back( psi_vFit_vertex_noMC2->position().x()-bestVtx.x());
            JPiPi_y->push_back( psi_vFit_vertex_noMC2->position().y()-bestVtx.y());
            JPiPi_z->push_back( psi_vFit_vertex_noMC2->position().z()-bestVtx.z());

            Pi1_hcalFraction->push_back(iTrack1->hcalFraction());
            Pi2_hcalFraction->push_back(iTrack2->hcalFraction());
            Pi1_vertexNdof->push_back(iTrack1->vertexNdof());
            Pi2_vertexNdof->push_back(iTrack2->vertexNdof());
            Pi1_vertexNchi2->push_back(iTrack1->vertexNormalizedChi2());
            Pi2_vertexNchi2->push_back(iTrack2->vertexNormalizedChi2());
            Pi1_lambda->push_back(iTrack1->pseudoTrack().lambda());
            Pi2_lambda->push_back(iTrack2->pseudoTrack().lambda());
            Pi1_lambdaError->push_back(iTrack1->pseudoTrack().lambdaError());
            Pi2_lambdaError->push_back(iTrack2->pseudoTrack().lambdaError());
            Pi1_qoverp->push_back(iTrack1->pseudoTrack().qoverp());
            Pi2_qoverp->push_back(iTrack2->pseudoTrack().qoverp());
            Pi1_qoverpError->push_back(iTrack1->pseudoTrack().qoverpError());
            Pi2_qoverpError->push_back(iTrack2->pseudoTrack().qoverpError());
            Pi1_validTkFraction->push_back(iTrack1->pseudoTrack().validFraction());
            Pi2_validTkFraction->push_back(iTrack2->pseudoTrack().validFraction());
            Pi1_numberOfMothers->push_back(iTrack1->numberOfMothers());
            Pi2_numberOfMothers->push_back(iTrack2->numberOfMothers());
            Pi1_numberOfSourceCandidatePtrs->push_back(iTrack1->numberOfSourceCandidatePtrs());
            Pi2_numberOfSourceCandidatePtrs->push_back(iTrack2->numberOfSourceCandidatePtrs());
            Pi1_pdgId->push_back(iTrack1->pdgId());
            Pi2_pdgId->push_back(iTrack2->pdgId());
            Pi1_numberOfValidHitsOnTrack->push_back(iTrack1->pseudoTrack().found());
            Pi2_numberOfValidHitsOnTrack->push_back(iTrack2->pseudoTrack().found());
            Pi1_innerDetId->push_back(iTrack1->pseudoTrack().innerDetId());
            Pi2_innerDetId->push_back(iTrack2->pseudoTrack().innerDetId());
            Pi1_innerOk->push_back(iTrack1->pseudoTrack().innerOk());
            Pi2_innerOk->push_back(iTrack2->pseudoTrack().innerOk());
            Pi1_isCaloMuon->push_back(iTrack1->isCaloMuon());
            Pi2_isCaloMuon->push_back(iTrack2->isCaloMuon());
            Pi1_isConvertedPhoton->push_back(iTrack1->isConvertedPhoton());
            Pi2_isConvertedPhoton->push_back(iTrack2->isConvertedPhoton());
            Pi1_isElectron->push_back(iTrack1->isElectron());
            Pi2_isElectron->push_back(iTrack2->isElectron());
            Pi1_isMuon ->push_back(iTrack1->isMuon());
            Pi2_isMuon->push_back(iTrack2->isMuon());
            Pi1_isPhoton->push_back(iTrack1->isPhoton());
            Pi2_isPhoton ->push_back(iTrack2->isPhoton());
            Pi1_isGlobalMuon->push_back(iTrack1->isGlobalMuon());
            Pi2_isGlobalMuon->push_back(iTrack2->isGlobalMuon());
            Pi1_isJet->push_back(iTrack1->isJet());
            Pi2_isJet->push_back(iTrack2->isJet());
            Pi1_isLonglived->push_back(iTrack1->longLived());
            Pi2_isLonglived->push_back(iTrack2->longLived());
            Pi1_massConstraint->push_back(iTrack1->massConstraint());
            Pi2_massConstraint->push_back(iTrack2->massConstraint());

	        nPiPair++;
	        JpsiPiPi_fit.clear();

	    }
	}
}


  //for (reco::PackedCandidate::const_iterator iTrack1= thePATTrackHandle->begin();  iTrack1 != thePATTrackHandle->end(); ++iTrack1){
  	
  
  
   if (nJ > 0 && nPiPair>0)
   { tree_->Fill();

    //std::cout << "filling tree" << endl;
   
    }

   
   nJ = 0; 
   nPiPair =0;
   

   J_mass->clear(); J_px->clear();   J_py->clear();  J_pz->clear();  
   J_px1->clear();  J_py1->clear();  J_pz1->clear(), J_charge1->clear();
   J_px2->clear();  J_py2->clear();  J_pz2->clear(), J_charge2->clear();

   
   mupC2->clear();
   mumC2->clear();
   mupNHits->clear(); mupNPHits->clear();
   mumNHits->clear(); mumNPHits->clear();
   mumdxy->clear(); mupdxy->clear(); mumdz->clear(); mupdz->clear(); muon_dca->clear();

   mu1soft->clear(); mu2soft->clear(); mu1tight->clear(); mu2tight->clear();
   mu1PF->clear(); mu2PF->clear(); mu1loose->clear(); mu2loose->clear(); 
   JpsiFTS.clear();
   muontt1.clear();
   muontt2.clear();
   J_lxy->clear();
   J_lxyErr->clear();
   J_vertexchi2->clear();
   Pi_dJP->clear();
   JPi_lxy->clear();
   JPi_lxyErr->clear();
   JPiPi_lxy->clear();
   JPiPi_lxyErr->clear();
   Pi_nhits1->clear();
   Pi_npixelhits1->clear();
   Pi_nhits2->clear();
   Pi_npixelhits2->clear();
   Pi_eta1->clear();
   Pi_eta2->clear();
   Pi_phi1->clear();
   Pi_phi2->clear();
   Pi_pt1->clear();
   Pi_pt2->clear();
   Pi_e1->clear();
   Pi_e2->clear();
   Pi_charge1->clear();
   Pi_charge2->clear();
   Pi_dxy1->clear();
   Pi_dxy2->clear();
   Pi_dxyerr1->clear();
   Pi_dxyerr2->clear();
   Pi_vertexchisq1->clear();
   Pi_vertexchisq2->clear();
   JPiPi_x->clear();
   JPiPi_y->clear();
   JPiPi_z->clear();
   Pi1_hcalFraction->clear();
   Pi2_hcalFraction->clear();
   Pi1_vertexNdof->clear();
   Pi2_vertexNdof->clear();
   Pi1_vertexNchi2->clear();
   Pi2_vertexNchi2->clear();
   Pi1_lambda->clear();
   Pi2_lambda->clear();
   Pi1_lambdaError->clear();
   Pi2_lambdaError->clear();
   Pi1_qoverp->clear();
   Pi2_qoverp->clear();
   Pi1_qoverpError->clear();
   Pi2_qoverpError->clear();
   Pi1_validTkFraction->clear();
   Pi2_validTkFraction->clear();
   Pi1_numberOfMothers->clear();
   Pi2_numberOfMothers->clear();
   Pi1_numberOfSourceCandidatePtrs->clear();
   Pi2_numberOfSourceCandidatePtrs->clear();
   Pi1_pdgId->clear();
   Pi2_pdgId->clear();
   Pi1_numberOfValidHitsOnTrack->clear();
   Pi2_numberOfValidHitsOnTrack->clear();
   Pi1_innerDetId->clear();
   Pi2_innerDetId->clear();
   Pi1_innerOk->clear();
   Pi2_innerOk->clear();
   Pi1_isCaloMuon->clear();
   Pi2_isCaloMuon->clear();
   Pi1_isConvertedPhoton->clear();
   Pi2_isConvertedPhoton->clear();
   Pi1_isElectron->clear();
   Pi2_isElectron->clear();
   Pi1_isMuon ->clear();
   Pi2_isMuon->clear();
   Pi1_isPhoton->clear();
   Pi2_isPhoton ->clear();
   Pi1_isGlobalMuon->clear();
   Pi2_isGlobalMuon->clear();
   Pi1_isJet->clear();
   Pi2_isJet->clear();
   Pi1_isLonglived->clear();
   Pi2_isLonglived->clear();
   Pi1_massConstraint->clear();
   Pi2_massConstraint->clear();
   J_Prob->clear();
   JPi_Prob->clear();
   JPiPi_Prob->clear();

}
bool jpsipipi::IsTheSame(const pat::GenericParticle& tk, const reco::Track  mu){
  double DeltaEta = fabs(mu.eta()-tk.eta());
  double DeltaP   = fabs(mu.p()-tk.p());
  if (DeltaEta < 0.02 && DeltaP < 0.02) return true;
  return false;
}

// ------------ method called once each job just before starting event loop  ------------

void 
jpsipipi::beginJob()
{

  std::cout << "Beginning analyzer job with value of isMC= " << isMC_ << std::endl;

  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("ntuple"," J/psi ntuple");

  tree_->Branch("nJ",&nJ,"nJ/i"); 
  tree_->Branch("nPiPair",&nPiPair,"nPiPair/i"); 
  tree_->Branch("J_mass", &J_mass);
  tree_->Branch("J_px", &J_px);
  tree_->Branch("J_py", &J_py);
  tree_->Branch("J_pz", &J_pz);

  tree_->Branch("J_px1", &J_px1);
  tree_->Branch("J_py1", &J_py1);
  tree_->Branch("J_pz1", &J_pz1);
  tree_->Branch("J_charge1", &J_charge1);

  tree_->Branch("J_px2", &J_px2);
  tree_->Branch("J_py2", &J_py2);
  tree_->Branch("J_pz2", &J_pz2);
  tree_->Branch("J_charge2", &J_charge2);

 
  tree_->Branch("mumC2",&mumC2);  
  tree_->Branch("mumNHits",&mumNHits);
  tree_->Branch("mumNPHits",&mumNPHits);
  tree_->Branch("mupC2",&mupC2);  
  tree_->Branch("mupNHits",&mupNHits);
  tree_->Branch("mupNPHits",&mupNPHits);
  tree_->Branch("mumdxy",&mumdxy);
  tree_->Branch("mupdxy",&mupdxy);
  tree_->Branch("mumdz",&mumdz);
  tree_->Branch("mupdz",&mupdz);
  tree_->Branch("muon_dca",&muon_dca);

  tree_->Branch("mu1soft",&mu1soft);
  tree_->Branch("mu2soft",&mu2soft);
  tree_->Branch("mu1tight",&mu1tight);
  tree_->Branch("mu2tight",&mu2tight);
  tree_->Branch("mu1PF",&mu1PF);
  tree_->Branch("mu2PF",&mu2PF);
  tree_->Branch("mu1loose",&mu1loose);
  tree_->Branch("mu2loose",&mu2loose);

  tree_->Branch("J_lxy",&J_lxy);
  tree_->Branch("J_lxyErr",&J_lxyErr);
  tree_->Branch("J_vertexchi2",&J_vertexchi2);
  tree_->Branch("Pi_dJP",&Pi_dJP);
  tree_->Branch("JPi_lxy",&JPi_lxy);
  tree_->Branch("JPi_lxyErr",&JPi_lxyErr);
  tree_->Branch("JPiPi_lxy",&JPiPi_lxy);
  tree_->Branch("JPiPi_lxyErr",&JPiPi_lxyErr);
  tree_->Branch("Pi_nhits1",&Pi_nhits1);
  tree_->Branch("Pi_npixelhits1",&Pi_npixelhits1);
  tree_->Branch("Pi_nhits2",&Pi_nhits2);
  tree_->Branch("Pi_npixelhits2",&Pi_npixelhits2);
  tree_->Branch("Pi_eta1",&Pi_eta1);
  tree_->Branch("Pi_eta2",&Pi_eta2);
  tree_->Branch("Pi_phi1",&Pi_phi1);
  tree_->Branch("Pi_phi2",&Pi_phi2);
  tree_->Branch("Pi_pt1",&Pi_pt1);
  tree_->Branch("Pi_pt2",&Pi_pt2);
  tree_->Branch("Pi_e1",&Pi_e1);
  tree_->Branch("Pi_e2",&Pi_e2);
  tree_->Branch("Pi_charge1",&Pi_charge1);
  tree_->Branch("Pi_charge2",&Pi_charge2);
  tree_->Branch("Pi_dxy1",&Pi_dxy1);
  tree_->Branch("Pi_dxy2",&Pi_dxy2);
  tree_->Branch("Pi_dxyerr1",&Pi_dxyerr1);
  tree_->Branch("Pi_dxyerr2",&Pi_dxyerr2);
  tree_->Branch("Pi_vertexchisq1",&Pi_vertexchisq1);
  tree_->Branch("Pi_vertexchisq2",&Pi_vertexchisq2);
  tree_->Branch("JPiPi_x",&JPiPi_x);
  tree_->Branch("JPiPi_y",&JPiPi_y);
  tree_->Branch("JPiPi_z",&JPiPi_z);
  tree_->Branch("Pi1_hcalFraction",&Pi1_hcalFraction);
  tree_->Branch("Pi2_hcalFraction",&Pi2_hcalFraction);
  tree_->Branch("Pi1_vertexNdof",&Pi1_vertexNdof);
  tree_->Branch("Pi2_vertexNdof",&Pi2_vertexNdof);
  tree_->Branch("Pi1_vertexNchi2",&Pi1_vertexNchi2);
  tree_->Branch("Pi2_vertexNchi2",&Pi2_vertexNchi2);
  tree_->Branch("Pi1_lambda",&Pi1_lambda);
  tree_->Branch("Pi2_lambda",&Pi2_lambda);
  tree_->Branch("Pi1_lambdaError",&Pi1_lambdaError);
  tree_->Branch("Pi2_lambdaError",&Pi2_lambdaError);
  tree_->Branch("Pi1_qoverp",&Pi1_qoverp);
  tree_->Branch("Pi2_qoverp",&Pi2_qoverp);
  tree_->Branch("Pi1_qoverpError",&Pi1_qoverpError);
  tree_->Branch("Pi2_qoverpError",&Pi2_qoverpError);
  tree_->Branch("Pi1_validTkFraction",&Pi1_validTkFraction);
  tree_->Branch("Pi2_validTkFraction",&Pi2_validTkFraction);
  tree_->Branch("Pi1_numberOfMothers",&Pi1_numberOfMothers);
  tree_->Branch("Pi2_numberOfMothers",&Pi2_numberOfMothers);
  tree_->Branch("Pi1_numberOfSourceCandidatePtrs",&Pi1_numberOfSourceCandidatePtrs);
  tree_->Branch("Pi2_numberOfSourceCandidatePtrs",&Pi2_numberOfSourceCandidatePtrs);
  tree_->Branch("Pi1_pdgId",&Pi1_pdgId);
  tree_->Branch("Pi2_pdgId",&Pi2_pdgId);
  tree_->Branch("Pi1_numberOfValidHitsOnTrack",&Pi1_numberOfValidHitsOnTrack);
  tree_->Branch("Pi2_numberOfValidHitsOnTrack",&Pi2_numberOfValidHitsOnTrack);
  tree_->Branch("Pi1_innerDetId",&Pi1_innerDetId);
  tree_->Branch("Pi2_innerDetId",&Pi2_innerDetId);
  tree_->Branch("Pi1_innerOk",&Pi1_innerOk);
  tree_->Branch("Pi2_innerOk",&Pi2_innerOk);
  tree_->Branch("Pi1_isCaloMuon",&Pi1_isCaloMuon);
  tree_->Branch("Pi2_isCaloMuon",&Pi2_isCaloMuon);
  tree_->Branch("Pi1_isConvertedPhoton",&Pi1_isConvertedPhoton);
  tree_->Branch("Pi2_isConvertedPhoton",&Pi2_isConvertedPhoton);
  tree_->Branch("Pi1_isElectron",&Pi1_isElectron);
  tree_->Branch("Pi2_isElectron",&Pi2_isElectron);
  tree_->Branch("Pi1_isMuon ",&Pi1_isMuon );
  tree_->Branch("Pi2_isMuon",&Pi2_isMuon);
  tree_->Branch("Pi1_isPhoton",&Pi1_isPhoton);
  tree_->Branch("Pi2_isPhoton ",&Pi2_isPhoton );
  tree_->Branch("Pi1_isGlobalMuon",&Pi1_isGlobalMuon);
  tree_->Branch("Pi2_isGlobalMuon",&Pi2_isGlobalMuon);
  tree_->Branch("Pi1_isJet",&Pi1_isJet);
  tree_->Branch("Pi2_isJet",&Pi2_isJet);
  tree_->Branch("Pi1_isLonglived",&Pi1_isLonglived);
  tree_->Branch("Pi2_isLonglived",&Pi2_isLonglived);
  tree_->Branch("Pi1_massConstraint",&Pi1_massConstraint);
  tree_->Branch("Pi2_massConstraint",&Pi2_massConstraint);
  tree_->Branch("J_Prob",&J_Prob);
  tree_->Branch("JPi_Prob",&JPi_Prob);
  tree_->Branch("JPiPi_Prob",&JPiPi_Prob);

}




// ------------ method called once each job just after ending the event loop  ------------
void jpsipipi::endJob() {
  tree_->GetDirectory()->cd();
  tree_->Write();
}

//define this as a plug-in
DEFINE_FWK_MODULE(jpsipipi);