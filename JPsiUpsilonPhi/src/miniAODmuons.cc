// -*- C++ -*-
//
// Package:    miniAODmuons
// Class:      miniAODmuons
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


#include "miniAODmuons.h"

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

miniAODmuons::miniAODmuons(const edm::ParameterSet& iConfig)
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
 
  
  nU(0),
  nJ(0),
  nPhi(0),
  
 
  J_mass(0), J_px(0), J_py(0), J_pz(0),J_energy(0),J_eta(0),J_phi(0),J_pt(0),

  J_px1(0), J_py1(0), J_pz1(0),
  J_px2(0), J_py2(0), J_pz2(0), 
  J_charge1(0), J_charge2(0),

  U_mumC2(0), U_mumNHits(0), U_mumNPHits(0),
  U_mupC2(0), U_mupNHits(0), U_mupNPHits(0),
  U_mumdxy(0), U_mupdxy(0), U_mumdz(0), U_mupdz(0),
  U_muon_dca(0),

  U_mu1soft(0), U_mu2soft(0), U_mu1tight(0), U_mu2tight(0), 
  U_mu1PF(0), U_mu2PF(0), U_mu1loose(0), U_mu2loose(0),
  U_mass(0), U_px(0), U_py(0), U_pz(0),U_energy(0),U_eta(0),U_phi(0),U_pt(0),

  U_px1(0), U_py1(0), U_pz1(0),
  U_px2(0), U_py2(0), U_pz2(0), 
  U_charge1(0), U_charge2(0),
  Phi_mass(0),Phi_px(0),Phi_py(0),Phi_pz(0),Phi_energy(0),Phi_eta(0),Phi_phi(0),Phi_pt(0)


{
   //now do what ever initialization is needed
}


miniAODmuons::~miniAODmuons()
{

}


//
// member functions
//

// ------------ method called to for each event  ------------
void miniAODmuons::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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

  //nVtx = primaryVertices_handle->size(); 
  
  //*****************************************
  //Let's begin by looking for J/psi

  //unsigned int nMu_tmp = thePATMuonHandle->size();
 
 for(View<pat::Muon>::const_iterator iMuon1 = thePATMuonHandle->begin(); iMuon1 != thePATMuonHandle->end(); ++iMuon1) 
    {
      
      for(View<pat::Muon>::const_iterator iMuon2 = iMuon1+1; iMuon2 != thePATMuonHandle->end(); ++iMuon2) 
	{
	  if(iMuon1==iMuon2) continue;
	  
	  //opposite charge 
	  if( (iMuon1->charge())*(iMuon2->charge()) == 1) continue;

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
	  ParticleMass muon_mass = 0.10565837; //pdg mass
	  //ParticleMass psi_mass = 3.096916;
	  float muon_sigma = muon_mass*1.e-6;
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
	  
	  if(psi_vFit_vertex_noMC->chiSquared()>50.) continue;
	  if(psi_vFit_noMC->currentState().mass()<2.9 || psi_vFit_noMC->currentState().mass()>3.3) continue;
	  
	  //fill variables?iMuon1->track()->pt()
	  TLorentzVector JPsiKinematics;
      JPsiKinematics.SetPxPyPzE(psi_vFit_noMC->currentState().globalMomentum().x(),psi_vFit_noMC->currentState().globalMomentum().y(),psi_vFit_noMC->currentState().globalMomentum().z(),psi_vFit_noMC->currentState().kinematicParameters().energy()) ;
	  J_mass->push_back( psi_vFit_noMC->currentState().mass());
	  J_px->push_back( JPsiKinematics.Px() );
	  J_py->push_back( JPsiKinematics.Py());
	  J_pz->push_back(  JPsiKinematics.Pz());
	  J_energy->push_back(JPsiKinematics.E());
	  J_eta->push_back( JPsiKinematics.Eta() );
	  J_phi->push_back( JPsiKinematics.Phi());
	  J_pt->push_back(  JPsiKinematics.Pt());


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

for(View<pat::Muon>::const_iterator iMuon3 = thePATMuonHandle->begin(); iMuon3 != thePATMuonHandle->end(); ++iMuon3) 
    { 
      //if(iMuon3==jpsiMuon1->at(0)||iMuon3==jpsiMuon2->at(0))continue;
      
      for(View<pat::Muon>::const_iterator iMuon4 = iMuon3+1; iMuon4 != thePATMuonHandle->end(); ++iMuon4) 
	{
	  //if(iMuon4==jpsiMuon1->at(0)||iMuon4==jpsiMuon2->at(0))continue;
	  if(iMuon3==iMuon4) continue;
	  
	  //opposite charge 
	  if( (iMuon3->charge())*(iMuon4->charge()) == 1) continue;

	  TrackRef glbTrackP;	  
	  TrackRef glbTrackM;	  
	  
	  if(iMuon3->charge() == 1){ glbTrackP = iMuon3->track();}
	  if(iMuon3->charge() == -1){ glbTrackM = iMuon3->track();}
	  
	  if(iMuon4->charge() == 1) { glbTrackP = iMuon4->track();}
	  if(iMuon4->charge() == -1){ glbTrackM = iMuon4->track();}
	  
	  if( glbTrackP.isNull() || glbTrackM.isNull() ) 
	    {
	      //std::cout << "continue due to no track ref" << endl;
	      continue;
	    }

	  if(iMuon3->track()->pt()<2.5) continue;
	  if(iMuon4->track()->pt()<2.5) continue;

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
	  //if (iMuon3->isTrackerMuon() || iMuon4->isTrackerMuon())
	  //if (muon::isHighPtMuon(*iMuon3,bestVtx) || muon::isHighPtMuon(*iMuon4,bestVtx))
	  if (muon::isGoodMuon(*iMuon3,muon::TMLastStationAngTight) || muon::isGoodMuon(*iMuon4,muon::TMLastStationAngTight))
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
	  ParticleMass muon_mass = 0.10565837; //pdg mass
	  //ParticleMass psi_mass = 3.096916;
	  float muon_sigma = muon_mass*1.e-6;
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
	  
	  RefCountedKinematicParticle upsilon_vFit_noMC = psiVertexFitTree->currentParticle();
	  RefCountedKinematicVertex upsilon_vFit_vertex_noMC = psiVertexFitTree->currentDecayVertex();
	  
	  if( upsilon_vFit_vertex_noMC->chiSquared() < 0 )
	    {
	      //std::cout << "negative chisq from psi fit" << endl;
	      continue;
	    }
	  
	  //some loose cuts go here
	  
	  if(upsilon_vFit_vertex_noMC->chiSquared()>50.) continue;
	  if(upsilon_vFit_noMC->currentState().mass()<9.4603-0.3 || upsilon_vFit_noMC->currentState().mass()>9.4603+0.3) continue;
	  
	  //fill variables?iMuon3->track()->pt()
	  
	  TLorentzVector UpsilonKinematics;
      UpsilonKinematics.SetPxPyPzE(upsilon_vFit_noMC->currentState().globalMomentum().x(),upsilon_vFit_noMC->currentState().globalMomentum().y(),upsilon_vFit_noMC->currentState().globalMomentum().z(),upsilon_vFit_noMC->currentState().kinematicParameters().energy()) ;
	  U_mass->push_back( upsilon_vFit_noMC->currentState().mass() );
	  U_px->push_back( UpsilonKinematics.Px() );
	  U_py->push_back( UpsilonKinematics.Py());
	  U_pz->push_back(  UpsilonKinematics.Pz());
	  U_energy->push_back(UpsilonKinematics.E());
	  U_eta->push_back( UpsilonKinematics.Eta() );
	  U_phi->push_back( UpsilonKinematics.Phi());
	  U_pt->push_back(  UpsilonKinematics.Pt());
	  
	  U_px1->push_back(iMuon3->track()->px());
	  U_py1->push_back(iMuon3->track()->py());
	  U_pz1->push_back(iMuon3->track()->pz());
	  U_charge1->push_back(iMuon3->charge());
	  
	  U_px2->push_back(iMuon4->track()->px());
	  U_py2->push_back(iMuon4->track()->py());
	  U_pz2->push_back(iMuon4->track()->pz());
	  U_charge2->push_back(iMuon4->charge());
	  
	   // ************
	  
	  U_mu1soft->push_back(iMuon3->isSoftMuon(bestVtx) );
	  U_mu2soft->push_back(iMuon4->isSoftMuon(bestVtx) );
	  U_mu1tight->push_back(iMuon3->isTightMuon(bestVtx) );
	  U_mu2tight->push_back(iMuon4->isTightMuon(bestVtx) );
	  U_mu1PF->push_back(iMuon3->isPFMuon());
	  U_mu2PF->push_back(iMuon4->isPFMuon());
	  U_mu1loose->push_back(muon::isLooseMuon(*iMuon3));
	  U_mu2loose->push_back(muon::isLooseMuon(*iMuon4));
	  
	  U_mumC2->push_back( glbTrackP->normalizedChi2() );
	  //mumAngT->push_back( muon::isGoodMuon(*iMuon3,muon::TMLastStationAngTight) ); // 
	  U_mumNHits->push_back( glbTrackP->numberOfValidHits() );
	  U_mumNPHits->push_back( glbTrackP->hitPattern().numberOfValidPixelHits() );	       
	  U_mupC2->push_back( glbTrackM->normalizedChi2() );
	  //mupAngT->push_back( muon::isGoodMuon(*iMuon4,muon::TMLastStationAngTight) );  // 
	  U_mupNHits->push_back( glbTrackM->numberOfValidHits() );
	  U_mupNPHits->push_back( glbTrackM->hitPattern().numberOfValidPixelHits() );
	  U_mumdxy->push_back(glbTrackP->dxy(bestVtx.position()) );
	  U_mupdxy->push_back(glbTrackM->dxy(bestVtx.position()) );
	  U_mumdz->push_back(glbTrackP->dz(bestVtx.position()) );
	  U_mupdz->push_back(glbTrackM->dz(bestVtx.position()) );
	  U_muon_dca->push_back(dca);          	  
	  
	  nU++;	       
	  muonParticles.clear();
	  
	}
    }
  for(View<pat::PackedCandidate>::const_iterator iTrack1= thePATTrackHandle->begin(); iTrack1 != thePATTrackHandle->end();++iTrack1){
  	if(!(iTrack1->bestTrack())) continue;
  	if(iTrack1->pt()<0.8)continue;
  	if(iTrack1->eta()>2||iTrack1->eta()<-2)continue;

  	for(View<pat::PackedCandidate>::const_iterator iTrack2= iTrack1+1; iTrack2 != thePATTrackHandle->end();++iTrack2){
      if(iTrack1==iTrack2) continue;
      if((iTrack1->charge())*(iTrack2->charge())==1) continue;
      if(iTrack2->pt()<0.8)continue;

      Track glbTrackp;
      Track glbTrackn;
      	  
	  
	  if(!(iTrack2->bestTrack())) continue;
	  if(iTrack1->charge() == 1){ glbTrackp = *(iTrack1->bestTrack());}
	  if(iTrack1->charge() == -1){ glbTrackn = *(iTrack1->bestTrack());}
	  
	  if(iTrack2->charge() == 1) { glbTrackp = *(iTrack2->bestTrack());}
	  if(iTrack2->charge() == -1){ glbTrackn = *(iTrack2->bestTrack());}
	  
	  //if(!glbTrackp || !glbTrackn )continue;

	  if(!(glbTrackn.quality(reco::TrackBase::highPurity))) continue;
	  if(!(glbTrackp.quality(reco::TrackBase::highPurity))) continue;	

	  if(glbTrackp.normalizedChi2()>2) continue;
	  if(glbTrackn.normalizedChi2()>2) continue;
      if(glbTrackp.pt()<0.8) continue;
	  if(glbTrackn.pt()<0.8) continue;
      if(glbTrackp.eta()>2||glbTrackp.eta()<-2) continue;
      if(glbTrackn.eta()>2||glbTrackn.eta()<-2) continue;
	  if(glbTrackp.numberOfValidHits()<5) continue;
	  if(glbTrackn.numberOfValidHits()<5) continue;
    
      reco::TransientTrack kaon1TT((*theB).build(glbTrackp));
	  reco::TransientTrack kaon2TT((*theB).build(glbTrackn));

	  FreeTrajectoryState ka1State = kaon1TT.impactPointTSCP().theState();
	  FreeTrajectoryState ka2State = kaon2TT.impactPointTSCP().theState();

	  if( !kaon1TT.impactPointTSCP().isValid() || !kaon2TT.impactPointTSCP().isValid() ) continue;

	  // Measure distance between tracks at their closest approach
	  ClosestApproachInRPhi cApp;
	  cApp.calculate(ka1State, ka2State);
	  if( !cApp.status() ) continue;
	  float dca = fabs( cApp.distance() );	  
	  if (dca < 0. || dca > 0.5) continue;

	  ParticleMass kaon_mass = 0.493677; //pdg mass
	  
	  float kaon_sigma = 0.000016;

	  KinematicParticleFactoryFromTransientTrack pFactory;
	  
	  //initial chi2 and ndf before kinematic fits.
	  float chi = 0.;
	  float ndf = 0.;
	  vector<RefCountedKinematicParticle> kaonParticles;
	  try {
	    kaonParticles.push_back(pFactory.particle(kaon1TT,kaon_mass,chi,ndf,kaon_sigma));
	    kaonParticles.push_back(pFactory.particle(kaon2TT,kaon_mass,chi,ndf,kaon_sigma));
	  }
	  catch(...) { 
	    std::cout<<" Exception caught ... continuing 1 "<<std::endl; 
	    continue;
	  }
	  
	  KinematicParticleVertexFitter fitter;   
	  
	  RefCountedKinematicTree psiVertexFitTree;
	  try {
	    psiVertexFitTree = fitter.fit(kaonParticles); 
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
	  //creating the particle fitter
	  /*KinematicParticleFitter csFitter;

	  // creating the constraint
      float phi_m_pdg = 1.01946;
	  float phi_m_sigma = 0.000016;
	  KinematicConstraint * phi_c2 = new MassKinematicConstraint(phi_m_pdg,phi_m_sigma);

	  //the constrained fit:
	  psiVertexFitTree = csFitter.fit(phi_c2,psiVertexFitTree);*/
	  psiVertexFitTree->movePointerToTheTop();
	  
	  RefCountedKinematicParticle upsilon_vFit_noMC = psiVertexFitTree->currentParticle();
	  RefCountedKinematicVertex upsilon_vFit_vertex_noMC = psiVertexFitTree->currentDecayVertex();
	  
	  if( upsilon_vFit_vertex_noMC->chiSquared() < 0 )
	    {
	      //std::cout << "negative chisq from psi fit" << endl;
	      continue;
	    }
	  
	  //some loose cuts go here
	  
	  if(upsilon_vFit_vertex_noMC->chiSquared()>10.) continue;
	  if(upsilon_vFit_noMC->currentState().mass()<1.01946-0.01 || upsilon_vFit_noMC->currentState().mass()>1.01946+0.01) continue;
	  
	  //fill variables?iMuon3->track()->pt()
	  //cout<<"mass is"<<(upsilon_vFit_noMC->currentState().mass());
	  TLorentzVector PhiKinematics;
      PhiKinematics.SetPxPyPzE(upsilon_vFit_noMC->currentState().globalMomentum().x(),upsilon_vFit_noMC->currentState().globalMomentum().y(),upsilon_vFit_noMC->currentState().globalMomentum().z(),upsilon_vFit_noMC->currentState().kinematicParameters().energy()) ;
	  Phi_mass  ->push_back(upsilon_vFit_noMC->currentState().mass());
	  Phi_px    ->push_back(PhiKinematics.Px() );
	  Phi_py    ->push_back(PhiKinematics.Py());
	  Phi_pz    ->push_back(PhiKinematics.Pz());
	  Phi_energy->push_back(PhiKinematics.E());
	  Phi_eta   ->push_back(PhiKinematics.Eta() );
	  Phi_phi   ->push_back(PhiKinematics.Phi());
	  Phi_pt    ->push_back(PhiKinematics.Pt());
      nPhi++;
	  kaonParticles.clear();
      

    }
  } 

  //for (reco::PackedCandidate::const_iterator iTrack1= thePATTrackHandle->begin();  iTrack1 != thePATTrackHandle->end(); ++iTrack1){
  	
  
  if (nJ > 0 && nPhi>0 ) 
    {

      //std::cout << "filling tree" << endl;
      tree_->Fill();
    }

   nU = 0;
   nJ = 0; 
   nPhi = 0;

   J_mass->clear();  J_px->clear();  J_py->clear();  J_pz->clear();  
   J_energy->clear();  
   J_eta->clear();  
   J_phi->clear();  
   J_pt->clear();  
   J_px1->clear();  J_py1->clear();  J_pz1->clear(), J_charge1->clear();
   J_px2->clear();  J_py2->clear();  J_pz2->clear(), J_charge2->clear();

   U_mass->clear();  U_px->clear();  U_py->clear();  U_pz->clear(); 
   U_energy->clear();
   U_eta->clear();  
   U_phi->clear(); 
   U_pt->clear(); 
   U_px1->clear();  U_py1->clear();  U_pz1->clear(), U_charge1->clear();
   U_px2->clear();  U_py2->clear();  U_pz2->clear(), U_charge2->clear();
   mumC2->clear();
   mumNHits->clear(); mumNPHits->clear();
   mupC2->clear();
   mupNHits->clear(); mupNPHits->clear();
   mumdxy->clear(); mupdxy->clear(); mumdz->clear(); mupdz->clear(); muon_dca->clear();

   mu1soft->clear(); mu2soft->clear(); mu1tight->clear(); mu2tight->clear();
   mu1PF->clear(); mu2PF->clear(); mu1loose->clear(); mu2loose->clear(); 

   U_mumC2->clear();
   U_mumNHits->clear(); U_mumNPHits->clear();
   U_mupC2->clear();
   U_mupNHits->clear(); U_mupNPHits->clear();
   U_mumdxy->clear(); U_mupdxy->clear(); U_mumdz->clear(); U_mupdz->clear(); U_muon_dca->clear();

   U_mu1soft->clear(); U_mu2soft->clear(); U_mu1tight->clear(); U_mu2tight->clear();
   U_mu1PF->clear(); U_mu2PF->clear(); U_mu1loose->clear(); U_mu2loose->clear(); 
   Phi_mass->clear();
   Phi_px->clear();
   Phi_py->clear();
   Phi_pz->clear();
   Phi_energy->clear();
   Phi_eta->clear();
   Phi_phi->clear();
   Phi_pt->clear();
}


// ------------ method called once each job just before starting event loop  ------------

void 
miniAODmuons::beginJob()
{

  std::cout << "Beginning analyzer job with value of isMC= " << isMC_ << std::endl;

  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("ntuple"," J/psi ntuple");

  tree_->Branch("nJ",&nJ,"nJ/i");
  tree_->Branch("nU",&nU,"nU/i"); 
  tree_->Branch("J_mass", &J_mass);
  tree_->Branch("J_px", &J_px);
  tree_->Branch("J_py", &J_py);
  tree_->Branch("J_pz", &J_pz);
  tree_->Branch("J_eta", &J_eta);
  tree_->Branch("J_phi", &J_phi);
  tree_->Branch("J_pt", &J_pt);

  tree_->Branch("J_px1", &J_px1);
  tree_->Branch("J_py1", &J_py1);
  tree_->Branch("J_pz1", &J_pz1);
  tree_->Branch("J_charge1", &J_charge1);

  tree_->Branch("J_px2", &J_px2);
  tree_->Branch("J_py2", &J_py2);
  tree_->Branch("J_pz2", &J_pz2);
  tree_->Branch("J_charge2", &J_charge2);

  tree_->Branch("U_mass", &U_mass);
  tree_->Branch("U_px", &U_px);
  tree_->Branch("U_py", &U_py);
  tree_->Branch("U_pz", &U_pz);
  tree_->Branch("U_energy", &U_energy);
  tree_->Branch("U_eta", &U_eta);
  tree_->Branch("U_phi", &U_phi);
  tree_->Branch("U_pt", &U_pt);

  tree_->Branch("U_px1", &U_px1);
  tree_->Branch("U_py1", &U_py1);
  tree_->Branch("U_pz1", &U_pz1);
  tree_->Branch("U_charge1", &U_charge1);

  tree_->Branch("U_px2", &U_px2);
  tree_->Branch("U_py2", &U_py2);
  tree_->Branch("U_pz2", &U_pz2);
  tree_->Branch("U_charge2", &U_charge2); 

 
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

  tree_->Branch("U_mumC2",&U_mumC2);  
  tree_->Branch("U_mumNHits",&U_mumNHits);
  tree_->Branch("U_mumNPHits",&U_mumNPHits);
  tree_->Branch("U_mupC2",&U_mupC2);  
  tree_->Branch("U_mupNHits",&U_mupNHits);
  tree_->Branch("U_mupNPHits",&U_mupNPHits);
  tree_->Branch("U_mumdxy",&U_mumdxy);
  tree_->Branch("U_mupdxy",&U_mupdxy);
  tree_->Branch("U_mumdz",&U_mumdz);
  tree_->Branch("U_mupdz",&U_mupdz);
  tree_->Branch("U_muon_dca",&U_muon_dca);

  tree_->Branch("U_mu1soft",&U_mu1soft);
  tree_->Branch("U_mu2soft",&U_mu2soft);
  tree_->Branch("U_mu1tight",&U_mu1tight);
  tree_->Branch("U_mu2tight",&U_mu2tight);
  tree_->Branch("U_mu1PF",&U_mu1PF);
  tree_->Branch("U_mu2PF",&U_mu2PF);
  tree_->Branch("U_mu1loose",&U_mu1loose);
  tree_->Branch("U_mu2loose",&U_mu2loose);
  tree_->Branch("Phi_mass",&Phi_mass);
  tree_->Branch("Phi_px",&Phi_px);
  tree_->Branch("Phi_py",&Phi_py);
  tree_->Branch("Phi_pz",&Phi_pz);
  tree_->Branch("Phi_eta",&Phi_eta);
  tree_->Branch("Phi_phi",&Phi_phi);
  tree_->Branch("Phi_pt",&Phi_pt);
}




// ------------ method called once each job just after ending the event loop  ------------
void miniAODmuons::endJob() {
  tree_->GetDirectory()->cd();
  tree_->Write();
}

//define this as a plug-in
DEFINE_FWK_MODULE(miniAODmuons);
