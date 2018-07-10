#ifndef _jpsipipi_h
#define _jpsipipi_h

// system include files
#include <memory>

// user include files

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/Handle.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"

#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackFromFTSFactory.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h" // muy importante para MiniAOD

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/V0Candidate/interface/V0Candidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"

#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"

#include "RecoVertex/VertexPrimitives/interface/BasicSingleVertexState.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"

#include "TFile.h"
#include "TTree.h"


//
// class decleration
//

class jpsipipi : public edm::EDAnalyzer {
public:
  explicit jpsipipi(const edm::ParameterSet&);
  ~jpsipipi();
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  void printout(const RefCountedKinematicVertex& myVertex) const;
  void printout(const RefCountedKinematicParticle& myParticle) const;
  void printout(const RefCountedKinematicTree& myTree) const;

    // ----------member data ---------------------------
  
  edm::EDGetTokenT<edm::View<pat::Muon>> dimuon_Label;
  edm::EDGetTokenT<edm::View<pat::PackedCandidate>> trakCollection_label;
  edm::EDGetTokenT<reco::VertexCollection> primaryVertices_Label;

  bool isMC_;

  TTree*      tree_;
  
  std::vector<double>       *mumC2;
  std::vector<int>         *mumNHits, *mumNPHits; 
  std::vector<double>       *mupC2;
  std::vector<int>         *mupNHits, *mupNPHits;
  std::vector<double>       *mumdxy, *mupdxy, *mumdz, *mupdz;
  std::vector<double>       *muon_dca;

  std::vector<bool>        *mu1soft, *mu2soft, *mu1tight, *mu2tight;  
  std::vector<bool>        *mu1PF, *mu2PF, *mu1loose, *mu2loose;  
 
  unsigned int              nJ;
    
  std::vector<double>       *J_mass, *J_px, *J_py, *J_pz, *J_energy;

  std::vector<double>       *J_px1, *J_py1, *J_pz1;
  std::vector<double>       *J_px2, *J_py2, *J_pz2;
  std::vector<int>         *J_charge1, *J_charge2;
  std::vector<float>       *J_vertexFitChi2;
  std::vector<float>       *J_vertexFitNdf;

  unsigned int              nPiPair;
  
  std::vector<int>          *Pi_Hits1;
  std::vector<int>          *Pi_Hits2;
  std::vector<int>          *Pi_pixelHits1;
  std::vector<int>          *Pi_pixelHits2;
  std::vector<double>       *Pi_Eta1;
  std::vector<double>       *Pi_Eta2;
  std::vector<double>       *Pi_Phi1;
  std::vector<double>       *Pi_Phi2;
  std::vector<double>       *Pi_Pt1;
  std::vector<double>       *Pi_Pt2;
  std::vector<double>       *Pi_E1;
  std::vector<double>       *Pi_E2;
  std::vector<double>       *Pi_VertexChi2_1;
  std::vector<double>       *Pi_VertexChi2_2;
  std::vector<double>       *Pi_Lxy1;
  std::vector<double>       *Pi_Lxy2;
  std::vector<double>       *Pi_LxyErr1;
  std::vector<double>       *Pi_LxyErr2;
  std::vector<double>       *Jpipi_mass;
  std::vector<double>       *Jpi1_mass;
};
#endif
