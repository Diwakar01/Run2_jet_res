// -*- C++ -*-
//
// Package:    edanalyzer_for_resolution/edanalyzer
// Class:      edanalyzer
//
/**\class edanalyzer edanalyzer.cc edanalyzer_for_resolution/edanalyzer/plugins/edanalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Diwakar Vats
//         Created:  Fri, 20 May 2022 12:44:02 GMT
//
//

/*
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include <string>
#include <map>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <memory>

#include <Math/Vector3D.h>
#include "Math/LorentzVector.h"
#include "Math/Point3D.h"

#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "DataFormats/Luminosity/interface/LumiSummary.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "FWCore/Common/interface/TriggerNames.h"

//#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"
//#include "DataFormats/L1Trigger/interface/L1MuonParticleFwd.h"
//#include "DataFormats/L1Trigger/interface/L1EmParticle.h"
//#include "DataFormats/L1Trigger/interface/L1EmParticleFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
//#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapRecord.h"
#include "CondFormats/L1TObjects/interface/L1GtPrescaleFactors.h"
#include "CondFormats/DataRecord/interface/L1GtPrescaleFactorsAlgoTrigRcd.h"
#include "CondFormats/DataRecord/interface/L1GtPrescaleFactorsTechTrigRcd.h"

#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
//#include "DataFormats/TrackReco/interface/DeDxData.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"
//#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/METReco/interface/METFwd.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETFwd.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/JetCollection.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
//#include "HiggsCPinTauDecays/TauRefit/interface/RefitVertex.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/PatCandidates/interface/PATTauDiscriminator.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/GenMETFwd.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
//#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "RecoEgamma/EgammaTools/interface/ConversionInfo.h"
#include "RecoEgamma/EgammaTools/interface/ConversionFinder.h"

#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include <Math/Functions.h>
#include <Math/SVector.h>
#include <Math/SMatrix.h>

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/GeometrySurface/interface/SimpleCylinderBounds.h"
#include "DataFormats/GeometrySurface/interface/SimpleDiskBounds.h"
#include "DataFormats/GeometrySurface/interface/Cylinder.h"
#include "DataFormats/GeometrySurface/interface/Plane.h"
#include "DataFormats/GeometrySurface/interface/BoundCylinder.h"
#include "DataFormats/GeometrySurface/interface/BoundDisk.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/MaterialEffects/interface/PropagatorWithMaterial.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
//#include "DesyTauAnalyses/CandidateTools/interface/NSVfitStandaloneAlgorithm.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

//#include "EgammaAnalysis/ElectronTools/interface/EGammaMvaEleEstimatorCSA14.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "DataFormats/L1Trigger/interface/Muon.h"
#include "DataFormats/L1Trigger/interface/EGamma.h"
#include "DataFormats/L1Trigger/interface/Tau.h"
#include "DataFormats/L1Trigger/interface/Jet.h"

//includes for Helix parameter calculations
#include "DataFormats/TrackReco/interface/TrackBase.h"

#include "SimDataFormats/HTXS/interface/HiggsTemplateCrossSections.h"

//needed for tauspinner
//#include "TauSpinner/SimpleParticle.h"
//#include "TauSpinner/tau_reweight_lib.h"







//#include "DesyTauAnalyses/NTupleMaker/plugins/NTupleMaker.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
//#include "AnalysisDataFormats/TauAnalysis/interface/CompositePtrCandidateT1T2MEtFwd.h"
//#include "TauAnalysis/CandidateTools/interface/NSVfitAlgorithmBase.h"
//#include "AnalysisDataFormats/TauAnalysis/interface/NSVfitEventHypothesisFwd.h"
//#include "AnalysisDataFormats/TauAnalysis/interface/NSVfitEventHypothesisBaseFwd.h"
//#include "AnalysisDataFormats/TauAnalysis/interface/NSVfitEventHypothesisByIntegration.h"
//#include "AnalysisDataFormats/TauAnalysis/interface/NSVfitResonanceHypothesisBase.h"
//#include "AnalysisDataFormats/TauAnalysis/interface/NSVfitSingleParticleHypothesis.h"
#include <DataFormats/RecoCandidate/interface/IsoDepositVetos.h>
#include <DataFormats/METReco/interface/GenMET.h>
#include <DataFormats/HLTReco/interface/TriggerTypeDefs.h>
#include "RecoBTag/BTagTools/interface/SignedImpactParameter3D.h"
#include <DataFormats/TrackReco/interface/Track.h>
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
//#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"

#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "PhysicsTools/JetMCUtils/interface/JetMCTag.h"
//#include "TauAnalysis/CandidateTools/interface/CompositePtrCandidateT1T2MEtProducer.h"
//#include "AnalysisDataFormats/TauAnalysis/interface/CompositePtrCandidateT1T2MEt.h"
//#include "AnalysisDataFormats/TauAnalysis/interface/CompositePtrCandidateT1T2MEtFwd.h"
#include "SimDataFormats/GeneratorProducts/interface/GenFilterInfo.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
//#include "DesyTauAnalyses/CandidateTools/interface/candidateAuxFunctions.h"
//#include "DesyTauAnalyses/NTupleMaker/interface/idAlgos.h"

//#include "ICTauSpinnerProducer.hh"


#include <TString.h>
#include <Compression.h>

*/


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"
#include "TTree.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/JetCollection.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

// using namespace reco;
// using namespace pat;
// using namespace std;
using reco::TrackCollection;

class edanalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit edanalyzer(const edm::ParameterSet&);
      ~edanalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      virtual float deltaR(float eta1, float phi1, float eta2, float phi2);
      // ----------member data ---------------------------
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> TriggerObjectCollectionToken_;  
  edm::EDGetTokenT<pat::JetCollection> JetCollectionToken_;
  double Pi = 3.14159265359;

  edm::Service<TFileService> fs;
  TH1F* h1;
  TH1F* h2;
  TH1F* h3;
  TH1F* h4;
  TH1F* h5;
  TH1F* h6;
  TH1F* h7;
  TH1F* h8;

  int count;
  /*
  std::vector<float> trigobject_phi;
  std::vector<float> trigobject_pt;
  std::vector<float> trigobject_eta;

  std::vector<float> pfjet_pt;
  std::vector<float> pfjet_eta;
  std::vector<float> pfjet_phi;

  edm::Service<TFileService> fs;
  TTree* t;
  */
  //  edm::EDGetTokenT<pat::JetCollection> JetCollectionToken_;
  //edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> TriggerObjectCollectionToken_;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
edanalyzer::edanalyzer(const edm::ParameterSet& iConfig)
  :
  
  TriggerObjectCollectionToken_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("TriggerObjectCollectionTag"))),
  
   JetCollectionToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("JetCollectionTag")))

{
   //now do what ever initialization is needed
  
  h1 = fs->make<TH1F>("h1","h1",50,-1000,1000);
  h2 = fs->make<TH1F>("h2","h2",50,-100,100);
  h3 = fs->make<TH1F>("h3","h3",50,-10,10);
  h4 = fs->make<TH1F>("h4","h4",50,-5,5);
  h5 = fs->make<TH1F>("h5","h5",50,-1000,1000);
  h6 = fs->make<TH1F>("h6","h6",50,-100,100);
  h7 = fs->make<TH1F>("h7","h7",50,-10,10);
  h8 = fs->make<TH1F>("h8","h8",50,-5,5);

  count = 0;
  /*
  t->Branch("pfjet_pt","std::vector<float>",&pfjet_pt);
  t->Branch("pfjet_eta","std::vector<float>",&pfjet_eta);
  t->Branch("pfjet_phi","std::vector<float>",&pfjet_phi);
  t->Branch("trigobject_pt","std::vector<float>",&trigobject_pt);
  t->Branch("trigobject_eta","std::vector<float>",&trigobject_eta);
  t->Branch("trigobject_phi","std::vector<float>",&trigobject_phi);
  */
}


edanalyzer::~edanalyzer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
edanalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace std;
   using namespace edm;
   /*
   edm::Handle<pat::JetCollection> pfjets;
   iEvent.getByToken(JetCollectionToken_, pfjets);
   edm::Handle<pat::TriggerObjectStandAloneCollection_> triggerObjects;
   iEvent.getByToken(TriggerObjectCollectionToken_, triggerObjects);
   */
   
   edm::Handle<pat::JetCollection> pfjets;
   iEvent.getByToken(JetCollectionToken_, pfjets);
   
   edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
   iEvent.getByToken(TriggerObjectCollectionToken_, triggerObjects);
   
    /*
   if(pfjets.isValid())
     {
       for(unsigned i = 0 ; i < pfjets->size() ; i++)
	 {
	   std::cout<<std::endl<<"i="<<i<<std::endl;
	   std::cout<<std::endl<<"(*pfjets)[i].pt()="<<(*pfjets)[i].pt()<<std::endl;
	   pfjet_pt.push_back((*pfjets)[i].pt());
	   pfjet_eta.push_back((*pfjets)[i].eta());
	   pfjet_phi.push_back((*pfjets)[i].phi());
	 }
     }


   for (unsigned int iTO=0; iTO<triggerObjects->size(); ++iTO) {
     trigobject_pt.push_back((*triggerObjects)[iTO].pt());
     trigobject_eta.push_back((*triggerObjects)[iTO].eta());
     trigobject_phi.push_back((*triggerObjects)[iTO].phi());
   }

   */

   for(unsigned int ijet = 0; ijet < pfjets->size(); ijet++){
     //     cout<<endl<<"jet no = "<<ijet<<endl;//<<"pfjet_pt->at(ijet)="<<pfjet_pt->at(ijet)<<endl;
     for(unsigned int itrigobject = 0; itrigobject < triggerObjects->size(); itrigobject++){
       //       cout<<endl<<"trig obj no = "<<itrigobject<<endl;//<<"trigobject_pt->at(itrigobject) = "<<trigobject_pt->at(itrigobject)<<endl;
       //       cout<<endl<<"deltaR="<<deltaR((*pfjets)[ijet].eta(), (*pfjets)[ijet].phi(), (*triggerObjects)[itrigobject].eta(), (*triggerObjects)[itrigobject].phi())<<endl;
       if(deltaR((*pfjets)[ijet].eta(), (*pfjets)[ijet].phi(), (*triggerObjects)[itrigobject].eta(), (*triggerObjects)[itrigobject].phi()) < 0.4){
	 //	 cout<<endl<<"filling h with jet, trigobj "<<ijet<<" ,"<<itrigobject<<endl;// with "<<pfjet_pt->at(ijet)-trigobject_pt->at(itrigobject)<<endl;
	 count++;
	   // h1->Fill((*pfjets)[ijet].pt()-(*triggerObjects)[itrigobject].pt());
	   // h2->Fill((*pfjets)[ijet].pt()-(*triggerObjects)[itrigobject].pt());
	   // h3->Fill((*pfjets)[ijet].pt()-(*triggerObjects)[itrigobject].pt());
	   // h4->Fill((*pfjets)[ijet].pt()-(*triggerObjects)[itrigobject].pt());
	   h1->Fill((*pfjets)[ijet].pt()/(*triggerObjects)[itrigobject].pt());
	   h2->Fill((*pfjets)[ijet].pt()/(*triggerObjects)[itrigobject].pt());
	   h3->Fill((*pfjets)[ijet].pt()/(*triggerObjects)[itrigobject].pt());
	   h4->Fill((*pfjets)[ijet].pt()/(*triggerObjects)[itrigobject].pt());
	   h5->Fill((*pfjets)[ijet].pt()-(*triggerObjects)[itrigobject].pt());
	   h6->Fill((*pfjets)[ijet].pt()-(*triggerObjects)[itrigobject].pt());
	   h7->Fill((*pfjets)[ijet].pt()-(*triggerObjects)[itrigobject].pt());
	   h8->Fill((*pfjets)[ijet].pt()-(*triggerObjects)[itrigobject].pt());

	  break;
	}
      }
    }
   //cout<<endl<<"count = "<<count<<endl;
}


   //   return trigobject_count;

//   t->Fill();

   /*
    Handle<TrackCollection> tracks;
    iEvent.getByToken(tracksToken_, tracks);
    for(TrackCollection::const_iterator itTrack = tracks->begin();
        itTrack != tracks->end();
        ++itTrack) {
      // do something with track parameters, e.g, plot the charge.
      // int charge = itTrack->charge();
    }

#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
   */


// ------------ method called once each job just before starting event loop  ------------
void
edanalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
edanalyzer::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
edanalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

float edanalyzer::deltaR(float eta1, float phi1, float eta2, float phi2){
  float deta = fabs(eta1 - eta2);
  float dphi = fabs(phi1 - phi2);
  if(dphi > Pi)dphi = 2*Pi - dphi;
  float dr = sqrt(deta*deta + dphi*dphi);
  return dr;
}

//define this as a plug-in
DEFINE_FWK_MODULE(edanalyzer);
