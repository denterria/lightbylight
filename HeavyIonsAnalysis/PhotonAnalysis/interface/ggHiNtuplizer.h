#ifndef ggHiNtuplizer_h
#define ggHiNtuplizer_h

#include "TTree.h"

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/HIPhotonIsolation.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "RecoEgamma/EgammaTools/interface/ConversionFinder.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "CondFormats/DataRecord/interface/GBRWrapperRcd.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalDetId/interface/ESDetId.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HepMCCandidate/interface/PdfInfo.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/METReco/interface/BeamHaloSummary.h"
#include "DataFormats/METReco/interface/HcalNoiseSummary.h"
#include "DataFormats/MuonReco/interface/MuonCocktails.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/Scalers/interface/DcsStatus.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/TruncatedPyramid.h"
#include "Geometry/EcalAlgo/interface/EcalPreshowerGeometry.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "PhysicsTools/CandUtils/interface/AddFourMomenta.h"
#include "RecoCaloTools/Navigation/interface/EcalPreshowerNavigator.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "RecoEcal/EgammaCoreTools/interface/Mustache.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"
#include "RecoEgamma/EgammaIsolationAlgos/plugins/EgammaTowerIsolationProducer.h"
#include "RecoEgamma/EgammaTools/interface/ConversionFinder.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "DataFormats/CaloTowers/interface/CaloTower.h"
#include "DataFormats/CaloTowers/interface/CaloTowerFwd.h"
#include <iostream>
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterFunctionFactory.h" 

using namespace std;


class ggHiNtuplizer : public edm::EDAnalyzer {

 public:

   ggHiNtuplizer(const edm::ParameterSet&);
   virtual ~ggHiNtuplizer() {};

 private:

   virtual void analyze(const edm::Event&, const edm::EventSetup&);

   void fillGenParticles (const edm::Event&);
   void fillGenPileupInfo(const edm::Event&);
   void fillElectrons    (const edm::Event&, const edm::EventSetup&, math::XYZPoint& pv);
   void fillPhotons      (const edm::Event&, const edm::EventSetup&, math::XYZPoint& pv);
   void fillGsfPhotons   (const edm::Event&, const edm::EventSetup&, math::XYZPoint& pv);
   void fillMuons        (const edm::Event&, const edm::EventSetup&, math::XYZPoint& pv);
   void fillCaloTower    (const edm::Event&, const edm::EventSetup&, math::XYZPoint& pv);
   void fillElectronTracks   (const edm::Event&, const edm::EventSetup&, math::XYZPoint& pv);
   void fillGeneralTracks    (const edm::Event&, const edm::EventSetup&, math::XYZPoint& pv);
   void fillSuperClusters    (const edm::Event&, const edm::EventSetup&, math::XYZPoint& pv);

   // Et and pT sums
   float getGenCalIso(edm::Handle<vector<reco::GenParticle> >&, reco::GenParticleCollection::const_iterator, float dRMax, bool removeMu, bool removeNu);
   float getGenTrkIso(edm::Handle<vector<reco::GenParticle> >&, reco::GenParticleCollection::const_iterator, float dRMax);

   // switches
   bool doGenParticles_;
   bool runOnParticleGun_;
   bool useValMapIso_;
   bool doPfIso_;
   bool doVsIso_; // also requires above boolean to make sense
   bool doVID_;

   // handles to collections of objects
   edm::EDGetTokenT<vector<PileupSummaryInfo> >    genPileupCollection_;
   edm::EDGetTokenT<vector<reco::GenParticle> >    genParticlesCollection_;
   edm::EDGetTokenT<edm::View<reco::GsfElectron> > gsfElectronsCollection_;
   edm::EDGetTokenT<reco::GsfElectronCollection>   input_electrons_;
   edm::EDGetTokenT<reco::GsfTrackCollection>  gsfTracks_;
   edm::EDGetTokenT<reco::TrackCollection>    genTracks_;
   edm::EDGetTokenT<edm::ValueMap<bool> > eleVetoIdMapToken_;
   edm::EDGetTokenT<edm::ValueMap<bool> > eleLooseIdMapToken_;
   edm::EDGetTokenT<edm::ValueMap<bool> > eleMediumIdMapToken_;
   edm::EDGetTokenT<edm::ValueMap<bool> > eleTightIdMapToken_;
   edm::EDGetTokenT<edm::View<reco::Photon> >      recoPhotonsCollection_;
   edm::EDGetTokenT<edm::View<reco::Photon> >      recoHyPhotonsCollection_;
   edm::EDGetTokenT<edm::ValueMap<reco::HIPhotonIsolation> > recoPhotonsHiIso_;
   edm::EDGetTokenT<edm::ValueMap<reco::HIPhotonIsolation> > recoHyPhotonsHiIso_;
   edm::EDGetTokenT<edm::View<reco::Muon> >        recoMuonsCollection_;
   //edm::EDGetTokenT<edm::SortedCollection<CaloTower,edm::StrictWeakOrdering<CaloTower>>> CaloTowerCollection_;
   edm::EDGetTokenT<edm::SortedCollection<CaloTower>>        CaloTowerCollection_;
   edm::EDGetTokenT<reco::SuperClusterCollection>      hybridsc_;
   edm::EDGetTokenT<reco::SuperClusterCollection>      mult55sc_;
   edm::EDGetTokenT<vector<reco::Vertex> >         vtxCollection_;
   edm::EDGetTokenT<double> rhoToken_;
   edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
   edm::EDGetTokenT<reco::ConversionCollection> conversionsToken_;
   edm::EDGetTokenT<edm::View<reco::PFCandidate> >    pfCollection_;
   edm::EDGetTokenT<edm::ValueMap<reco::VoronoiBackground> > voronoiBkgCalo_;
   edm::EDGetTokenT<edm::ValueMap<reco::VoronoiBackground> > voronoiBkgPF_;

   EffectiveAreas effectiveAreas_;

   TTree*         tree_;

   // variables associated with tree branches
   UInt_t          run_;
   ULong64_t       event_;
   UInt_t          lumis_;
   Bool_t         isData_;

   // Vtx info
   Int_t          nVtx_;
   vector<int>    isFakeVtx_;
   vector<float>  xVtx_;
   vector<float>  yVtx_;
   vector<float>  zVtx_;
   
   // PileupSummaryInfo
   Int_t          nPUInfo_;
   vector<int>    nPU_;
   vector<int>    puBX_;
   vector<float>  puTrue_;

   // reco::GenParticle
   Int_t          nMC_;
   vector<int>    mcPID_;
   vector<int>    mcStatus_;
   vector<float>  mcVtx_x_;
   vector<float>  mcVtx_y_;
   vector<float>  mcVtx_z_;
   vector<float>  mcPt_;
   vector<float>  mcP_;
   vector<float>  mcEta_;
   vector<float>  mcPhi_;
   vector<float>  mcE_;
   vector<float>  mcEt_;
   vector<float>  mcMass_;
   vector<int>    mcParentage_;
   vector<int>    mcMomPID_;
   vector<float>  mcMomPt_;
   vector<float>  mcMomEta_;
   vector<float>  mcMomPhi_;
   vector<float>  mcMomMass_;
   vector<int>    mcGMomPID_;
   vector<int>    mcIndex_;
   vector<float>  mcCalIsoDR03_;
   vector<float>  mcCalIsoDR04_;
   vector<float>  mcTrkIsoDR03_;
   vector<float>  mcTrkIsoDR04_;

   // reco::GedGsfElectron
   Int_t          nEle_;
   vector<int>    eleCharge_;
   vector<int>    eleChargeConsistent_;
   vector<int>    eleSCPixCharge_;
   vector<int>    eleCtfCharge_;
   vector<float>  eleEn_;
   vector<float>  eleD0_;
   vector<float>  eleDz_;
   vector<float>  eleD0Err_;
   vector<float>  eleDzErr_;
   vector<float>  eleTrkPt_;
   vector<float>  eleTrkEta_;
   vector<float>  eleTrkPhi_;
   vector<int>    eleTrkCharge_;
   vector<float>  eleTrkChi2_;
   vector<float>  eleTrkNdof_;
   vector<float>  eleTrkNormalizedChi2_;
   vector<int>    eleTrkValidHits_;
   vector<int>    eleTrkLayers_;
   vector<float>  eleP_;
   vector<float>  eleP_atVtx_;
   vector<float>  elePt_;
   vector<float>  eleEta_;
   vector<float>  elePhi_;
   vector<float>  eleSCx_;
   vector<float>  eleSCy_;
   vector<float>  eleSCz_;
   vector<float>  eleSCEn_;
   vector<float>  eleESEn_;
   vector<float>  eleSCEta_;
   vector<float>  eleSCPhi_;
   vector<float>  eleSCRawEn_;
   vector<float>  eleSCEtaWidth_;
   vector<float>  eleSCPhiWidth_;
   vector<float>  eleHoverE_;
   vector<float>  eleHoverEBc_;
   vector<float>  eleEoverP_;
   vector<float>  eleEoverPInv_;
   vector<float>  eleBrem_;
   vector<float>  eledEtaAtVtx_;
   vector<float>  eledPhiAtVtx_;
   vector<float>  eledEtaSeedAtVtx_;
   vector<float>  eleSigmaIEtaIEta_;
   vector<float>  eleSigmaIEtaIEta_2012_;
   vector<float>  eleSigmaIPhiIPhi_;
// vector<int>    eleConvVeto_;     // TODO: not available in reco::
   vector<int>    eleMissHits_;
   vector<float>  eleTrackIso_;
   vector<float>  eleECalIso_;
   vector<float>  eleHCalIso_;
   vector<float>  eleECalDriven_;
   vector<float>  eleESEffSigmaRR_;
   vector<float>  elePFChIso_;
   vector<float>  elePFPhoIso_;
   vector<float>  elePFNeuIso_;
   vector<float>  elePFPUIso_;
   vector<float>  elePFChIso03_;
   vector<float>  elePFPhoIso03_;
   vector<float>  elePFNeuIso03_;
   vector<float>  elePFChIso04_;
   vector<float>  elePFPhoIso04_;
   vector<float>  elePFNeuIso04_;
   vector<float>  eleR9_;
   vector<float>  eleE3x3_;
   vector<float>  eleE5x5_;
   vector<float>  eleR9Full5x5_;
   vector<float>  eleE3x3Full5x5_;
   vector<float>  eleE5x5Full5x5_;
   vector<int>	  NClusters_;
   vector<int>    NEcalClusters_;
   vector<float>  eleSeedEn_;
   vector<float>  eleSeedEta_;
   vector<float>  eleSeedPhi_;
   vector<float>  eleSeedCryEta_;
   vector<float>  eleSeedCryPhi_;
   vector<float>  eleSeedCryIeta_;
   vector<float>  eleSeedCryIphi_;
   vector<float>  eleBC1E_;
   vector<float>  eleBC1Eta_;
   vector<float>  eleBC2E_;
   vector<float>  eleBC2Eta_;
   vector<int>    eleIDVeto_; //50nsV1 is depreacated; updated in 76X
   vector<int>    eleIDLoose_;
   vector<int>    eleIDMedium_;
   vector<int>    eleIDTight_;
   vector<int>    elepassConversionVeto_;
   vector<float>    eleEffAreaTimesRho_;

   // conversion variables
   Int_t           nConv_;
   vector<int>     conv_eleind_;
   vector<int>     conv_phoind_;
   vector<float>   conv_vtxprob_;
   vector<float>   conv_lxy_;
   vector<int>     conv_nhits_bvtx_;

   

   // reco::Photon
   Int_t          nPho_;
   vector<float>  phoE_;
   vector<float>  phoEt_;
   vector<float>  phoEta_;
   vector<float>  phoPhi_;
   vector<int>    phoSCSize_;
   vector<float>  phoSCE_;
   vector<float>  phoSCEt_;
   vector<float>  phoSCRawE_;
   vector<float>  phoSCRawEt_;
   vector<float>  phoESEn_;
   vector<float>  phoSCEta_;
   vector<float>  phoSCPhi_;
   vector<float>  phoSCx_;
   vector<float>  phoSCy_;
   vector<float>  phoSCz_;
   vector<float>  phoSCEtaWidth_;
   vector<float>  phoSCPhiWidth_;
   vector<float>  phoSCBrem_;
   vector<int>    phohasPixelSeed_;
   vector<int>    phopassConversionVeto_;
   vector<int>    phoEleVeto_;         // TODO: not available in reco::
   vector<float>  phoR9_;
   vector<float>  phoHadTowerOverEm_;
   vector<float>  phoHoverE_;
   vector<float>  phoSigmaIEtaIEta_;
// vector<float>  phoSigmaIEtaIPhi_;   // TODO: not available in reco::
// vector<float>  phoSigmaIPhiIPhi_;   // TODO: not available in reco::
   vector<float>  phoE1x3_;
   vector<float>  phoE2x2_;
   vector<float>  phoE3x3_;
   vector<float>  phoE2x5Max_;
   vector<float>  phoE1x5_;
   vector<float>  phoE2x5_;
   vector<float>  phoE5x5_;
   vector<float>  phoMaxEnergyXtal_;
   vector<float>  phoSigmaEtaEta_;
   vector<float>  phoR1x5_;
   vector<float>  phoR2x5_;
   vector<float>  phoESEffSigmaRR_;
   vector<float>  phoSigmaIEtaIEta_2012_;
   vector<float>  phoSigmaIEtaIPhi_2012_;
   vector<float>  phoSigmaIPhiIPhi_2012_;
   vector<float>  phoE1x3_2012_;
   vector<float>  phoE2x2_2012_;
   vector<float>  phoE3x3_2012_;
   vector<float>  phoE2x5Max_2012_;
   vector<float>  phoE5x5_2012_;
   vector<float>  phoBC1E_;
   vector<float>  phoBC1Eta_;
   vector<float>  phoBC2E_;
   vector<float>  phoBC2Eta_;
   vector<float>  pho_ecalClusterIsoR2_;
   vector<float>  pho_ecalClusterIsoR3_;
   vector<float>  pho_ecalClusterIsoR4_;
   vector<float>  pho_ecalClusterIsoR5_;
   vector<float>  pho_hcalRechitIsoR1_;
   vector<float>  pho_hcalRechitIsoR2_;
   vector<float>  pho_hcalRechitIsoR3_;
   vector<float>  pho_hcalRechitIsoR4_;
   vector<float>  pho_hcalRechitIsoR5_;
   vector<float>  pho_trackIsoR1PtCut20_;
   vector<float>  pho_trackIsoR2PtCut20_;
   vector<float>  pho_trackIsoR3PtCut20_;
   vector<float>  pho_trackIsoR4PtCut20_;
   vector<float>  pho_trackIsoR5PtCut20_;
   vector<float>  pho_swissCrx_;
   vector<float>  pho_seedTime_;
   vector<int>    pho_genMatchedIndex_;

  vector<int>    phoIsConv_;
  vector<int>    phoNConv_;
  vector<int>    phoConvNTrks_;
  vector<float>  phoConvInvMass_;
  vector<float>  phoConvCotTheta_;
  vector<float>  phoConvEoverP_;
  vector<float>  phoConvZofPVfromTrks_;
  vector<float>  phoConvMinDist_;
  vector<float>  phoConvdPhiAtVtx_;
  vector<float>  phoConvdPhiAtCalo_;
  vector<float>  phoConvdEtaAtCalo_;
  vector<float>  phoConvTrkd0_x_;
  vector<float>  phoConvTrkd0_y_;
  vector<float>  phoConvTrkPin_x_;
  vector<float>  phoConvTrkPin_y_;
  vector<float>  phoConvTrkPout_x_;
  vector<float>  phoConvTrkPout_y_;
  vector<float>  phoConvTrkdz_x_;
  vector<float>  phoConvTrkdz_y_;
  vector<float>  phoConvTrkdzErr_x_;
  vector<float>  phoConvTrkdzErr_y_;
  vector<float>  phoConvChi2_;
  vector<float>  phoConvChi2Prob_;
  vector<float>  phoConvCharge1_;
  vector<float>  phoConvCharge2_;
  vector<int>    phoConvValidVtx_;
  vector<float>  phoConvLikeLihood_;
  vector<float>  phoConvP4_0_;
  vector<float>  phoConvP4_1_;
  vector<float>  phoConvP4_2_;
  vector<float>  phoConvP4_3_;
  vector<float>  phoConvVtx_x_;
  vector<float>  phoConvVtx_y_;
  vector<float>  phoConvVtx_z_;
  vector<float>  phoConvVtxErr_x_;
  vector<float>  phoConvVtxErr_y_;
  vector<float>  phoConvVtxErr_z_;
  vector<float>  phoConvPairMomentum_x_;
  vector<float>  phoConvPairMomentum_y_;
  vector<float>  phoConvPairMomentum_z_;
  vector<float>  phoConvRefittedMomentum_x_;
  vector<float>  phoConvRefittedMomentum_y_;
  vector<float>  phoConvRefittedMomentum_z_;
  vector<int>    SingleLegConv_;


   // reco::Muon
   Int_t          nMu_;
   vector<float>  muPt_;
   vector<float>  muEta_;
   vector<float>  muPhi_;
   vector<int>    muCharge_;
   vector<int>    muType_;
   vector<int>    muIsGood_;
   vector<float>  muD0_;
   vector<float>  muDz_;
   vector<float>  muChi2NDF_;
   vector<float>  muInnerD0_;
   vector<float>  muInnerDz_;
   vector<int>    muTrkLayers_;
   vector<int>    muPixelLayers_;
   vector<int>    muPixelHits_;
   vector<int>    muMuonHits_;
   vector<int>    muTrkQuality_;
   vector<int>    muStations_;
   vector<float>  muIsoTrk_;
   vector<float>  muPFChIso_;
   vector<float>  muPFPhoIso_;
   vector<float>  muPFNeuIso_;
   vector<float>  muPFPUIso_;

   // reco::calotower
   Int_t          nTower_;
   Int_t          nTower_barrel_;
   vector<float> CaloTower_hadE_;
   vector<float> CaloTower_emE_;
   vector<float> CaloTower_e_;
   vector<float> CaloTower_et_;
   vector<float> CaloTower_eta_;
   vector<float> CaloTower_phi_;

  // reco::Photon
   Int_t          nHyPho_;
   vector<float>  hyphoE_;
   vector<float>  hyphoEt_;
   vector<float>  hyphoEta_;
   vector<float>  hyphoPhi_;
   vector<int>    hyphoSCSize_;
   vector<float>  hyphoSCE_;
   vector<float>  hyphoSCEt_;
   vector<float>  hyphoSCRawE_;
   vector<float>  hyphoSCRawEt_;
   vector<float>  hyphoESEn_;
   vector<float>  hyphoSCEta_;
   vector<float>  hyphoSCPhi_;
   vector<float>  hyphoSCx_;
   vector<float>  hyphoSCy_;
   vector<float>  hyphoSCz_;
   vector<float>  hyphoSCEtaWidth_;
   vector<float>  hyphoSCPhiWidth_;
   vector<float>  hyphoSCBrem_;
   vector<int>    hyphohasPixelSeed_;
   vector<int>    hyphopassConversionVeto_;
   vector<int>    hyphoEleVeto_;         // TODO: not available in reco::
   vector<float>  hyphoR9_;
   vector<float>  hyphoHadTowerOverEm_;
   vector<float>  hyphoHoverE_;
   vector<float>  hyphoSigmaIEtaIEta_;
// vector<float>  hyphoSigmaIEtaIPhi_;   // TODO: not available in reco::
// vector<float>  hyphoSigmaIPhiIPhi_;   // TODO: not available in reco::
   vector<float>  hyphoE1x3_;
   vector<float>  hyphoE2x2_;
   vector<float>  hyphoE3x3_;
   vector<float>  hyphoE2x5Max_;
   vector<float>  hyphoE1x5_;
   vector<float>  hyphoE2x5_;
   vector<float>  hyphoE5x5_;
   vector<float>  hyphoMaxEnergyXtal_;
   vector<float>  hyphoSigmaEtaEta_;
   vector<float>  hyphoR1x5_;
   vector<float>  hyphoR2x5_;
   vector<float>  hyphoESEffSigmaRR_;
   vector<float>  hyphoSigmaIEtaIEta_2012_;
   vector<float>  hyphoSigmaIEtaIPhi_2012_;
   vector<float>  hyphoSigmaIPhiIPhi_2012_;
   vector<float>  hyphoE1x3_2012_;
   vector<float>  hyphoE2x2_2012_;
   vector<float>  hyphoE3x3_2012_;
   vector<float>  hyphoE2x5Max_2012_;
   vector<float>  hyphoE5x5_2012_;
   vector<float>  hyphoBC1E_;
   vector<float>  hyphoBC1Eta_;
   vector<float>  hyphoBC2E_;
   vector<float>  hyphoBC2Eta_;
   vector<float>  hypho_ecalClusterIsoR2_;
   vector<float>  hypho_ecalClusterIsoR3_;
   vector<float>  hypho_ecalClusterIsoR4_;
   vector<float>  hypho_ecalClusterIsoR5_;
   vector<float>  hypho_hcalRechitIsoR1_;
   vector<float>  hypho_hcalRechitIsoR2_;
   vector<float>  hypho_hcalRechitIsoR3_;
   vector<float>  hypho_hcalRechitIsoR4_;
   vector<float>  hypho_hcalRechitIsoR5_;
   vector<float>  hypho_trackIsoR1PtCut20_;
   vector<float>  hypho_trackIsoR2PtCut20_;
   vector<float>  hypho_trackIsoR3PtCut20_;
   vector<float>  hypho_trackIsoR4PtCut20_;
   vector<float>  hypho_trackIsoR5PtCut20_;
   vector<float>  hypho_swissCrx_;
   vector<float>  hypho_seedTime_;
   vector<int>    hypho_genMatchedIndex_;


   // reco::GsfTrack elctron tracks
    Int_t         ngsfEle_;
    vector<float> elegsfTrkPt_;            
    vector<float> elegsfTrkP_;            
    vector<float> elegsfTrkEta_;           
    vector<float> elegsfTrkPhi_;           
    vector<int>   elegsfTrkCharge_;        
    vector<float> elegsfTrkChi2_;          
    vector<float> elegsfTrkNdof_;          
    vector<float> elegsfTrkNormalizedChi2_;
    vector<int>   elegsfTrkValidHits_;   
    vector<int>   elegsfTrkMissHits_;       
    vector<int>   elegsfTrkLayers_;    
    vector<float> elegsfD0_;
    vector<float> elegsfDz_;              
    vector<float> elegsfD0Err_;            
    vector<float> elegsfDzErr_;       

   // reco::general track     
    Int_t         ngenTrk_;
    vector<float> gentrkPt_;            
    vector<float> gentrkP_;            
    vector<float> gentrkEta_;           
    vector<float> gentrkPhi_;     
    vector<int>   gentrkcharge_; 
    vector<float> gentrkvx_;            
    vector<float> gentrkvy_;           
    vector<float> gentrkvz_;     
    vector<float> gentrknormchi2_;     
    vector<float> gentrkchi2_;  
    vector<float> gentrkd0_;                         
    vector<float> gentrkdxy_;                      
    vector<float> gentrkdz_;  
    vector<float> gentrkdxyError_;  
    vector<float> gentrkdzError_;   
    vector<int>   gentrkValidHits_;           
    vector<int>   gentrkMissHits_;   
    vector<int>   gentrkPurity_;     
    
    //reco::hybrid supercluster
    Int_t          nsc_hybrid_;
    vector<float>  sc_hybrid_E_;
    vector<float>  sc_hybrid_Et_;
    vector<float>  sc_hybrid_Eta_;
    vector<float>  sc_hybrid_Phi_;
    vector<float>  sc_hybrid_x_;
    vector<float>  sc_hybrid_y_;
    vector<float>  sc_hybrid_z_;
    vector<float>  sc_hybrid_EtaWidth_;
    vector<float>  sc_hybrid_PhiWidth_;
    vector<float>  sc_hybrid_RawE_;
    vector<float>  sc_hybrid_RawEt_;
    
    //reco::mult55 supercluster
    Int_t          nsc_mult55_;
    vector<float>  sc_mult55_E_;
    vector<float>  sc_mult55_Et_;
    vector<float>  sc_mult55_Eta_;
    vector<float>  sc_mult55_Phi_;
    vector<float>  sc_mult55_x_;
    vector<float>  sc_mult55_y_;
    vector<float>  sc_mult55_z_;
    vector<float>  sc_mult55_EtaWidth_;
    vector<float>  sc_mult55_PhiWidth_;
    vector<float>  sc_mult55_RawE_;
    vector<float>  sc_mult55_RawEt_;
    
    
        
};

#endif
