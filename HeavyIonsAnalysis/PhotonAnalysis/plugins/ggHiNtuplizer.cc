#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "EgammaAnalysis/ElectronTools/interface/SuperClusterHelper.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "HeavyIonsAnalysis/PhotonAnalysis/src/pfIsoCalculator.h"

#include "HeavyIonsAnalysis/PhotonAnalysis/interface/ggHiNtuplizer.h"
#include "HeavyIonsAnalysis/PhotonAnalysis/interface/GenParticleParentage.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"


using namespace std;

ggHiNtuplizer::ggHiNtuplizer(const edm::ParameterSet& ps):
  effectiveAreas_( (ps.getParameter<edm::FileInPath>("effAreasConfigFile")).fullPath() )

{

  // class instance configuration
  doGenParticles_         = ps.getParameter<bool>("doGenParticles");
  runOnParticleGun_       = ps.getParameter<bool>("runOnParticleGun");
  useValMapIso_           = ps.getParameter<bool>("useValMapIso");
  doPfIso_                = ps.getParameter<bool>("doPfIso");
  doVsIso_                = ps.getParameter<bool>("doVsIso");
  genPileupCollection_    = consumes<vector<PileupSummaryInfo> >   (ps.getParameter<edm::InputTag>("pileupCollection"));
  genParticlesCollection_ = consumes<vector<reco::GenParticle> >   (ps.getParameter<edm::InputTag>("genParticleSrc"));
  gsfElectronsCollection_ = consumes<edm::View<reco::GsfElectron> >(ps.getParameter<edm::InputTag>("gsfElectronLabel"));
  input_electrons_        = consumes<reco::GsfElectronCollection>  (ps.getParameter<edm::InputTag>("elecoll"));
  gsfTracks_              = consumes<reco::GsfTrackCollection>     (ps.getParameter<edm::InputTag>("eletrk"));
  genTracks_              = consumes<reco::TrackCollection>        (ps.getParameter<edm::InputTag>("gentrk"));
  recoPhotonsCollection_  = consumes<edm::View<reco::Photon> >     (ps.getParameter<edm::InputTag>("recoPhotonSrc"));
  recoHyPhotonsCollection_ = consumes<edm::View<reco::Photon> >     (ps.getParameter<edm::InputTag>("recoHyPhotonSrc"));
  recoMuonsCollection_    = consumes<edm::View<reco::Muon> >       (ps.getParameter<edm::InputTag>("recoMuonSrc"));
  //CaloTowerCollection_    = consumes<edm::SortedCollection<CaloTower,edm::StrictWeakOrdering<CaloTower>>>(ps.getParameter<edm::InputTag>("recoCaloTower"));
  CaloTowerCollection_    = consumes<edm::SortedCollection<CaloTower>>(ps.getParameter<edm::InputTag>("recoCaloTower"));
  hybridsc_               = consumes<reco::SuperClusterCollection>    (ps.getParameter<edm::InputTag>("hybridsc"));
  mult55sc_               = consumes<reco::SuperClusterCollection>    (ps.getParameter<edm::InputTag>("mult55sc"));
  vtxCollection_          = consumes<vector<reco::Vertex> >        (ps.getParameter<edm::InputTag>("VtxLabel"));
  doVID_                  = ps.getParameter<bool>("doElectronVID");
  //if(doVID_){
  rhoToken_               = consumes<double> (ps.getParameter <edm::InputTag>("rho"));
  beamSpotToken_          = consumes<reco::BeamSpot>(ps.getParameter <edm::InputTag>("beamSpot"));
  conversionsToken_       = consumes< reco::ConversionCollection >(ps.getParameter<edm::InputTag>("conversions"));
  eleVetoIdMapToken_      = consumes<edm::ValueMap<bool> >(ps.getParameter<edm::InputTag>("electronVetoID"));
  eleLooseIdMapToken_     = consumes<edm::ValueMap<bool> >(ps.getParameter<edm::InputTag>("electronLooseID"));
  eleMediumIdMapToken_    = consumes<edm::ValueMap<bool> >(ps.getParameter<edm::InputTag>("electronMediumID"));
  eleTightIdMapToken_     = consumes<edm::ValueMap<bool> >(ps.getParameter<edm::InputTag>("electronTightID"));
  //}
  if(useValMapIso_){
    recoPhotonsHiIso_ = consumes<edm::ValueMap<reco::HIPhotonIsolation> > (ps.getParameter<edm::InputTag>("recoPhotonHiIsolationMap"));
    recoHyPhotonsHiIso_ = consumes<edm::ValueMap<reco::HIPhotonIsolation> > (ps.getParameter<edm::InputTag>("recoHyPhotonHiIsolationMap"));
  }
   if(doPfIso_){
    pfCollection_ = consumes<edm::View<reco::PFCandidate> > (ps.getParameter<edm::InputTag>("particleFlowCollection"));
    if(doVsIso_){
      voronoiBkgCalo_ = consumes<edm::ValueMap<reco::VoronoiBackground> > (ps.getParameter<edm::InputTag>("voronoiBackgroundCalo"));
      voronoiBkgPF_ = consumes<edm::ValueMap<reco::VoronoiBackground> > (ps.getParameter<edm::InputTag>("voronoiBackgroundPF"));
    }
  }
 
 
  
  
  // initialize output TTree
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("EventTree", "Event data");

  tree_->Branch("run",    &run_);
  tree_->Branch("event",  &event_);
  tree_->Branch("lumis",  &lumis_);
  tree_->Branch("isData", &isData_);
  tree_->Branch("nVtx",   &nVtx_);
  tree_->Branch("isFakeVtx",   &isFakeVtx_);
  tree_->Branch("xVtx",   &xVtx_);
  tree_->Branch("yVtx",   &yVtx_);
  tree_->Branch("zVtx",   &zVtx_);
  
  if (doGenParticles_) {
    tree_->Branch("nPUInfo",      &nPUInfo_);
    tree_->Branch("nPU",          &nPU_);
    tree_->Branch("puBX",         &puBX_);
    tree_->Branch("puTrue",       &puTrue_);

    tree_->Branch("nMC",          &nMC_);
    tree_->Branch("mcPID",        &mcPID_);
    tree_->Branch("mcStatus",     &mcStatus_);
    tree_->Branch("mcVtx_x",      &mcVtx_x_);
    tree_->Branch("mcVtx_y",      &mcVtx_y_);
    tree_->Branch("mcVtx_z",      &mcVtx_z_);
    tree_->Branch("mcPt",         &mcPt_);
    tree_->Branch("mcP",          &mcP_);
    tree_->Branch("mcEta",        &mcEta_);
    tree_->Branch("mcPhi",        &mcPhi_);
    tree_->Branch("mcE",          &mcE_);
    tree_->Branch("mcEt",         &mcEt_);
    tree_->Branch("mcMass",       &mcMass_);
    tree_->Branch("mcParentage",  &mcParentage_);
    tree_->Branch("mcMomPID",     &mcMomPID_);
    tree_->Branch("mcMomPt",      &mcMomPt_);
    tree_->Branch("mcMomEta",     &mcMomEta_);
    tree_->Branch("mcMomPhi",     &mcMomPhi_);
    tree_->Branch("mcMomMass",    &mcMomMass_);
    tree_->Branch("mcGMomPID",    &mcGMomPID_);
    tree_->Branch("mcIndex",      &mcIndex_);
    tree_->Branch("mcCalIsoDR03", &mcCalIsoDR03_);
    tree_->Branch("mcCalIsoDR04", &mcCalIsoDR04_);
    tree_->Branch("mcTrkIsoDR03", &mcTrkIsoDR03_);
    tree_->Branch("mcTrkIsoDR04", &mcTrkIsoDR04_);
  }

  tree_->Branch("nEle",                  &nEle_);
  tree_->Branch("eleCharge",             &eleCharge_);
  tree_->Branch("eleChargeConsistent",   &eleChargeConsistent_);
  tree_->Branch("eleSCPixCharge",	 &eleSCPixCharge_);
  tree_->Branch("eleCtfCharge",		 &eleCtfCharge_);
  tree_->Branch("eleEn",                 &eleEn_);
  tree_->Branch("eleD0",                 &eleD0_);
  tree_->Branch("eleDz",                 &eleDz_);
  tree_->Branch("eleD0Err",              &eleD0Err_);
  tree_->Branch("eleDzErr",              &eleDzErr_);
  tree_->Branch("eleTrkPt",              &eleTrkPt_);
  tree_->Branch("eleTrkEta",             &eleTrkEta_);
  tree_->Branch("eleTrkPhi",             &eleTrkPhi_);
  tree_->Branch("eleTrkCharge",          &eleTrkCharge_);
  tree_->Branch("eleTrkChi2",            &eleTrkChi2_);
  tree_->Branch("eleTrkNdof",            &eleTrkNdof_);
  tree_->Branch("eleTrkNormalizedChi2",  &eleTrkNormalizedChi2_);
  tree_->Branch("eleTrkValidHits",       &eleTrkValidHits_);
  tree_->Branch("eleTrkLayers",          &eleTrkLayers_);
  tree_->Branch("eleP",                  &eleP_);
  tree_->Branch("eleP_atVtx",            &eleP_atVtx_);
  tree_->Branch("elePt",                 &elePt_);
  tree_->Branch("eleEta",                &eleEta_);
  tree_->Branch("elePhi",                &elePhi_);
  tree_->Branch("eleSCx",                &eleSCx_);
  tree_->Branch("eleSCy",                &eleSCy_);
  tree_->Branch("eleSCz",                &eleSCz_);
  tree_->Branch("eleSCEn",               &eleSCEn_);
  tree_->Branch("eleESEn",               &eleESEn_);
  tree_->Branch("eleSCEta",              &eleSCEta_);
  tree_->Branch("eleSCPhi",              &eleSCPhi_);
  tree_->Branch("eleSCRawEn",            &eleSCRawEn_);
  tree_->Branch("eleSCEtaWidth",         &eleSCEtaWidth_);
  tree_->Branch("eleSCPhiWidth",         &eleSCPhiWidth_);
  tree_->Branch("eleHoverE",             &eleHoverE_);
  tree_->Branch("eleHoverEBc",           &eleHoverEBc_);
  tree_->Branch("eleEoverP",             &eleEoverP_);
  tree_->Branch("eleEoverPInv",          &eleEoverPInv_);
  tree_->Branch("eleBrem",               &eleBrem_);
  tree_->Branch("eledEtaAtVtx",          &eledEtaAtVtx_);
  tree_->Branch("eledPhiAtVtx",          &eledPhiAtVtx_);
  tree_->Branch("eledEtaSeedAtVtx",      &eledEtaSeedAtVtx_);
  tree_->Branch("eleSigmaIEtaIEta",      &eleSigmaIEtaIEta_);
  tree_->Branch("eleSigmaIEtaIEta_2012", &eleSigmaIEtaIEta_2012_);
  tree_->Branch("eleSigmaIPhiIPhi",      &eleSigmaIPhiIPhi_);
  // tree_->Branch("eleConvVeto",           &eleConvVeto_);  // TODO: not available in reco::
  tree_->Branch("eleMissHits",           &eleMissHits_);
  tree_->Branch("eleTrackIso",           &eleTrackIso_);
  tree_->Branch("eleECalIso",            &eleECalIso_);
  tree_->Branch("eleHCalIso",            &eleHCalIso_);
  tree_->Branch("eleECalDriven",         &eleECalDriven_);
  tree_->Branch("eleESEffSigmaRR",       &eleESEffSigmaRR_);
  tree_->Branch("elePFChIso",            &elePFChIso_);
  tree_->Branch("elePFPhoIso",           &elePFPhoIso_);
  tree_->Branch("elePFNeuIso",           &elePFNeuIso_);
  tree_->Branch("elePFPUIso",            &elePFPUIso_);
  tree_->Branch("elePFChIso03",          &elePFChIso03_);
  tree_->Branch("elePFPhoIso03",         &elePFPhoIso03_);
  tree_->Branch("elePFNeuIso03",         &elePFNeuIso03_);
  tree_->Branch("elePFChIso04",          &elePFChIso04_);
  tree_->Branch("elePFPhoIso04",         &elePFPhoIso04_);
  tree_->Branch("elePFNeuIso04",         &elePFNeuIso04_);

  tree_->Branch("eleR9",		 &eleR9_);
  tree_->Branch("eleE3x3",		 &eleE3x3_);
  tree_->Branch("eleE5x5",		 &eleE5x5_);
  tree_->Branch("eleR9Full5x5",		 &eleR9Full5x5_);
  tree_->Branch("eleE3x3Full5x5",	 &eleE3x3Full5x5_);
  tree_->Branch("eleE5x5Full5x5",	 &eleE5x5Full5x5_);
  tree_->Branch("NClusters",		 &NClusters_);
  tree_->Branch("NEcalClusters",	 &NEcalClusters_);
  tree_->Branch("eleSeedEn",		 &eleSeedEn_);
  tree_->Branch("eleSeedEta",		 &eleSeedEta_);
  tree_->Branch("eleSeedPhi",		 &eleSeedPhi_);
  tree_->Branch("eleSeedCryEta",	 &eleSeedCryEta_);
  tree_->Branch("eleSeedCryPhi",	 &eleSeedCryPhi_);
  tree_->Branch("eleSeedCryIeta",	 &eleSeedCryIeta_);
  tree_->Branch("eleSeedCryIphi",	 &eleSeedCryIphi_);

  tree_->Branch("eleBC1E",               &eleBC1E_);
  tree_->Branch("eleBC1Eta",             &eleBC1Eta_);
  tree_->Branch("eleBC2E",               &eleBC2E_);
  tree_->Branch("eleBC2Eta",             &eleBC2Eta_);
  tree_->Branch("eleIDVeto",             &eleIDVeto_);
  tree_->Branch("eleIDLoose",            &eleIDLoose_);
  tree_->Branch("eleIDMedium",           &eleIDMedium_);
  tree_->Branch("eleIDTight",            &eleIDTight_);
  tree_->Branch("elepassConversionVeto", &elepassConversionVeto_);
  tree_->Branch("eleEffAreaTimesRho",    &eleEffAreaTimesRho_);

  tree_->Branch("nConv",                 &nConv_);
  tree_->Branch("conv_eleind",           &conv_eleind_);
  tree_->Branch("conv_phoind",           &conv_phoind_);
  tree_->Branch("conv_vtxprob",          &conv_vtxprob_);
  tree_->Branch("conv_lxy",              &conv_lxy_);
  tree_->Branch("conv_nhits_bvtx",       &conv_nhits_bvtx_);

  // ged photon branches
  tree_->Branch("nPho",                  &nPho_);
  tree_->Branch("phoE",                  &phoE_);
  tree_->Branch("phoEt",                 &phoEt_);
  tree_->Branch("phoEta",                &phoEta_);
  tree_->Branch("phoPhi",                &phoPhi_);
  tree_->Branch("phoSCSize",             &phoSCSize_);
  tree_->Branch("phoSCE",                &phoSCE_);
  tree_->Branch("phoSCEt",               &phoSCEt_);
  tree_->Branch("phoSCRawE",             &phoSCRawE_);
  tree_->Branch("phoSCRawEt",            &phoSCRawEt_);
  tree_->Branch("phoESEn",               &phoESEn_);
  tree_->Branch("phoSCEta",              &phoSCEta_);
  tree_->Branch("phoSCPhi",              &phoSCPhi_);
  tree_->Branch("phoSCx",                &phoSCx_);
  tree_->Branch("phoSCy",                &phoSCy_);
  tree_->Branch("phoSCz",                &phoSCz_);
  tree_->Branch("phoSCEtaWidth",         &phoSCEtaWidth_);
  tree_->Branch("phoSCPhiWidth",         &phoSCPhiWidth_);
  tree_->Branch("phoSCBrem",             &phoSCBrem_);
  tree_->Branch("phohasPixelSeed",       &phohasPixelSeed_);
  tree_->Branch("phopassConversionVeto", &phopassConversionVeto_);
  //tree_->Branch("phoEleVeto",            &phoEleVeto_);        // TODO: not available in reco::
  tree_->Branch("phoR9",                 &phoR9_);
  tree_->Branch("phoHadTowerOverEm",     &phoHadTowerOverEm_);
  tree_->Branch("phoHoverE",             &phoHoverE_);
  tree_->Branch("phoSigmaIEtaIEta",      &phoSigmaIEtaIEta_);
  // tree_->Branch("phoSigmaIEtaIPhi",      &phoSigmaIEtaIPhi_);  // TODO: not available in reco::
  // tree_->Branch("phoSigmaIPhiIPhi",      &phoSigmaIPhiIPhi_);  // TODO: not available in reco::
  tree_->Branch("phoE1x3",               &phoE1x3_);
  tree_->Branch("phoE2x2",               &phoE2x2_);
  tree_->Branch("phoE3x3",               &phoE3x3_);
  tree_->Branch("phoE2x5Max",            &phoE2x5Max_);
  tree_->Branch("phoE1x5",               &phoE1x5_);
  tree_->Branch("phoE2x5",               &phoE2x5_);
  tree_->Branch("phoE5x5",               &phoE5x5_);
  tree_->Branch("phoMaxEnergyXtal",      &phoMaxEnergyXtal_);
  tree_->Branch("phoSigmaEtaEta",        &phoSigmaEtaEta_);
  tree_->Branch("phoR1x5",               &phoR1x5_);
  tree_->Branch("phoR2x5",               &phoR2x5_);
  tree_->Branch("phoESEffSigmaRR",       &phoESEffSigmaRR_);
  tree_->Branch("phoSigmaIEtaIEta_2012", &phoSigmaIEtaIEta_2012_);
  tree_->Branch("phoSigmaIEtaIPhi_2012", &phoSigmaIEtaIPhi_2012_);
  tree_->Branch("phoSigmaIPhiIPhi_2012", &phoSigmaIPhiIPhi_2012_);
  tree_->Branch("phoE1x3_2012",          &phoE1x3_2012_);
  tree_->Branch("phoE2x2_2012",          &phoE2x2_2012_);
  tree_->Branch("phoE3x3_2012",          &phoE3x3_2012_);
  tree_->Branch("phoE2x5Max_2012",       &phoE2x5Max_2012_);
  tree_->Branch("phoE5x5_2012",          &phoE5x5_2012_);
  tree_->Branch("phoBC1E",               &phoBC1E_);
  tree_->Branch("phoBC1Eta",             &phoBC1Eta_);
  tree_->Branch("phoBC2E",               &phoBC2E_);
  tree_->Branch("phoBC2Eta",             &phoBC2Eta_);
  tree_->Branch("pho_ecalClusterIsoR2", &pho_ecalClusterIsoR2_);
  tree_->Branch("pho_ecalClusterIsoR3", &pho_ecalClusterIsoR3_);
  tree_->Branch("pho_ecalClusterIsoR4", &pho_ecalClusterIsoR4_);
  tree_->Branch("pho_ecalClusterIsoR5", &pho_ecalClusterIsoR5_);
  tree_->Branch("pho_hcalRechitIsoR1", &pho_hcalRechitIsoR1_);
  tree_->Branch("pho_hcalRechitIsoR2", &pho_hcalRechitIsoR2_);
  tree_->Branch("pho_hcalRechitIsoR3", &pho_hcalRechitIsoR3_);
  tree_->Branch("pho_hcalRechitIsoR4", &pho_hcalRechitIsoR4_);
  tree_->Branch("pho_hcalRechitIsoR5", &pho_hcalRechitIsoR5_);
  tree_->Branch("pho_trackIsoR1PtCut20", &pho_trackIsoR1PtCut20_);
  tree_->Branch("pho_trackIsoR2PtCut20", &pho_trackIsoR2PtCut20_);
  tree_->Branch("pho_trackIsoR3PtCut20", &pho_trackIsoR3PtCut20_);
  tree_->Branch("pho_trackIsoR4PtCut20", &pho_trackIsoR4PtCut20_);
  tree_->Branch("pho_trackIsoR5PtCut20", &pho_trackIsoR5PtCut20_);
  tree_->Branch("pho_swissCrx", &pho_swissCrx_);
  tree_->Branch("pho_seedTime", &pho_seedTime_);
  if (doGenParticles_) {
    tree_->Branch("pho_genMatchedIndex", &pho_genMatchedIndex_);
  }

 
  //ged photon conversion branches 
  tree_->Branch("phoIsConv", &phoIsConv_);
  tree_->Branch("phoNConv", &phoNConv_);
  tree_->Branch("phoConvInvMass", &phoConvInvMass_);
  tree_->Branch("phoConvCotTheta", &phoConvCotTheta_);
  tree_->Branch("phoConvEoverP", &phoConvEoverP_);
  tree_->Branch("phoConvZofPVfromTrks", &phoConvZofPVfromTrks_);
  tree_->Branch("phoConvMinDist", &phoConvMinDist_);
  tree_->Branch("phoConvdPhiAtVtx", &phoConvdPhiAtVtx_);
  tree_->Branch("phoConvdPhiAtCalo", &phoConvdPhiAtCalo_);
  tree_->Branch("phoConvdEtaAtCalo", &phoConvdEtaAtCalo_);
  tree_->Branch("phoConvTrkd0_x", &phoConvTrkd0_x_);
  tree_->Branch("phoConvTrkd0_y", &phoConvTrkd0_y_);
  tree_->Branch("phoConvTrkPin_x", &phoConvTrkPin_x_);
  tree_->Branch("phoConvTrkPin_y", &phoConvTrkPin_y_);
  tree_->Branch("phoConvTrkPout_x", &phoConvTrkPout_x_);
  tree_->Branch("phoConvTrkPout_y", &phoConvTrkPout_y_);
  tree_->Branch("phoConvTrkdz_x", &phoConvTrkdz_x_);
  tree_->Branch("phoConvTrkdz_y", &phoConvTrkdz_y_);
  tree_->Branch("phoConvTrkdzErr_x", &phoConvTrkdzErr_x_);
  tree_->Branch("phoConvTrkdzErr_y", &phoConvTrkdzErr_y_);
  tree_->Branch("phoConvChi2", &phoConvChi2_);
  tree_->Branch("phoConvChi2Prob", &phoConvChi2Prob_);
  tree_->Branch("phoConvNTrks", &phoConvNTrks_);
  tree_->Branch("phoConvCharge1", &phoConvCharge1_);
  tree_->Branch("phoConvCharge2", &phoConvCharge2_);
  tree_->Branch("phoConvValidVtx", &phoConvValidVtx_);
  tree_->Branch("phoConvLikeLihood", &phoConvLikeLihood_);
  tree_->Branch("phoConvP4_0", &phoConvP4_0_);
  tree_->Branch("phoConvP4_1", &phoConvP4_1_);
  tree_->Branch("phoConvP4_2", &phoConvP4_2_);
  tree_->Branch("phoConvP4_3", &phoConvP4_3_);
  tree_->Branch("phoConvVtx_x", &phoConvVtx_x_);
  tree_->Branch("phoConvVtx_y", &phoConvVtx_y_);
  tree_->Branch("phoConvVtx_z", &phoConvVtx_z_);
  tree_->Branch("phoConvVtxErr_x", &phoConvVtxErr_x_);
  tree_->Branch("phoConvVtxErr_y", &phoConvVtxErr_y_);
  tree_->Branch("phoConvVtxErr_z", &phoConvVtxErr_z_);
  tree_->Branch("phoConvPairMomentum_x", &phoConvPairMomentum_x_);
  tree_->Branch("phoConvPairMomentum_y", &phoConvPairMomentum_y_);
  tree_->Branch("phoConvPairMomentum_z", &phoConvPairMomentum_z_);
  tree_->Branch("phoConvRefittedMomentum_x", &phoConvRefittedMomentum_x_);
  tree_->Branch("phoConvRefittedMomentum_y", &phoConvRefittedMomentum_y_);
  tree_->Branch("phoConvRefittedMomentum_z", &phoConvRefittedMomentum_z_);
  tree_->Branch("SingleLegConv", &SingleLegConv_);

  //muon branches
  tree_->Branch("nMu",                   &nMu_);
  tree_->Branch("muPt",                  &muPt_);
  tree_->Branch("muEta",                 &muEta_);
  tree_->Branch("muPhi",                 &muPhi_);
  tree_->Branch("muCharge",              &muCharge_);
  tree_->Branch("muType",                &muType_);
  tree_->Branch("muIsGood",              &muIsGood_);
  tree_->Branch("muD0",                  &muD0_);
  tree_->Branch("muDz",                  &muDz_);
  tree_->Branch("muChi2NDF",             &muChi2NDF_);
  tree_->Branch("muInnerD0",             &muInnerD0_);
  tree_->Branch("muInnerDz",             &muInnerDz_);
  tree_->Branch("muTrkLayers",           &muTrkLayers_);
  tree_->Branch("muPixelLayers",         &muPixelLayers_);
  tree_->Branch("muPixelHits",           &muPixelHits_);
  tree_->Branch("muMuonHits",            &muMuonHits_);
  tree_->Branch("muTrkQuality",          &muTrkQuality_);
  tree_->Branch("muStations",            &muStations_);
  tree_->Branch("muIsoTrk",              &muIsoTrk_);
  tree_->Branch("muPFChIso",             &muPFChIso_);
  tree_->Branch("muPFPhoIso",            &muPFPhoIso_);
  tree_->Branch("muPFNeuIso",            &muPFNeuIso_);
  tree_->Branch("muPFPUIso",             &muPFPUIso_);

  tree_->Branch("nTower",                &nTower_);
  tree_->Branch("CaloTower_hadE",        &CaloTower_hadE_);
  tree_->Branch("CaloTower_emE",         &CaloTower_emE_);
  tree_->Branch("CaloTower_e",           &CaloTower_e_);
  tree_->Branch("CaloTower_et",          &CaloTower_et_);
  tree_->Branch("CaloTower_eta",         &CaloTower_eta_);
  tree_->Branch("CaloTower_phi",         &CaloTower_phi_);

  tree_->Branch("nHyPho",                &nHyPho_);
  tree_->Branch("hyphoE",                  &hyphoE_);
  tree_->Branch("hyphoEt",                 &hyphoEt_);
  tree_->Branch("hyphoEta",                &hyphoEta_);
  tree_->Branch("hyphoPhi",                &hyphoPhi_);
  tree_->Branch("hyphoSCSize",             &hyphoSCSize_);
  tree_->Branch("hyphoSCE",                &hyphoSCE_);
  tree_->Branch("hyphoSCEt",               &hyphoSCEt_);
  tree_->Branch("hyphoSCRawE",             &hyphoSCRawE_);
  tree_->Branch("hyphoSCRawEt",            &hyphoSCRawEt_);
  tree_->Branch("hyphoESEn",               &hyphoESEn_);
  tree_->Branch("hyphoSCEta",              &hyphoSCEta_);
  tree_->Branch("hyphoSCPhi",              &hyphoSCPhi_);
  tree_->Branch("hyphoSCx",                &hyphoSCx_);
  tree_->Branch("hyphoSCy",                &hyphoSCy_);
  tree_->Branch("hyphoSCz",                &hyphoSCz_);
  tree_->Branch("hyphoSCEtaWidth",         &hyphoSCEtaWidth_);
  tree_->Branch("hyphoSCPhiWidth",         &hyphoSCPhiWidth_);
  tree_->Branch("hyphoSCBrem",             &hyphoSCBrem_);
  tree_->Branch("hyphohasPixelSeed",       &hyphohasPixelSeed_);
  tree_->Branch("hyphopassConversionVeto", &hyphopassConversionVeto_);
  //tree_->Branch("hyphoEleVeto",            &hyphoEleVeto_);        // TODO: not available in reco::
  tree_->Branch("hyphoR9",                 &hyphoR9_);
  tree_->Branch("hyphoHadTowerOverEm",     &hyphoHadTowerOverEm_);
  tree_->Branch("hyphoHoverE",             &hyphoHoverE_);
  tree_->Branch("hyphoSigmaIEtaIEta",      &hyphoSigmaIEtaIEta_);
  // tree_->Branch("hyphoSigmaIEtaIPhi",      &hyphoSigmaIEtaIPhi_);  // TODO: not available in reco::
  // tree_->Branch("hyphoSigmaIPhiIPhi",      &hyphoSigmaIPhiIPhi_);  // TODO: not available in reco::
  tree_->Branch("hyphoE1x3",               &hyphoE1x3_);
  tree_->Branch("hyphoE2x2",               &hyphoE2x2_);
  tree_->Branch("hyphoE3x3",               &hyphoE3x3_);
  tree_->Branch("hyphoE2x5Max",            &hyphoE2x5Max_);
  tree_->Branch("hyphoE1x5",               &hyphoE1x5_);
  tree_->Branch("hyphoE2x5",               &hyphoE2x5_);
  tree_->Branch("hyphoE5x5",               &hyphoE5x5_);
  tree_->Branch("hyphoMaxEnergyXtal",      &hyphoMaxEnergyXtal_);
  tree_->Branch("hyphoSigmaEtaEta",        &hyphoSigmaEtaEta_);
  tree_->Branch("hyphoR1x5",               &hyphoR1x5_);
  tree_->Branch("hyphoR2x5",               &hyphoR2x5_);
  tree_->Branch("hyphoESEffSigmaRR",       &hyphoESEffSigmaRR_);
  tree_->Branch("hyphoSigmaIEtaIEta_2012", &hyphoSigmaIEtaIEta_2012_);
  tree_->Branch("hyphoSigmaIEtaIPhi_2012", &hyphoSigmaIEtaIPhi_2012_);
  tree_->Branch("hyphoSigmaIPhiIPhi_2012", &hyphoSigmaIPhiIPhi_2012_);
  tree_->Branch("hyphoE1x3_2012",          &hyphoE1x3_2012_);
  tree_->Branch("hyphoE2x2_2012",          &hyphoE2x2_2012_);
  tree_->Branch("hyphoE3x3_2012",          &hyphoE3x3_2012_);
  tree_->Branch("hyphoE2x5Max_2012",       &hyphoE2x5Max_2012_);
  tree_->Branch("hyphoE5x5_2012",          &hyphoE5x5_2012_);
  tree_->Branch("hyphoBC1E",               &hyphoBC1E_);
  tree_->Branch("hyphoBC1Eta",             &hyphoBC1Eta_);
  tree_->Branch("hyphoBC2E",               &hyphoBC2E_);
  tree_->Branch("hyphoBC2Eta",             &hyphoBC2Eta_);
  tree_->Branch("hypho_ecalClusterIsoR2", &hypho_ecalClusterIsoR2_);
  tree_->Branch("hypho_ecalClusterIsoR3", &hypho_ecalClusterIsoR3_);
  tree_->Branch("hypho_ecalClusterIsoR4", &hypho_ecalClusterIsoR4_);
  tree_->Branch("hypho_ecalClusterIsoR5", &hypho_ecalClusterIsoR5_);
  tree_->Branch("hypho_hcalRechitIsoR1", &hypho_hcalRechitIsoR1_);
  tree_->Branch("hypho_hcalRechitIsoR2", &hypho_hcalRechitIsoR2_);
  tree_->Branch("hypho_hcalRechitIsoR3", &hypho_hcalRechitIsoR3_);
  tree_->Branch("hypho_hcalRechitIsoR4", &hypho_hcalRechitIsoR4_);
  tree_->Branch("hypho_hcalRechitIsoR5", &hypho_hcalRechitIsoR5_);
  tree_->Branch("hypho_trackIsoR1PtCut20", &hypho_trackIsoR1PtCut20_);
  tree_->Branch("hypho_trackIsoR2PtCut20", &hypho_trackIsoR2PtCut20_);
  tree_->Branch("hypho_trackIsoR3PtCut20", &hypho_trackIsoR3PtCut20_);
  tree_->Branch("hypho_trackIsoR4PtCut20", &hypho_trackIsoR4PtCut20_);
  tree_->Branch("hypho_trackIsoR5PtCut20", &hypho_trackIsoR5PtCut20_);
  tree_->Branch("hypho_swissCrx", &hypho_swissCrx_);
  tree_->Branch("hypho_seedTime", &hypho_seedTime_);
  if (doGenParticles_) {
    tree_->Branch("hypho_genMatchedIndex", &hypho_genMatchedIndex_);
  }

  tree_->Branch("ngsfEle",               &ngsfEle_);
  tree_->Branch("elegsfTrkPt",           &elegsfTrkPt_);
  tree_->Branch("elegsfTrkP",            &elegsfTrkP_);
  tree_->Branch("elegsfTrkEta",          &elegsfTrkEta_);
  tree_->Branch("elegsfTrkPhi",          &elegsfTrkPhi_);
  tree_->Branch("elegsfTrkCharge",       &elegsfTrkCharge_);
  tree_->Branch("elegsfTrkChi2",         &elegsfTrkChi2_);
  tree_->Branch("elegsfTrkNdof",         &elegsfTrkNdof_);
  tree_->Branch("elegsfTrkNormalizedChi2",&elegsfTrkNormalizedChi2_);
  tree_->Branch("elegsfTrkValidHits",    &elegsfTrkValidHits_);
  tree_->Branch("elegsfTrkMissHits",     &elegsfTrkMissHits_);
  tree_->Branch("elegsfTrkLayers",       &elegsfTrkLayers_);
  tree_->Branch("elegsfD0",              &elegsfD0_);
  tree_->Branch("elegsfDz",              &elegsfDz_);
  tree_->Branch("elegsfD0Err",           &elegsfD0Err_);
  tree_->Branch("elegsfDzErr",           &elegsfDzErr_);

  tree_->Branch("ngenTrk",            &ngenTrk_);
  tree_->Branch("gentrkPt",           &gentrkPt_);
  tree_->Branch("gentrkP",            &gentrkP_);
  tree_->Branch("gentrkEta",          &gentrkEta_);
  tree_->Branch("gentrkPhi",          &gentrkPhi_);
  tree_->Branch("gentrkcharge",       &gentrkcharge_);
  tree_->Branch("gentrkvx",           &gentrkvx_);
  tree_->Branch("gentrkvy",           &gentrkvy_);
  tree_->Branch("gentrkvz",           &gentrkvz_);
  tree_->Branch("gentrknormchi2",     &gentrknormchi2_);
  tree_->Branch("gentrkchi2",         &gentrkchi2_);
  tree_->Branch("gentrkd0",           &gentrkd0_);
  tree_->Branch("gentrkdxy",          &gentrkdxy_);
  tree_->Branch("gentrkdz",           &gentrkdz_);
  tree_->Branch("gentrkdxyError",     &gentrkdxyError_);
  tree_->Branch("gentrkdzError",      &gentrkdzError_);
  tree_->Branch("gentrkValidHits",    &gentrkValidHits_);
  tree_->Branch("gentrkMissHits",     &gentrkMissHits_);
  tree_->Branch("gentrkPurity",       &gentrkPurity_);
  
  
  // hybrid cluster collection
  tree_->Branch("nsc_hybrid",         &nsc_hybrid_);
  tree_->Branch("sc_hybrid_E",        &sc_hybrid_E_);
  tree_->Branch("sc_hybrid_Et",       &sc_hybrid_Et_);
  tree_->Branch("sc_hybrid_Eta",      &sc_hybrid_Eta_);
  tree_->Branch("sc_hybrid_Phi",      &sc_hybrid_Phi_);
  tree_->Branch("sc_hybrid_x",        &sc_hybrid_x_);
  tree_->Branch("sc_hybrid_y",        &sc_hybrid_y_);
  tree_->Branch("sc_hybrid_z",        &sc_hybrid_z_);
  tree_->Branch("sc_hybrid_EtaWidth", &sc_hybrid_EtaWidth_);
  tree_->Branch("sc_hybrid_PhiWidth", &sc_hybrid_PhiWidth_);
  tree_->Branch("sc_hybrid_RawE",     &sc_hybrid_RawE_);
  tree_->Branch("sc_hybrid_RawEt",    &sc_hybrid_RawEt_);
  
  // mult55 cluster collection
  tree_->Branch("nsc_mult55",         &nsc_mult55_);
  tree_->Branch("sc_mult55_E",        &sc_mult55_E_);
  tree_->Branch("sc_mult55_Et",       &sc_mult55_Et_);
  tree_->Branch("sc_mult55_Eta",      &sc_mult55_Eta_);
  tree_->Branch("sc_mult55_Phi",      &sc_mult55_Phi_);
  tree_->Branch("sc_mult55_x",        &sc_mult55_x_);
  tree_->Branch("sc_mult55_y",        &sc_mult55_y_);
  tree_->Branch("sc_mult55_z",        &sc_mult55_z_);
  tree_->Branch("sc_mult55_EtaWidth", &sc_mult55_EtaWidth_);
  tree_->Branch("sc_mult55_PhiWidth", &sc_mult55_PhiWidth_);
  tree_->Branch("sc_mult55_RawE",     &sc_mult55_RawE_);
  tree_->Branch("sc_mult55_RawEt",    &sc_mult55_RawEt_);
  
  
}

void ggHiNtuplizer::analyze(const edm::Event& e, const edm::EventSetup& es)
{
  
  // cleanup from previous event
  nVtx_ = 0;
  nPUInfo_ = 0;
  nMC_ = 0;
  nEle_ = 0;
  nPho_ = 0;
  nHyPho_ = 0;
  nConv_ = 0;
  nMu_ = 0;
  nTower_= 0;
  ngsfEle_ = 0;
  ngenTrk_ = 0;
  
  
  isFakeVtx_            .clear();
  xVtx_                 .clear();
  yVtx_                 .clear();
  zVtx_                 .clear();
  
  nPU_                  .clear();
  puBX_                 .clear();
  puTrue_               .clear();

  mcPID_                .clear();
  mcStatus_             .clear();
  mcVtx_x_              .clear();
  mcVtx_y_              .clear();
  mcVtx_z_              .clear();
  mcPt_                 .clear();
  mcP_                  .clear();
  mcEta_                .clear();
  mcPhi_                .clear();
  mcE_                  .clear();
  mcEt_                 .clear();
  mcMass_               .clear();
  mcParentage_          .clear();
  mcMomPID_             .clear();
  mcMomPt_              .clear();
  mcMomEta_             .clear();
  mcMomPhi_             .clear();
  mcMomMass_            .clear();
  mcGMomPID_            .clear();
  mcIndex_              .clear();
  mcCalIsoDR03_         .clear();
  mcCalIsoDR04_         .clear();
  mcTrkIsoDR03_         .clear();
  mcTrkIsoDR04_         .clear();

  eleCharge_            .clear();
  eleChargeConsistent_  .clear();
  eleSCPixCharge_	.clear();
  eleCtfCharge_		.clear();
  eleEn_                .clear();
  eleD0_                .clear();
  eleDz_                .clear();
  eleD0Err_             .clear();
  eleDzErr_             .clear();
  eleTrkPt_             .clear();
  eleTrkEta_            .clear();
  eleTrkPhi_            .clear();
  eleTrkCharge_         .clear();
  eleTrkChi2_           .clear();
  eleTrkNdof_           .clear();
  eleTrkNormalizedChi2_ .clear();
  eleTrkValidHits_      .clear();
  eleTrkLayers_         .clear();
  eleP_                 .clear();
  eleP_atVtx_           .clear();
  elePt_                .clear();
  eleEta_               .clear();
  elePhi_               .clear();
  eleSCx_               .clear();
  eleSCy_               .clear();
  eleSCz_               .clear();
  eleSCEn_              .clear();
  eleESEn_              .clear();
  eleSCEta_             .clear();
  eleSCPhi_             .clear();
  eleSCRawEn_           .clear();
  eleSCEtaWidth_        .clear();
  eleSCPhiWidth_        .clear();
  eleHoverE_            .clear();
  eleHoverEBc_          .clear();
  eleEoverP_            .clear();
  eleEoverPInv_         .clear();
  eleBrem_              .clear();
  eledEtaAtVtx_         .clear();
  eledPhiAtVtx_         .clear();
  eledEtaSeedAtVtx_     .clear();
  eleSigmaIEtaIEta_     .clear();
  eleSigmaIEtaIEta_2012_.clear();
  eleSigmaIPhiIPhi_     .clear();
  // eleConvVeto_          .clear();  // TODO: not available in reco::
  eleMissHits_          .clear();
  eleTrackIso_          .clear();
  eleECalIso_           .clear();
  eleHCalIso_           .clear();
  eleECalDriven_        .clear();
  eleESEffSigmaRR_      .clear();
  elePFChIso_           .clear();
  elePFPhoIso_          .clear();
  elePFNeuIso_          .clear();
  elePFPUIso_           .clear();
  elePFChIso03_         .clear();
  elePFPhoIso03_        .clear();
  elePFNeuIso03_        .clear();
  elePFChIso04_         .clear();
  elePFPhoIso04_        .clear();
  elePFNeuIso04_        .clear();
  eleR9_		.clear();
  eleE3x3_		.clear();
  eleE5x5_		.clear();
  eleR9Full5x5_		.clear();
  eleE3x3Full5x5_	.clear();
  eleE5x5Full5x5_	.clear();
  NClusters_		.clear();
  NEcalClusters_	.clear();
  eleSeedEn_		.clear();
  eleSeedEta_		.clear();
  eleSeedPhi_		.clear();
  eleSeedCryEta_	.clear();
  eleSeedCryPhi_	.clear();
  eleSeedCryIeta_	.clear();
  eleSeedCryIphi_	.clear();
  eleBC1E_              .clear();
  eleBC1Eta_            .clear();
  eleBC2E_              .clear();
  eleBC2Eta_            .clear();
  eleIDVeto_            .clear();
  eleIDLoose_           .clear();
  eleIDMedium_          .clear();
  eleIDTight_           .clear();
  elepassConversionVeto_.clear();
  eleEffAreaTimesRho_   .clear();
 
  conv_vtxprob_         .clear();
  conv_lxy_             .clear();
  conv_nhits_bvtx_      .clear();

  phoE_                 .clear();
  phoEt_                .clear();
  phoEta_               .clear();
  phoPhi_               .clear();
  phoSCSize_            .clear();
  phoSCE_               .clear();
  phoSCEt_              .clear();
  phoSCRawE_            .clear();
  phoSCRawEt_           .clear();
  phoESEn_              .clear();
  phoSCEta_             .clear();
  phoSCPhi_             .clear();
  phoSCx_               .clear();
  phoSCy_               .clear();
  phoSCz_               .clear();
  phoSCEtaWidth_        .clear();
  phoSCPhiWidth_        .clear();
  phoSCBrem_            .clear();
  phohasPixelSeed_      .clear();
  phopassConversionVeto_.clear();
  //phoEleVeto_           .clear();  // TODO: not available in reco::
  phoR9_                .clear();
  phoHadTowerOverEm_    .clear();
  phoHoverE_            .clear();
  phoSigmaIEtaIEta_     .clear();
  // phoSigmaIEtaIPhi_     .clear();  // TODO: not available in reco::
  // phoSigmaIPhiIPhi_     .clear();  // TODO: not available in reco::
  phoE1x3_              .clear();
  phoE2x2_              .clear();
  phoE3x3_              .clear();
  phoE2x5Max_           .clear();
  phoE1x5_              .clear();
  phoE2x5_              .clear();
  phoE5x5_              .clear();
  phoMaxEnergyXtal_     .clear();
  phoSigmaEtaEta_       .clear();
  phoR1x5_              .clear();
  phoR2x5_              .clear();
  phoESEffSigmaRR_      .clear();
  phoSigmaIEtaIEta_2012_.clear();
  phoSigmaIEtaIPhi_2012_.clear();
  phoSigmaIPhiIPhi_2012_.clear();
  phoE1x3_2012_         .clear();
  phoE2x2_2012_         .clear();
  phoE3x3_2012_         .clear();
  phoE2x5Max_2012_      .clear();
  phoE5x5_2012_         .clear();
  phoBC1E_              .clear();
  phoBC1Eta_            .clear();
  phoBC2E_              .clear();
  phoBC2Eta_            .clear();
  pho_ecalClusterIsoR2_.clear();
  pho_ecalClusterIsoR3_.clear();
  pho_ecalClusterIsoR4_.clear();
  pho_ecalClusterIsoR5_.clear();
  pho_hcalRechitIsoR1_.clear();
  pho_hcalRechitIsoR2_.clear();
  pho_hcalRechitIsoR3_.clear();
  pho_hcalRechitIsoR4_.clear();
  pho_hcalRechitIsoR5_.clear();
  pho_trackIsoR1PtCut20_.clear();
  pho_trackIsoR2PtCut20_.clear();
  pho_trackIsoR3PtCut20_.clear();
  pho_trackIsoR4PtCut20_.clear();
  pho_trackIsoR5PtCut20_.clear();
  pho_swissCrx_.clear();
  pho_seedTime_.clear();
  pho_genMatchedIndex_.clear();

  
  phoIsConv_.clear();
  phoNConv_.clear();
  phoConvNTrks_.clear();
  phoConvInvMass_.clear();
  phoConvCotTheta_.clear();
  phoConvEoverP_.clear();
  phoConvZofPVfromTrks_.clear();
  phoConvMinDist_.clear();
  phoConvdPhiAtVtx_.clear();
  phoConvdPhiAtCalo_.clear();
  phoConvdEtaAtCalo_.clear();
  phoConvTrkd0_x_.clear();
  phoConvTrkd0_y_.clear();
  phoConvTrkPin_x_.clear();
  phoConvTrkPin_y_.clear();
  phoConvTrkPout_x_.clear();
  phoConvTrkPout_y_.clear();
  phoConvTrkdz_x_.clear();
  phoConvTrkdz_y_.clear();
  phoConvTrkdzErr_x_.clear();
  phoConvTrkdzErr_y_.clear();
  phoConvChi2_.clear();
  phoConvChi2Prob_.clear();
  phoConvCharge1_.clear();
  phoConvCharge2_.clear();
  phoConvValidVtx_.clear();
  phoConvLikeLihood_.clear();
  phoConvP4_0_.clear();
  phoConvP4_1_.clear();
  phoConvP4_2_.clear();
  phoConvP4_3_.clear();
  phoConvVtx_x_.clear();
  phoConvVtx_y_.clear();
  phoConvVtx_z_.clear();
  phoConvVtxErr_x_.clear();
  phoConvVtxErr_y_.clear();
  phoConvVtxErr_z_.clear();
  phoConvPairMomentum_x_.clear();
  phoConvPairMomentum_y_.clear();
  phoConvPairMomentum_z_.clear();
  phoConvRefittedMomentum_x_.clear();
  phoConvRefittedMomentum_y_.clear();
  phoConvRefittedMomentum_z_.clear();
  SingleLegConv_.clear();


  muPt_                 .clear();
  muEta_                .clear();
  muPhi_                .clear();
  muCharge_             .clear();
  muType_               .clear();
  muIsGood_             .clear();
  muD0_                 .clear();
  muDz_                 .clear();
  muChi2NDF_            .clear();
  muInnerD0_            .clear();
  muInnerDz_            .clear();
  muTrkLayers_          .clear();
  muPixelLayers_        .clear();
  muPixelHits_          .clear();
  muMuonHits_           .clear();
  muTrkQuality_         .clear();
  muStations_           .clear();
  muIsoTrk_             .clear();
  muPFChIso_            .clear();
  muPFPhoIso_           .clear();
  muPFNeuIso_           .clear();
  muPFPUIso_            .clear();
  
  CaloTower_hadE_       .clear();
  CaloTower_emE_        .clear();
  CaloTower_e_          .clear();
  CaloTower_et_         .clear();
  CaloTower_eta_        .clear();
  CaloTower_phi_        .clear();
  

  hyphoE_                 .clear();
  hyphoEt_                .clear();
  hyphoEta_               .clear();
  hyphoPhi_               .clear();
  hyphoSCSize_            .clear();
  hyphoSCE_               .clear();
  hyphoSCEt_              .clear();
  hyphoSCRawE_            .clear();
  hyphoSCRawEt_           .clear();
  hyphoESEn_              .clear();
  hyphoSCEta_             .clear();
  hyphoSCPhi_             .clear();
  hyphoSCx_               .clear();
  hyphoSCy_               .clear();
  hyphoSCz_               .clear();
  hyphoSCEtaWidth_        .clear();
  hyphoSCPhiWidth_        .clear();
  hyphoSCBrem_            .clear();
  hyphohasPixelSeed_      .clear();
  hyphopassConversionVeto_.clear();
  //hyphoEleVeto_           .clear();  // TODO: not available in reco::
  hyphoR9_                .clear();
  hyphoHadTowerOverEm_    .clear();
  hyphoHoverE_            .clear();
  hyphoSigmaIEtaIEta_     .clear();
  // hyphoSigmaIEtaIPhi_     .clear();  // TODO: not available in reco::
  // hyphoSigmaIPhiIPhi_     .clear();  // TODO: not available in reco::
  hyphoE1x3_              .clear();
  hyphoE2x2_              .clear();
  hyphoE3x3_              .clear();
  hyphoE2x5Max_           .clear();
  hyphoE1x5_              .clear();
  hyphoE2x5_              .clear();
  hyphoE5x5_              .clear();
  hyphoMaxEnergyXtal_     .clear();
  hyphoSigmaEtaEta_       .clear();
  hyphoR1x5_              .clear();
  hyphoR2x5_              .clear();
  hyphoESEffSigmaRR_      .clear();
  hyphoSigmaIEtaIEta_2012_.clear();
  hyphoSigmaIEtaIPhi_2012_.clear();
  hyphoSigmaIPhiIPhi_2012_.clear();
  hyphoE1x3_2012_         .clear();
  hyphoE2x2_2012_         .clear();
  hyphoE3x3_2012_         .clear();
  hyphoE2x5Max_2012_      .clear();
  hyphoE5x5_2012_         .clear();
  hyphoBC1E_              .clear();
  hyphoBC1Eta_            .clear();
  hyphoBC2E_              .clear();
  hyphoBC2Eta_            .clear();
  hypho_ecalClusterIsoR2_.clear();
  hypho_ecalClusterIsoR3_.clear();
  hypho_ecalClusterIsoR4_.clear();
  hypho_ecalClusterIsoR5_.clear();
  hypho_hcalRechitIsoR1_.clear();
  hypho_hcalRechitIsoR2_.clear();
  hypho_hcalRechitIsoR3_.clear();
  hypho_hcalRechitIsoR4_.clear();
  hypho_hcalRechitIsoR5_.clear();
  hypho_trackIsoR1PtCut20_.clear();
  hypho_trackIsoR2PtCut20_.clear();
  hypho_trackIsoR3PtCut20_.clear();
  hypho_trackIsoR4PtCut20_.clear();
  hypho_trackIsoR5PtCut20_.clear();
  hypho_swissCrx_.clear();
  hypho_seedTime_.clear();
  hypho_genMatchedIndex_.clear();


  elegsfD0_             .clear();
  elegsfDz_             .clear();
  elegsfD0Err_          .clear();
  elegsfDzErr_          .clear();
  elegsfTrkPt_          .clear();
  elegsfTrkP_           .clear();
  elegsfTrkEta_         .clear();
  elegsfTrkPhi_         .clear();
  elegsfTrkCharge_      .clear();
  elegsfTrkChi2_        .clear();
  elegsfTrkNdof_        .clear();
  elegsfTrkNormalizedChi2_ .clear();
  elegsfTrkValidHits_   .clear();
  elegsfTrkMissHits_    .clear();
  elegsfTrkLayers_      .clear();

  gentrkPt_             .clear();
  gentrkP_              .clear();
  gentrkEta_            .clear();
  gentrkPhi_            .clear();
  gentrkcharge_         .clear();
  gentrkvx_	        .clear();
  gentrkvy_		.clear(); 
  gentrkvz_		.clear();
  gentrknormchi2_       .clear();                 
  gentrkchi2_           .clear();
  gentrkd0_             .clear(); 
  gentrkdxy_            .clear();
  gentrkdz_             .clear();
  gentrkdxyError_       .clear();     
  gentrkdzError_        .clear();
  gentrkValidHits_      .clear();                     
  gentrkMissHits_       .clear();  


  nsc_hybrid_ = 0;
  sc_hybrid_E_          .clear();
  sc_hybrid_Et_         .clear();
  sc_hybrid_Eta_        .clear();
  sc_hybrid_Phi_        .clear();   
  sc_hybrid_x_          .clear();      
  sc_hybrid_y_          .clear();      
  sc_hybrid_z_          .clear();      
  sc_hybrid_EtaWidth_   .clear();         
  sc_hybrid_PhiWidth_   .clear();         
  sc_hybrid_RawE_       .clear();         
  sc_hybrid_RawEt_      .clear();         

  nsc_mult55_ = 0;
  sc_mult55_E_          .clear();
  sc_mult55_Et_         .clear();
  sc_mult55_Eta_        .clear();
  sc_mult55_Phi_        .clear();
  sc_mult55_x_          .clear();         
  sc_mult55_y_          .clear();         
  sc_mult55_z_          .clear();         
  sc_mult55_EtaWidth_   .clear();         
  sc_mult55_PhiWidth_   .clear();         
  sc_mult55_RawE_       .clear();         
  sc_mult55_RawEt_      .clear();         



  
  
  run_    = e.id().run();
  event_  = e.id().event();
  lumis_  = e.luminosityBlock();
  isData_ = e.isRealData();

  // MC truth
  if (doGenParticles_ && !isData_) {
    fillGenPileupInfo(e);
    fillGenParticles(e);
  }
  
  edm::Handle<vector<reco::Vertex> > vtxHandle;
  e.getByToken(vtxCollection_, vtxHandle);
  
  // best-known primary vertex coordinates
  math::XYZPoint pv(0, 0, 0);
  for (vector<reco::Vertex>::const_iterator v = vtxHandle->begin(); v != vtxHandle->end(); ++v){
    //if (!v->isFake()) {
      pv.SetXYZ(v->x(), v->y(), v->z());
      isFakeVtx_    .push_back(v->isFake());
      xVtx_    .push_back(v->x());
      yVtx_    .push_back(v->y());
      zVtx_    .push_back(v->z());
      nVtx_++;
       break;
    }
  
  
  fillElectrons(e, es, pv);
  fillPhotons(e, es, pv);
  fillMuons(e, es, pv);
  fillCaloTower(e, es, pv);
  fillGsfPhotons(e, es, pv);
  fillElectronTracks(e, es, pv);
  fillGeneralTracks(e, es, pv);
  fillSuperClusters(e, es, pv);

  
  tree_->Fill();
}

void ggHiNtuplizer::fillGenPileupInfo(const edm::Event& e)
{
  // Fills information about pileup from MC truth.

  edm::Handle<vector<PileupSummaryInfo> > genPileupHandle;
  e.getByToken(genPileupCollection_, genPileupHandle);

  for (vector<PileupSummaryInfo>::const_iterator pu = genPileupHandle->begin(); pu != genPileupHandle->end(); ++pu) {
    nPU_   .push_back(pu->getPU_NumInteractions());
    puTrue_.push_back(pu->getTrueNumInteractions());
    puBX_  .push_back(pu->getBunchCrossing());

    nPUInfo_++;
  }

}

void ggHiNtuplizer::fillGenParticles(const edm::Event& e)
{
  // Fills tree branches with generated particle info.

  edm::Handle<vector<reco::GenParticle> > genParticlesHandle;
  e.getByToken(genParticlesCollection_, genParticlesHandle);

  int genIndex = 0;

  // loop over MC particles
  for (vector<reco::GenParticle>::const_iterator p = genParticlesHandle->begin(); p != genParticlesHandle->end(); ++p) {
    genIndex++;

    // skip all primary particles if not particle gun MC
    if (!runOnParticleGun_ && !p->mother()) continue;

    // stable particles with pT > 5 GeV
    bool isStableFast = (p->status() == 1 && p->pt() > 0.0);

    // stable leptons
    bool isStableLepton = (p->status() == 1 && abs(p->pdgId()) >= 11 && abs(p->pdgId()) <= 16);

    // (unstable) Z, W, H, top, bottom
    bool isHeavy = (p->pdgId() == 23 || abs(p->pdgId()) == 24 || p->pdgId() == 25 ||
		    abs(p->pdgId()) == 6 || abs(p->pdgId()) == 5);

    // reduce size of output root file
    if (!isStableFast && !isStableLepton && !isHeavy)
      continue;

    mcPID_   .push_back(p->pdgId());
    mcStatus_.push_back(p->status());
    mcVtx_x_ .push_back(p->vx());
    mcVtx_y_ .push_back(p->vy());
    mcVtx_z_ .push_back(p->vz());
    mcPt_    .push_back(p->pt());
    mcP_     .push_back(p->p());
    mcEta_   .push_back(p->eta());
    mcPhi_   .push_back(p->phi());
    mcE_     .push_back(p->energy());
    mcEt_    .push_back(p->et());
    mcMass_  .push_back(p->mass());
   cout <<" gen P:" << p->p() << endl; 
    reco::GenParticleRef partRef = reco::GenParticleRef(
							genParticlesHandle, p - genParticlesHandle->begin());
    genpartparentage::GenParticleParentage particleHistory(partRef);
    
    mcParentage_.push_back(particleHistory.hasLeptonParent()*16   +
			   particleHistory.hasBosonParent()*8     +
			   particleHistory.hasNonPromptParent()*4 +
			   particleHistory.hasQCDParent()*2       +
			   particleHistory.hasExoticParent());
    
    int   momPID  = -999;
    float momPt   = -999;
    float momEta  = -999;
    float momPhi  = -999;
    float momMass = -999;
    int   gmomPID = -999;
    
    if (particleHistory.hasRealParent()) {
      reco::GenParticleRef momRef = particleHistory.parent();
      
      // mother
      if (momRef.isNonnull() && momRef.isAvailable()) {
	momPID  = momRef->pdgId();
	momPt   = momRef->pt();
	momEta  = momRef->eta();
	momPhi  = momRef->phi();
	momMass = momRef->mass();

	// granny
	genpartparentage::GenParticleParentage motherParticle(momRef);
	if (motherParticle.hasRealParent()) {
	  reco::GenParticleRef granny = motherParticle.parent();
	  gmomPID = granny->pdgId();
	}
      }
    }

    mcMomPID_ .push_back(momPID);
    mcMomPt_  .push_back(momPt);
    mcMomEta_ .push_back(momEta);
    mcMomPhi_ .push_back(momPhi);
    mcMomMass_.push_back(momMass);
    mcGMomPID_.push_back(gmomPID);

    mcIndex_  .push_back(genIndex - 1);

    mcCalIsoDR03_.push_back(getGenCalIso(genParticlesHandle, p, 0.3, false, false));
    mcCalIsoDR04_.push_back(getGenCalIso(genParticlesHandle, p, 0.4, false, false));
    mcTrkIsoDR03_.push_back(getGenTrkIso(genParticlesHandle, p, 0.3) );
    mcTrkIsoDR04_.push_back(getGenTrkIso(genParticlesHandle, p, 0.4) );

    nMC_++;

  } // gen-level particles loop

}

float ggHiNtuplizer::getGenCalIso(edm::Handle<vector<reco::GenParticle> > &handle,
				  reco::GenParticleCollection::const_iterator thisPart,
				  float dRMax, bool removeMu, bool removeNu)
{
  // Returns Et sum.

  float etSum = 0;

  for (reco::GenParticleCollection::const_iterator p = handle->begin(); p != handle->end(); ++p) {
    if (p == thisPart) continue;
    if (p->status() != 1) continue;

    // has to come from the same collision
    if (thisPart->collisionId() != p->collisionId())
      continue;

    int pdgCode = abs(p->pdgId());

    // skip muons/neutrinos, if requested
    if (removeMu && pdgCode == 13) continue;
    if (removeNu && (pdgCode == 12 || pdgCode == 14 || pdgCode == 16)) continue;

    // must be within deltaR cone
    float dR = reco::deltaR(thisPart->momentum(), p->momentum());
    if (dR > dRMax) continue;

    etSum += p->et();
  }

  return etSum;
}

float ggHiNtuplizer::getGenTrkIso(edm::Handle<vector<reco::GenParticle> > &handle,
				  reco::GenParticleCollection::const_iterator thisPart, float dRMax)
{
  // Returns pT sum without counting neutral particles.

  float ptSum = 0;

  for (reco::GenParticleCollection::const_iterator p = handle->begin(); p != handle->end(); ++p) {
    if (p == thisPart) continue;
    if (p->status() != 1) continue;
    if (p->charge() == 0) continue;  // do not count neutral particles

    // has to come from the same collision
    if (thisPart->collisionId() != p->collisionId())
      continue;

    // must be within deltaR cone
    float dR = reco::deltaR(thisPart->momentum(), p->momentum());
    if (dR > dRMax) continue;

    ptSum += p->pt();
  }

  return ptSum;
}

void ggHiNtuplizer::fillElectrons(const edm::Event& e, const edm::EventSetup& es, math::XYZPoint& pv)
{
  // Fills tree branches with reco GSF electrons.

  edm::Handle<edm::View<reco::GsfElectron> > gsfElectronsHandle;
  e.getByToken(gsfElectronsCollection_, gsfElectronsHandle);

  edm::Handle<reco::ConversionCollection> conversions;
  edm::Handle<reco::BeamSpot> theBeamSpot;
  edm::Handle<double> rhoH;
  edm::Handle<edm::ValueMap<bool> > veto_id_decisions; 
  edm::Handle<edm::ValueMap<bool> > loose_id_decisions;
  edm::Handle<edm::ValueMap<bool> > medium_id_decisions;
  edm::Handle<edm::ValueMap<bool> > tight_id_decisions; 
  if(doVID_){
    // Get the conversions collection
    e.getByToken(conversionsToken_, conversions);
    
    // Get the beam spot
    e.getByToken(beamSpotToken_,theBeamSpot);
    
    // Get rho value
    e.getByToken(rhoToken_,rhoH);

    // Get the electron ID data from the event stream.
    // Note: this implies that the VID ID modules have been run upstream.
    // If you need more info, check with the EGM group.
    e.getByToken(eleVetoIdMapToken_ , veto_id_decisions);
    e.getByToken(eleLooseIdMapToken_ ,loose_id_decisions);
    e.getByToken(eleMediumIdMapToken_,medium_id_decisions);
    e.getByToken(eleTightIdMapToken_,tight_id_decisions);
  }
  
  
  // loop over electrons
  for (edm::View<reco::GsfElectron>::const_iterator ele = gsfElectronsHandle->begin(); ele != gsfElectronsHandle->end(); ++ele) {
    eleCharge_           .push_back(ele->charge());
    eleChargeConsistent_ .push_back((int)ele->isGsfCtfScPixChargeConsistent());
    eleSCPixCharge_      .push_back(ele->scPixCharge());
    if(!(ele->closestTrack().isNull())) {
      eleCtfCharge_        .push_back(ele->closestTrack()->charge());
    } else {
      eleCtfCharge_        .push_back(-99.);
    }
    //cout << "event "<< e.id().event() << " electron pt : "<< ele->pt()<< endl;
    eleEn_               .push_back(ele->energy());
    eleD0_               .push_back(ele->gsfTrack()->dxy(pv));
    eleDz_               .push_back(ele->gsfTrack()->dz(pv));
    eleD0Err_            .push_back(ele->gsfTrack()->dxyError());
    eleDzErr_            .push_back(ele->gsfTrack()->dzError());
    eleTrkPt_            .push_back(ele->gsfTrack()->pt());
    eleTrkEta_           .push_back(ele->gsfTrack()->eta());
    eleTrkPhi_           .push_back(ele->gsfTrack()->phi());
    eleTrkCharge_        .push_back(ele->gsfTrack()->charge());
    eleTrkChi2_          .push_back(ele->gsfTrack()->chi2());
    eleTrkNdof_          .push_back(ele->gsfTrack()->ndof());
    eleTrkNormalizedChi2_.push_back(ele->gsfTrack()->normalizedChi2());
    eleTrkValidHits_     .push_back(ele->gsfTrack()->numberOfValidHits());
    eleTrkLayers_        .push_back(ele->gsfTrack()->hitPattern().trackerLayersWithMeasurement());
    eleP_                .push_back(ele->p());
    eleP_atVtx_          .push_back(ele->trackMomentumAtVtx().R());
    elePt_               .push_back(ele->pt());
    eleEta_              .push_back(ele->eta());
    elePhi_              .push_back(ele->phi());
    eleSCx_              .push_back(ele->superCluster()->x());
    eleSCy_              .push_back(ele->superCluster()->y());
    eleSCz_              .push_back(ele->superCluster()->z());
    eleSCEn_             .push_back(ele->superCluster()->energy());
    eleESEn_             .push_back(ele->superCluster()->preshowerEnergy());
    eleSCEta_            .push_back(ele->superCluster()->eta());
    eleSCPhi_            .push_back(ele->superCluster()->phi());
    eleSCRawEn_          .push_back(ele->superCluster()->rawEnergy());
    eleSCEtaWidth_       .push_back(ele->superCluster()->etaWidth());
    eleSCPhiWidth_       .push_back(ele->superCluster()->phiWidth());
    eleHoverE_           .push_back(ele->hcalOverEcal());
    eleHoverEBc_         .push_back(ele->hcalOverEcalBc());
    eleEoverP_           .push_back(ele->eSuperClusterOverP());
    eleEoverPInv_        .push_back(fabs(1./ele->ecalEnergy()-1./ele->trackMomentumAtVtx().R()));
    eleBrem_             .push_back(ele->fbrem());
    eledEtaAtVtx_        .push_back(ele->deltaEtaSuperClusterTrackAtVtx());
    eledPhiAtVtx_        .push_back(ele->deltaPhiSuperClusterTrackAtVtx());
    eledEtaSeedAtVtx_    .push_back(ele->deltaEtaSeedClusterTrackAtVtx());
    eleSigmaIEtaIEta_    .push_back(ele->sigmaIetaIeta());
    eleSigmaIPhiIPhi_    .push_back(ele->sigmaIphiIphi());
    //    eleConvVeto_         .push_back((int)ele->passConversionVeto()); // TODO: not available in reco::
    eleMissHits_         .push_back(ele->gsfTrack()->numberOfLostHits());
    eleTrackIso_         .push_back(ele->dr03TkSumPt());
    eleECalIso_          .push_back(ele->dr03EcalRecHitSumEt());
    eleHCalIso_          .push_back(ele->dr03HcalTowerSumEt());
    eleECalDriven_       .push_back(ele->ecalDrivenSeed());
    //      eleESEffSigmaRR_     .push_back(lazyTool.eseffsirir(*(ele->superCluster())));
    
    // full 5x5
    // vector<float> vCovEle = lazyTool_noZS.localCovariances(*(ele->superCluster()->seed()));
    // eleSigmaIEtaIEta_2012_.push_back(isnan(vCovEle[0]) ? 0. : sqrt(vCovEle[0]));
    eleSigmaIEtaIEta_2012_.push_back(ele->full5x5_sigmaIetaIeta() );
    
    // isolation
    reco::GsfElectron::PflowIsolationVariables pfIso = ele->pfIsolationVariables();
    elePFChIso_          .push_back(pfIso.sumChargedHadronPt);
    elePFPhoIso_         .push_back(pfIso.sumPhotonEt);
    elePFNeuIso_         .push_back(pfIso.sumNeutralHadronEt);
    elePFPUIso_          .push_back(pfIso.sumPUPt);

    //calculation on-fly
    pfIsoCalculator pfIsoCal(e,es, pfCollection_, voronoiBkgPF_, pv);
    if (fabs(ele->superCluster()->eta()) > 1.566) {
      elePFChIso03_          .push_back(pfIsoCal.getPfIso(*ele, 1, 0.3, 0.015, 0.));
      elePFChIso04_          .push_back(pfIsoCal.getPfIso(*ele, 1, 0.4, 0.015, 0.));

      elePFPhoIso03_         .push_back(pfIsoCal.getPfIso(*ele, 4, 0.3, 0.08, 0.));
      elePFPhoIso04_         .push_back(pfIsoCal.getPfIso(*ele, 4, 0.4, 0.08, 0.));
    }
    else {
      elePFChIso03_          .push_back(pfIsoCal.getPfIso(*ele, 1, 0.3, 0.0, 0.));
      elePFChIso04_          .push_back(pfIsoCal.getPfIso(*ele, 1, 0.4, 0.0, 0.));
      elePFPhoIso03_         .push_back(pfIsoCal.getPfIso(*ele, 4, 0.3, 0.0, 0.));
      elePFPhoIso04_         .push_back(pfIsoCal.getPfIso(*ele, 4, 0.4, 0.0, 0.));
    }

    elePFNeuIso03_         .push_back(pfIsoCal.getPfIso(*ele, 5, 0.3, 0., 0.));
    elePFNeuIso04_         .push_back(pfIsoCal.getPfIso(*ele, 5, 0.4, 0., 0.));

    eleR9_               .push_back(ele->r9());
    eleE3x3_             .push_back(ele->r9()*ele->superCluster()->energy());
    eleE5x5_             .push_back(ele->e5x5());
    eleR9Full5x5_        .push_back(ele->full5x5_r9());
    eleE3x3Full5x5_      .push_back(ele->full5x5_r9()*ele->superCluster()->energy());
    eleE5x5Full5x5_      .push_back(ele->full5x5_e5x5());

    //seed                                                                                                                               
    NClusters_            .push_back(ele->superCluster()->clustersSize());
    const int N_ECAL      = ele->superCluster()->clustersEnd() - ele->superCluster()->clustersBegin();
    NEcalClusters_        .push_back(std::max(0, N_ECAL-1));
    eleSeedEn_            .push_back(ele->superCluster()->seed()->energy());
    eleSeedEta_           .push_back(ele->superCluster()->seed()->eta());
    eleSeedPhi_           .push_back(ele->superCluster()->seed()->phi());
    //local coordinates
    edm::Ptr<reco::CaloCluster> theseed = ele->superCluster()->seed();
    EcalClusterLocal ecalLocal;
    if(theseed->hitsAndFractions().at(0).first.subdetId()==EcalBarrel) {
        float cryPhi, cryEta, thetatilt, phitilt;
	int ieta, iphi;
        ecalLocal.localCoordsEB(*(theseed), es, cryEta, cryPhi, ieta, iphi, thetatilt, phitilt);
        eleSeedCryEta_  .push_back(cryEta);
        eleSeedCryPhi_  .push_back(cryPhi);
        eleSeedCryIeta_ .push_back(ieta);
        eleSeedCryIphi_ .push_back(iphi);
    }
    else {
        float cryX, cryY, thetatilt, phitilt;
        int ix, iy;
        ecalLocal.localCoordsEE(*(theseed), es, cryX, cryY, ix, iy, thetatilt, phitilt);
        eleSeedCryEta_  .push_back(0);
        eleSeedCryPhi_  .push_back(0);
        eleSeedCryIeta_ .push_back(0);
        eleSeedCryIphi_ .push_back(0);
    }

    if(doVID_){
      float rho = -999 ;
      if (rhoH.isValid())
      rho = *rhoH;

      float eA = effectiveAreas_.getEffectiveArea(fabs(ele->superCluster()->eta()));
      eleEffAreaTimesRho_.push_back(eA*rho);
      
      bool passConvVeto = !ConversionTools::hasMatchedConversion(*ele, 
								 conversions,
								 theBeamSpot->position());
      elepassConversionVeto_.push_back( (int) passConvVeto );

      // seed
      // eleBC1E_             .push_back(ele->superCluster()->seed()->energy());
      // eleBC1Eta_           .push_back(ele->superCluster()->seed()->eta());

      //parameters of the very first PFCluster
      // reco::CaloCluster_iterator bc = ele->superCluster()->clustersBegin();
      // if (bc != ele->superCluster()->clustersEnd()) {
      //    eleBC2E_  .push_back((*bc)->energy());
      //    eleBC2Eta_.push_back((*bc)->eta());
      // }
      // else {
      //    eleBC2E_  .push_back(-99);
      //    eleBC2Eta_.push_back(-99);
      // } 

      const edm::Ptr<reco::GsfElectron> elePtr(gsfElectronsHandle,ele-gsfElectronsHandle->begin()); //value map is keyed of edm::Ptrs so we need to make one
      bool passVetoID   = false;
      bool passLooseID  = false;
      bool passMediumID = false;
      bool passTightID  = false;
      if(veto_id_decisions.isValid()) passVetoID   = (*veto_id_decisions)[elePtr];
      if(veto_id_decisions.isValid()) passLooseID  = (*loose_id_decisions)[elePtr];
      if(veto_id_decisions.isValid()) passMediumID = (*medium_id_decisions)[elePtr];
      if(veto_id_decisions.isValid()) passTightID  = (*tight_id_decisions)[elePtr];
      eleIDVeto_            .push_back((int)passVetoID);
      eleIDLoose_           .push_back((int)passLooseID);
      eleIDMedium_          .push_back((int)passMediumID);
      eleIDTight_           .push_back((int)passTightID);
    }
    //cout << "event:"<< e.id().event()<< "  Supercluster x:"<< ele->superCluster()->x() << endl;
    //cout << "event:"<< e.id().event()<< "  Supercluster x:"<< ele->superCluster()->position().x() << endl;

    nEle_++;

  } // electrons loop
}

void ggHiNtuplizer::fillPhotons(const edm::Event& e, const edm::EventSetup& es, math::XYZPoint& pv)
{
  // Fills tree branches with photons.
  
  edm::Handle<edm::View<reco::Photon> > recoPhotonsHandle;
  e.getByToken(recoPhotonsCollection_, recoPhotonsHandle);
  edm::Handle<edm::ValueMap<reco::HIPhotonIsolation> > recoPhotonHiIsoHandle;
  edm::ValueMap<reco::HIPhotonIsolation> isoMap;
  
  edm::Handle<reco::GsfElectronCollection> elecs_Handle;
  //edm::Handle<edm::View<reco::GsfElectron> > gsfElectronsHandle;
  e.getByToken(input_electrons_, elecs_Handle);
  
  edm::Handle<reco::ConversionCollection> conversions;
  e.getByToken(conversionsToken_, conversions);
  
  edm::Handle<reco::BeamSpot> theBeamSpot;
  e.getByToken(beamSpotToken_,theBeamSpot);
  
  int localIte = -1;
  for (reco::ConversionCollection::const_iterator conv = conversions->begin(); conv != conversions->end(); ++conv) {
    
    localIte++;
    reco::Vertex vtx = conv->conversionVertex();
    reco::ConversionRef refConv(conversions,localIte);
    if (vtx.isValid()) {
      int iel=-1;
      bool foundAnElec = false;
      for(reco::GsfElectronCollection::const_iterator gsfEle = elecs_Handle->begin(); gsfEle!=elecs_Handle->end(); ++gsfEle) {
	iel++;
	if (ConversionTools::matchesConversion(*gsfEle, *conv)) {
	  foundAnElec = true;
	  break;
	}
      } //electron collection
      if (!foundAnElec) iel=-1;
      int ipho=-1;
      bool foundAPhoton = false;
      
      for (edm::View<reco::Photon>::const_iterator pho = recoPhotonsHandle->begin(); pho != recoPhotonsHandle->end(); ++pho) {
	ipho++;
	reco::SuperCluster theSC = *(pho->superCluster());
	reco::ConversionRef theConv = ConversionTools::matchedConversion(theSC,conversions, theBeamSpot->position());
	if (refConv == theConv) { foundAPhoton=true; break;}
      } // photon collection
      if (!foundAPhoton) ipho=-1;
      if ((ipho==-1)&&(iel==-1)) continue;
      
      conv_eleind_  .push_back(iel);
      conv_phoind_  .push_back(ipho);
      conv_vtxprob_ .push_back(TMath::Prob( vtx.chi2(), vtx.ndof()));
      
      math::XYZVector mom(conv->refittedPairMomentum());
      double dbsx = vtx.x() - theBeamSpot->position().x();
      double dbsy = vtx.y() - theBeamSpot->position().y();
      conv_lxy_    .push_back((mom.x()*dbsx + mom.y()*dbsy)/mom.rho());

      int nHitsMax = 0;
      for (std::vector<uint8_t>::const_iterator it = conv->nHitsBeforeVtx().begin(); it!=conv->nHitsBeforeVtx().end(); ++it) {
	if ((*it) > nHitsMax) nHitsMax = (*it);
	//cout << " no. hits before vertex " << nHitsMax<< endl;
      } 
      conv_nhits_bvtx_   .push_back((int)nHitsMax);
    }// vertex valid
    nConv_++ ;
  } // conversion collection
  
  if(useValMapIso_){
    e.getByToken(recoPhotonsHiIso_, recoPhotonHiIsoHandle);
    isoMap = * recoPhotonHiIsoHandle;
  }
  
  // loop over photons
  for (edm::View<reco::Photon>::const_iterator pho = recoPhotonsHandle->begin(); pho != recoPhotonsHandle->end(); ++pho) {
    phoE_             .push_back(pho->energy());
    phoEt_            .push_back(pho->et());
    phoEta_           .push_back(pho->eta());
    phoPhi_           .push_back(pho->phi());
    phoSCSize_        .push_back(pho->superCluster()->clustersSize());
    phoSCE_           .push_back(pho->superCluster()->energy());
    phoSCEt_          .push_back(pho->superCluster()->energy()/cosh(pho->superCluster()->eta()));
    phoSCRawE_        .push_back(pho->superCluster()->rawEnergy());
    phoSCRawEt_       .push_back(pho->superCluster()->rawEnergy()/cosh(pho->superCluster()->eta()));
    phoESEn_          .push_back(pho->superCluster()->preshowerEnergy());
    phoSCEta_         .push_back(pho->superCluster()->eta());
    phoSCPhi_         .push_back(pho->superCluster()->phi());
    phoSCx_           .push_back(pho->superCluster()->x());
    phoSCy_           .push_back(pho->superCluster()->y());
    phoSCz_           .push_back(pho->superCluster()->z());
    phoSCEtaWidth_    .push_back(pho->superCluster()->etaWidth());
    phoSCPhiWidth_    .push_back(pho->superCluster()->phiWidth());
    phoSCBrem_        .push_back(pho->superCluster()->phiWidth()/pho->superCluster()->etaWidth());
    phohasPixelSeed_  .push_back((int)pho->hasPixelSeed());
     


    bool phopassConvVeto = !ConversionTools::hasMatchedPromptElectron(pho->superCluster(), elecs_Handle, conversions, theBeamSpot->position());
    phopassConversionVeto_.push_back( (int) phopassConvVeto );

    //phoEleVeto_       .push_back((int)pho->passElectronVeto());   // TODO: not available in reco::
    phoR9_            .push_back(pho->r9());
    phoHadTowerOverEm_.push_back(pho->hadTowOverEm());
    phoHoverE_        .push_back(pho->hadronicOverEm());
    phoSigmaIEtaIEta_ .push_back(pho->sigmaIetaIeta());
    //phoSigmaIEtaIPhi_ .push_back(pho->sep());   // TODO: not available in reco::
    //phoSigmaIPhiIPhi_ .push_back(pho->spp());   // TODO: not available in reco::

    // additional shower shape variables
    phoE3x3_   .push_back(pho->e3x3());
    phoE1x5_   .push_back(pho->e1x5());
    phoE2x5_   .push_back(pho->e2x5());
    phoE5x5_   .push_back(pho->e5x5());
    phoMaxEnergyXtal_.push_back(pho->maxEnergyXtal());
    phoSigmaEtaEta_.push_back(pho->sigmaEtaEta());
    phoR1x5_.push_back(pho->r1x5());
    phoR2x5_.push_back(pho->r2x5());

    // phoE1x3_          .push_back(lazyTool.e1x3(      *(pho->superCluster()->seed())));
    // phoE2x2_          .push_back(lazyTool.e2x2(      *(pho->superCluster()->seed())));
    // phoE2x5Max_       .push_back(lazyTool.e2x5Max(   *(pho->superCluster()->seed())));
    // phoESEffSigmaRR_  .push_back(lazyTool.eseffsirir(*(pho->superCluster())));

    // full 5x5
    // vector<float> vCov = lazyTool_noZS.localCovariances(*(pho->superCluster()->seed()));
    // phoSigmaIEtaIEta_2012_ .push_back(isnan(vCov[0]) ? 0. : sqrt(vCov[0]));
    // phoSigmaIEtaIPhi_2012_ .push_back(vCov[1]);
    // phoSigmaIPhiIPhi_2012_ .push_back(isnan(vCov[2]) ? 0. : sqrt(vCov[2]));
    phoSigmaIEtaIEta_2012_.push_back(pho->full5x5_sigmaIetaIeta() );

    phoE3x3_2012_.push_back(pho->full5x5_e3x3());
    // phoE1x3_2012_          .push_back(lazyTool_noZS.e1x3(   *(pho->superCluster()->seed())));
    // phoE2x2_2012_          .push_back(lazyTool_noZS.e2x2(   *(pho->superCluster()->seed())));
    // phoE2x5Max_2012_       .push_back(lazyTool_noZS.e2x5Max(*(pho->superCluster()->seed())));
    // phoE5x5_2012_          .push_back(lazyTool_noZS.e5x5(   *(pho->superCluster()->seed())));

    // seed
    // phoBC1E_     .push_back(pho->superCluster()->seed()->energy());
    // phoBC1Eta_   .push_back(pho->superCluster()->seed()->eta());

    // parameters of the very first PFCluster
    // reco::CaloCluster_iterator bc = pho->superCluster()->clustersBegin();
    // if (bc != pho->superCluster()->clustersEnd()) {
    //    phoBC2E_  .push_back((*bc)->energy());
    //    phoBC2Eta_.push_back((*bc)->eta());
    // }
    // else {
    //    phoBC2E_  .push_back(-99);
    //    phoBC2Eta_.push_back(-99);
    // }

    if(useValMapIso_)
    {
      unsigned int idx = pho - recoPhotonsHandle->begin();
      edm::RefToBase<reco::Photon> photonRef = recoPhotonsHandle->refAt(idx);

      pho_ecalClusterIsoR2_.push_back(isoMap[photonRef].ecalClusterIsoR2());
      pho_ecalClusterIsoR3_.push_back(isoMap[photonRef].ecalClusterIsoR3());
      pho_ecalClusterIsoR4_.push_back(isoMap[photonRef].ecalClusterIsoR4());
      pho_ecalClusterIsoR5_.push_back(isoMap[photonRef].ecalClusterIsoR5());
      pho_hcalRechitIsoR1_.push_back(isoMap[photonRef].hcalRechitIsoR1());
      pho_hcalRechitIsoR2_.push_back(isoMap[photonRef].hcalRechitIsoR2());
      pho_hcalRechitIsoR3_.push_back(isoMap[photonRef].hcalRechitIsoR3());
      pho_hcalRechitIsoR4_.push_back(isoMap[photonRef].hcalRechitIsoR4());
      pho_hcalRechitIsoR5_.push_back(isoMap[photonRef].hcalRechitIsoR5());
      pho_trackIsoR1PtCut20_.push_back(isoMap[photonRef].trackIsoR1PtCut20());
      pho_trackIsoR2PtCut20_.push_back(isoMap[photonRef].trackIsoR2PtCut20());
      pho_trackIsoR3PtCut20_.push_back(isoMap[photonRef].trackIsoR3PtCut20());
      pho_trackIsoR4PtCut20_.push_back(isoMap[photonRef].trackIsoR4PtCut20());
      pho_trackIsoR5PtCut20_.push_back(isoMap[photonRef].trackIsoR5PtCut20());
      pho_swissCrx_.push_back(isoMap[photonRef].swissCrx());
      pho_seedTime_.push_back(isoMap[photonRef].seedTime());
    }

   // Conversion
    phoIsConv_.push_back(pho->hasConversionTracks());
    phoNConv_.push_back(pho->conversions().size());
   
    //cout << " conversion size " << pho->conversions().size() << endl; 
    float phoConvInvMass       = 0;
    float phoConvCotTheta      = 0;
    float phoConvEoverP        = 0;
    float phoConvZofPVfromTrks = 0;
    float phoConvMinDist       = 0;
    float phoConvdPhiAtVtx     = 0;
    float phoConvdPhiAtCalo    = 0;
    float phoConvdEtaAtCalo    = 0;
    float phoConvTrkd0_x       = 0;
    float phoConvTrkd0_y       = 0;
    float phoConvTrkPin_x      = 0;
    float phoConvTrkPin_y      = 0;
    float phoConvTrkPout_x     = 0;
    float phoConvTrkPout_y     = 0;
    float phoConvTrkdz_x       = 0;
    float phoConvTrkdz_y       = 0;
    float phoConvTrkdzErr_x    = 0;
    float phoConvTrkdzErr_y    = 0;
    float phoConvChi2          = 0;
    float phoConvChi2Prob      = 0;
    float phoConvNTrks         = 0;
    float phoConvLikeLihood    = 0;
    float phoConvP4_0 = 0.;
    float phoConvP4_1 = 0.;
    float phoConvP4_2 = 0.;
    float phoConvP4_3 = 0.;
    float phoConvVtx_x = 0.;
    float phoConvVtx_y = 0.;
    float phoConvVtx_z = 0.;
    float phoConvVtxErr_x = 0.;
    float phoConvVtxErr_y = 0.;
    float phoConvVtxErr_z = 0.;
    float phoConvPairMomentum_x = 0.;
    float phoConvPairMomentum_y = 0.;
    float phoConvPairMomentum_z = 0.;
    float phoConvRefittedMomentum_x = 0.;
    float phoConvRefittedMomentum_y = 0.;
    float phoConvRefittedMomentum_z = 0.;
    bool  phoConvValidVtx = 0;
    
    float phoConvCharge1 = 0;
    float phoConvCharge2 = 0;

    if (pho->hasConversionTracks()) {
      //cout << " coming in loop " << endl; 
      reco::ConversionRefVector pconversions = pho->conversions();
      //reco::ConversionRef pconv = pconversions[0];
      //cout << " coming in loop " << endl; 
      
      
      
      
      for (unsigned int j=0; j<pconversions.size(); ++j) {
        reco::ConversionRef pconv = pconversions[j];
	//pconv = pconversions[j];
        phoConvValidVtx = pconv->conversionVertex().isValid();
	reco::Vertex vtx=pconv->conversionVertex();
	

	phoConvP4_0 = pconv->refittedPair4Momentum().px();
	phoConvP4_1 = pconv->refittedPair4Momentum().py();
	phoConvP4_2 = pconv->refittedPair4Momentum().pz();
	phoConvP4_3 = pconv->refittedPair4Momentum().energy();
	  
	phoConvVtx_x = vtx.x();
	phoConvVtx_y = vtx.y();
	phoConvVtx_z = vtx.z();
	phoConvVtxErr_x = vtx.xError();
	phoConvVtxErr_y = vtx.yError();
	phoConvVtxErr_z = vtx.zError();
	phoConvPairMomentum_x = pconv->pairMomentum().x();
	phoConvPairMomentum_y = pconv->pairMomentum().y();
	phoConvPairMomentum_z = pconv->pairMomentum().z();
	phoConvRefittedMomentum_x = pconv->refittedPairMomentum().x();
	phoConvRefittedMomentum_y = pconv->refittedPairMomentum().y();
	phoConvRefittedMomentum_z = pconv->refittedPairMomentum().z();
	
	phoConvChi2       = vtx.chi2();
	phoConvChi2Prob   = ChiSquaredProbability(vtx.chi2(), vtx.ndof());
	phoConvNTrks      = pconv->nTracks();
	phoConvLikeLihood = pconv->MVAout();
	
	const std::vector<edm::RefToBase<reco::Track> > tracks = pconv->tracks();
	for (unsigned int k=0; k<tracks.size(); ++k) {
	  if (k==0) {
	    phoConvTrkdz_x    = tracks[k]->dz();
	    phoConvTrkdzErr_x = tracks[k]->dzError();
	    phoConvCharge1    = tracks[k]->charge();
	  } 
	  if (k==1) {
	    phoConvTrkdz_y    = tracks[k]->dz();
	    phoConvTrkdzErr_y = tracks[k]->dzError();
	    phoConvCharge2    = tracks[k]->charge();
	  }
	}

      if (pconv->tracks().size() > 0) {
	phoConvTrkd0_x   = pconv->tracksSigned_d0()[0];
	phoConvTrkPin_x  = sqrt(pconv->tracksPin()[0].Mag2());
	phoConvTrkPout_x = sqrt(pconv->tracksPout()[0].Mag2());
      }
      
      if (pconv->tracks().size() > 1) {
	phoConvTrkd0_y   = pconv->tracksSigned_d0()[1];
	phoConvTrkPin_y  = sqrt(pconv->tracksPin()[1].Mag2());
	phoConvTrkPout_y = sqrt(pconv->tracksPout()[1].Mag2());
      }
      phoConvInvMass       = pconv->pairInvariantMass();
      phoConvCotTheta      = pconv->pairCotThetaSeparation();

      //float eoverp= -99999.;
    
      //phoConvEoverP = phoE_/sqrt(pconv->refittedPairMomentum().Mag2());
      //phoConvEoverP        = pconv->EoverPrefittedTracks();
      //phoConvZofPVfromTrks = pconv->zOfPrimaryVertexFromTracks();
      //phoConvMinDist       = pconv->distOfMinimumApproach();
      //phoConvdPhiAtVtx     = pconv->dPhiTracksAtVtx();
      //phoConvdPhiAtCalo    = pconv->dPhiTracksAtEcal();
      //phoConvdEtaAtCalo    = pconv->dEtaTracksAtEcal();
      
      }// conversion size
    

      
    }
    // fill the vectors
    phoConvInvMass_      .push_back(phoConvInvMass);
    phoConvCotTheta_     .push_back(phoConvCotTheta);
    phoConvEoverP_       .push_back(phoConvEoverP);
    phoConvZofPVfromTrks_.push_back(phoConvZofPVfromTrks);
    phoConvMinDist_      .push_back(phoConvMinDist);
    phoConvdPhiAtVtx_    .push_back(phoConvdPhiAtVtx);
    phoConvdPhiAtCalo_   .push_back(phoConvdPhiAtCalo);
    phoConvdEtaAtCalo_   .push_back(phoConvdEtaAtCalo);
    phoConvTrkd0_x_      .push_back(phoConvTrkd0_x);
    phoConvTrkd0_y_      .push_back(phoConvTrkd0_y);
    phoConvTrkPin_x_     .push_back(phoConvTrkPin_x);
    phoConvTrkPin_y_     .push_back(phoConvTrkPin_y);
    phoConvTrkPout_x_    .push_back(phoConvTrkPout_x);
    phoConvTrkPout_y_    .push_back(phoConvTrkPout_y);
    phoConvTrkdz_x_      .push_back(phoConvTrkdz_x);
    phoConvTrkdz_y_      .push_back(phoConvTrkdz_y);
    phoConvTrkdzErr_x_   .push_back(phoConvTrkdzErr_x);
    phoConvTrkdzErr_y_   .push_back(phoConvTrkdzErr_y);
    phoConvChi2_         .push_back(phoConvChi2);
    phoConvChi2Prob_     .push_back(phoConvChi2Prob);
    phoConvNTrks_        .push_back(phoConvNTrks);
    phoConvCharge1_      .push_back(phoConvCharge1);
    phoConvCharge2_      .push_back(phoConvCharge2);
    phoConvValidVtx_     .push_back(phoConvValidVtx);
    phoConvLikeLihood_   .push_back(phoConvLikeLihood);
    phoConvP4_0_         .push_back(phoConvP4_0);
    phoConvP4_1_         .push_back(phoConvP4_1);
    phoConvP4_2_         .push_back(phoConvP4_2);
    phoConvP4_3_         .push_back(phoConvP4_3);
    phoConvVtx_x_        .push_back(phoConvVtx_x);
    phoConvVtx_y_        .push_back(phoConvVtx_y);
    phoConvVtx_z_        .push_back(phoConvVtx_z);
    phoConvVtxErr_x_     .push_back(phoConvVtxErr_x);
    phoConvVtxErr_y_     .push_back(phoConvVtxErr_y);
    phoConvVtxErr_z_     .push_back(phoConvVtxErr_z);
    phoConvPairMomentum_x_.push_back(phoConvPairMomentum_x);
    phoConvPairMomentum_y_.push_back(phoConvPairMomentum_y);
    phoConvPairMomentum_z_.push_back(phoConvPairMomentum_z);
    phoConvRefittedMomentum_x_.push_back(phoConvRefittedMomentum_x);
    phoConvRefittedMomentum_y_.push_back(phoConvRefittedMomentum_y);
    phoConvRefittedMomentum_z_.push_back(phoConvRefittedMomentum_z);
    
    
    //////////////////////////////////    // MC matching /////////////////////////////////////////
    if (doGenParticles_ && !isData_) {
      float delta2 = 0.15*0.15;
      
      bool gpTemp(false);
      Float_t currentMaxPt(-1);
      int matchedIndex = -1;
      
      for (unsigned igen = 0; igen < mcEt_.size(); ++igen){
	if (mcStatus_[igen] != 1 || (mcPID_[igen]) != 22 ) continue;
	if(reco::deltaR2(pho->eta(), pho->phi(), mcEta_[igen], mcPhi_[igen])<delta2
	   && mcPt_[igen] > currentMaxPt ) {

	  gpTemp = true;
	  currentMaxPt = mcPt_[igen];
	  matchedIndex = igen;
	}
      }

      // if no matching photon was found try with other particles
      std::vector<int> otherPdgIds_(1,11);
      if( !gpTemp ) {
	currentMaxPt = -1;
	for (unsigned igen = 0; igen < mcEt_.size(); ++igen){
	  if (mcStatus_[igen] != 1 || find(otherPdgIds_.begin(),otherPdgIds_.end(),fabs(mcPID_[igen])) == otherPdgIds_.end() ) continue;
	  if(reco::deltaR2(pho->eta(), pho->phi(), mcEta_[igen], mcPhi_[igen])<delta2
	     && mcPt_[igen] > currentMaxPt ) {

	    gpTemp = true;
	    currentMaxPt = mcPt_[igen];
	    matchedIndex = igen;
	  }

	} // end of loop over gen particles
      } // if not matched to gen photon

      pho_genMatchedIndex_.push_back(matchedIndex);
    } // if it's a MC

    nPho_++;

  } // photons loop
}

void ggHiNtuplizer::fillMuons(const edm::Event& e, const edm::EventSetup& es, math::XYZPoint& pv)
{
  // Fills tree branches with reco muons.
  
  edm::Handle<edm::View<reco::Muon> > recoMuonsHandle;
  e.getByToken(recoMuonsCollection_, recoMuonsHandle);
  
  for (edm::View<reco::Muon>::const_iterator mu = recoMuonsHandle->begin(); mu != recoMuonsHandle->end(); ++mu) {
    // if (mu->pt() < 5) continue;
    if (!(mu->isPFMuon() || mu->isGlobalMuon() || mu->isTrackerMuon())) continue;
    
    muPt_    .push_back(mu->pt());
    muEta_   .push_back(mu->eta());
    muPhi_   .push_back(mu->phi());
    muCharge_.push_back(mu->charge());
    muType_  .push_back(mu->type());
    muIsGood_.push_back((int) muon::isGoodMuon(*mu, muon::selectionTypeFromString("TMOneStationTight")));
    muD0_    .push_back(mu->muonBestTrack()->dxy(pv));
    muDz_    .push_back(mu->muonBestTrack()->dz(pv));
    
    const reco::TrackRef glbMu = mu->globalTrack();
    const reco::TrackRef innMu = mu->innerTrack();
    
    if (glbMu.isNull()) {
      muChi2NDF_ .push_back(-99);
      muMuonHits_.push_back(-99);
    } else {
      muChi2NDF_.push_back(glbMu->normalizedChi2());
      muMuonHits_.push_back(glbMu->hitPattern().numberOfValidMuonHits());
    }
    
    if (innMu.isNull()) {
      muInnerD0_     .push_back(-99);
      muInnerDz_     .push_back(-99);
      muTrkLayers_   .push_back(-99);
      muPixelLayers_ .push_back(-99);
      muPixelHits_   .push_back(-99);
      muTrkQuality_  .push_back(-99);
    } else {
      muInnerD0_     .push_back(innMu->dxy(pv));
      muInnerDz_     .push_back(innMu->dz(pv));
      muTrkLayers_   .push_back(innMu->hitPattern().trackerLayersWithMeasurement());
      muPixelLayers_ .push_back(innMu->hitPattern().pixelLayersWithMeasurement());
      muPixelHits_   .push_back(innMu->hitPattern().numberOfValidPixelHits());
      muTrkQuality_  .push_back(innMu->quality(reco::TrackBase::highPurity));
    }
    
    muStations_ .push_back(mu->numberOfMatchedStations());
    muIsoTrk_   .push_back(mu->isolationR03().sumPt);
    muPFChIso_  .push_back(mu->pfIsolationR04().sumChargedHadronPt);
    muPFPhoIso_ .push_back(mu->pfIsolationR04().sumPhotonEt);
    muPFNeuIso_ .push_back(mu->pfIsolationR04().sumNeutralHadronEt);
    muPFPUIso_  .push_back(mu->pfIsolationR04().sumPUPt);
    
    nMu_++;
  }
} // muons loop


void ggHiNtuplizer::fillCaloTower(const edm::Event& e, const edm::EventSetup& es, math::XYZPoint& pv)
{
  // Fills tree branches with reco muons.
  
  //edm::Handle<edm::SortedCollection<CaloTower,edm::StrictWeakOrdering<CaloTower>>> CaloTowerHandle;
  edm::Handle<edm::SortedCollection<CaloTower>> CaloTowerHandle;
  e.getByToken(CaloTowerCollection_, CaloTowerHandle);
  
  //for (edm::SortedCollection<CaloTower,edm::StrictWeakOrdering<CaloTower>> calo = CaloTowerHandle->begin(); calo != CaloTowerHandle->end(); ++calo) {
  for (edm::SortedCollection<CaloTower>::const_iterator calo = CaloTowerHandle->begin(); calo != CaloTowerHandle->end(); ++calo) {

       CaloTower_emE_  .push_back(calo->emEnergy());
       CaloTower_hadE_ .push_back(calo->hadEnergy());
       CaloTower_e_    .push_back(calo->energy());
       CaloTower_et_   .push_back(calo->et());
       CaloTower_phi_  .push_back(calo->phi());
       CaloTower_eta_  .push_back(calo->eta());
    
    nTower_++;
  }
} // calo tower loop

void ggHiNtuplizer::fillGsfPhotons(const edm::Event& e, const edm::EventSetup& es, math::XYZPoint& pv)
{
  // Fills tree branches with photons.
  
  edm::Handle<edm::View<reco::Photon> > recoHyPhotonsHandle;
  e.getByToken(recoHyPhotonsCollection_, recoHyPhotonsHandle);

  edm::Handle<edm::ValueMap<reco::HIPhotonIsolation> > recoHyPhotonHiIsoHandle;
  edm::ValueMap<reco::HIPhotonIsolation> isoMap2;
  
  edm::Handle<reco::BeamSpot> theBeamSpot;
  e.getByToken(beamSpotToken_,theBeamSpot);
  
  
  if(useValMapIso_){
    e.getByToken(recoHyPhotonsHiIso_, recoHyPhotonHiIsoHandle);
    isoMap2 = * recoHyPhotonHiIsoHandle;
  }
  
  // loop over photons
  for (edm::View<reco::Photon>::const_iterator hypho = recoHyPhotonsHandle->begin(); hypho != recoHyPhotonsHandle->end(); ++hypho) {
    hyphoE_             .push_back(hypho->energy());
    hyphoEt_            .push_back(hypho->et());
    hyphoEta_           .push_back(hypho->eta());
    hyphoPhi_           .push_back(hypho->phi());
    hyphoSCSize_        .push_back(hypho->superCluster()->clustersSize());
    hyphoSCE_           .push_back(hypho->superCluster()->energy());
    hyphoSCEt_          .push_back(hypho->superCluster()->energy()/cosh(hypho->superCluster()->eta()));
    hyphoSCRawE_        .push_back(hypho->superCluster()->rawEnergy());
    hyphoSCRawEt_       .push_back(hypho->superCluster()->rawEnergy()/cosh(hypho->superCluster()->eta()));
    hyphoESEn_          .push_back(hypho->superCluster()->preshowerEnergy());
    hyphoSCEta_         .push_back(hypho->superCluster()->eta());
    hyphoSCPhi_         .push_back(hypho->superCluster()->phi());
    hyphoSCx_           .push_back(hypho->superCluster()->x());
    hyphoSCy_           .push_back(hypho->superCluster()->y());
    hyphoSCz_           .push_back(hypho->superCluster()->z());
    hyphoSCEtaWidth_    .push_back(hypho->superCluster()->etaWidth());
    hyphoSCPhiWidth_    .push_back(hypho->superCluster()->phiWidth());
    hyphoSCBrem_        .push_back(hypho->superCluster()->phiWidth()/hypho->superCluster()->etaWidth());
    hyphohasPixelSeed_  .push_back((int)hypho->hasPixelSeed());
     
    //hyphoEleVeto_       .push_back((int)hypho->passElectronVeto());   // TODO: not available in reco::
    hyphoR9_            .push_back(hypho->r9());
    hyphoHadTowerOverEm_.push_back(hypho->hadTowOverEm());
    hyphoHoverE_        .push_back(hypho->hadronicOverEm());
    hyphoSigmaIEtaIEta_ .push_back(hypho->sigmaIetaIeta());
    
    // additional shower shape variables
    hyphoE3x3_   .push_back(hypho->e3x3());
    hyphoE1x5_   .push_back(hypho->e1x5());
    hyphoE2x5_   .push_back(hypho->e2x5());
    hyphoE5x5_   .push_back(hypho->e5x5());
    hyphoMaxEnergyXtal_.push_back(hypho->maxEnergyXtal());
    hyphoSigmaEtaEta_.push_back(hypho->sigmaEtaEta());
    hyphoR1x5_.push_back(hypho->r1x5());
    hyphoR2x5_.push_back(hypho->r2x5());

    hyphoSigmaIEtaIEta_2012_.push_back(hypho->full5x5_sigmaIetaIeta() );
    hyphoE3x3_2012_.push_back(hypho->full5x5_e3x3());

    if(useValMapIso_)
    {
      unsigned int idx = hypho - recoHyPhotonsHandle->begin();
      edm::RefToBase<reco::Photon> photonRef = recoHyPhotonsHandle->refAt(idx);

      hypho_ecalClusterIsoR2_.push_back(isoMap2[photonRef].ecalClusterIsoR2());
      hypho_ecalClusterIsoR3_.push_back(isoMap2[photonRef].ecalClusterIsoR3());
      hypho_ecalClusterIsoR4_.push_back(isoMap2[photonRef].ecalClusterIsoR4());
      hypho_ecalClusterIsoR5_.push_back(isoMap2[photonRef].ecalClusterIsoR5());
      hypho_hcalRechitIsoR1_.push_back(isoMap2[photonRef].hcalRechitIsoR1());
      hypho_hcalRechitIsoR2_.push_back(isoMap2[photonRef].hcalRechitIsoR2());
      hypho_hcalRechitIsoR3_.push_back(isoMap2[photonRef].hcalRechitIsoR3());
      hypho_hcalRechitIsoR4_.push_back(isoMap2[photonRef].hcalRechitIsoR4());
      hypho_hcalRechitIsoR5_.push_back(isoMap2[photonRef].hcalRechitIsoR5());
      hypho_trackIsoR1PtCut20_.push_back(isoMap2[photonRef].trackIsoR1PtCut20());
      hypho_trackIsoR2PtCut20_.push_back(isoMap2[photonRef].trackIsoR2PtCut20());
      hypho_trackIsoR3PtCut20_.push_back(isoMap2[photonRef].trackIsoR3PtCut20());
      hypho_trackIsoR4PtCut20_.push_back(isoMap2[photonRef].trackIsoR4PtCut20());
      hypho_trackIsoR5PtCut20_.push_back(isoMap2[photonRef].trackIsoR5PtCut20());
      hypho_swissCrx_.push_back(isoMap2[photonRef].swissCrx());
      hypho_seedTime_.push_back(isoMap2[photonRef].seedTime());
    }

       
    
    //////////////////////////////////    // MC matching /////////////////////////////////////////
    if (doGenParticles_ && !isData_) {
      float delta2 = 0.15*0.15;
      
      bool gpTemp(false);
      Float_t currentMaxPt(-1);
      int matchedIndex = -1;
      
      for (unsigned igen = 0; igen < mcEt_.size(); ++igen){
	if (mcStatus_[igen] != 1 || (mcPID_[igen]) != 22 ) continue;
	if(reco::deltaR2(hypho->eta(), hypho->phi(), mcEta_[igen], mcPhi_[igen])<delta2
	   && mcPt_[igen] > currentMaxPt ) {

	  gpTemp = true;
	  currentMaxPt = mcPt_[igen];
	  matchedIndex = igen;
	}
      }

      // if no matching photon was found try with other particles
      std::vector<int> otherPdgIds_(1,11);
      if( !gpTemp ) {
	currentMaxPt = -1;
	for (unsigned igen = 0; igen < mcEt_.size(); ++igen){
	  if (mcStatus_[igen] != 1 || find(otherPdgIds_.begin(),otherPdgIds_.end(),fabs(mcPID_[igen])) == otherPdgIds_.end() ) continue;
	  if(reco::deltaR2(hypho->eta(), hypho->phi(), mcEta_[igen], mcPhi_[igen])<delta2
	     && mcPt_[igen] > currentMaxPt ) {

	    gpTemp = true;
	    currentMaxPt = mcPt_[igen];
	    matchedIndex = igen;
	  }

	} // end of loop over gen particles
      } // if not matched to gen photon

      hypho_genMatchedIndex_.push_back(matchedIndex);
    } // if it's a MC

    nHyPho_++;

  } // photons loop
}

void ggHiNtuplizer::fillElectronTracks(const edm::Event& e, const edm::EventSetup& es, math::XYZPoint& pv)
{
  // Fills tree branches with reco GSF electrons tracks.

  edm::Handle<reco::GsfTrackCollection> gsfTracksHandle;
  e.getByToken(gsfTracks_, gsfTracksHandle);
//      for(reco::GsfElectronCollection::const_iterator gsfEle = elecs_Handle->begin(); gsfEle!=elecs_Handle->end(); ++gsfEle) {

 for (reco::GsfTrackCollection::const_iterator gtrk = gsfTracksHandle->begin(); gtrk != gsfTracksHandle->end(); ++gtrk) {
     //cout << "event " << e.id().event() << " electron track pt : " << gtrk->pt() << "  inner momentum:" << gtrk->innerMomentum() << "  outer mom:" << gtrk->outerMomentum() << endl;   

    elegsfD0_               .push_back(gtrk->dxy(pv));
    elegsfDz_               .push_back(gtrk->dz(pv));
    elegsfD0Err_            .push_back(gtrk->dxyError());
    elegsfDzErr_            .push_back(gtrk->dzError());
    elegsfTrkPt_            .push_back(gtrk->pt());
    elegsfTrkP_             .push_back(gtrk->p());
    elegsfTrkEta_           .push_back(gtrk->eta());
    elegsfTrkPhi_           .push_back(gtrk->phi());
    elegsfTrkCharge_        .push_back(gtrk->charge());
    elegsfTrkChi2_          .push_back(gtrk->chi2());
    elegsfTrkNdof_          .push_back(gtrk->ndof());
    elegsfTrkNormalizedChi2_.push_back(gtrk->normalizedChi2());
    elegsfTrkValidHits_     .push_back(gtrk->numberOfValidHits());
    elegsfTrkMissHits_      .push_back(gtrk->numberOfLostHits());
    elegsfTrkLayers_        .push_back(gtrk->hitPattern().trackerLayersWithMeasurement());
    ngsfEle_++;
  }  



}

void ggHiNtuplizer::fillGeneralTracks(const edm::Event& e, const edm::EventSetup& es, math::XYZPoint& pv)
{
  // Fills tree branches with reco GSF electrons tracks.

  edm::Handle<reco::TrackCollection >  genTracksHandle;
  e.getByToken(genTracks_, genTracksHandle);

//      for(reco::GsfElectronCollection::const_iterator gsfEle = elecs_Handle->begin(); gsfEle!=elecs_Handle->end(); ++gsfEle) {

 for (reco::TrackCollection::const_iterator gentrk = genTracksHandle->begin(); gentrk != genTracksHandle->end(); ++gentrk) {
     //cout << "event " << e.id().event() << " electron track pt : " << gtrk->pt() << endl;   
    // cout << "event " << e.id().event() << " gen track pt : " << gentrk->pt() << "  inner momentum:" << gentrk->innerMomentum() << "  outer mom:" << gentrk->outerMomentum() << endl;   
    gentrkPt_             .push_back(gentrk->pt());            
    gentrkP_              .push_back(gentrk->p());            
    gentrkEta_            .push_back(gentrk->eta());           
    gentrkPhi_            .push_back(gentrk->phi());     
    gentrkcharge_         .push_back(gentrk->charge()); 
    gentrkvx_	          .push_back(gentrk->vx());       
    gentrkvy_		  .push_back(gentrk->vy());             
    gentrkvz_		  .push_back(gentrk->vz());  
    gentrknormchi2_       .push_back(gentrk->normalizedChi2());                   
    gentrkchi2_           .push_back(gentrk->chi2());    
    gentrkd0_             .push_back(gentrk->d0());               
    gentrkdxy_            .push_back(gentrk->dxy());             
    gentrkdz_             .push_back(gentrk->dz());     
    gentrkdxyError_       .push_back(gentrk->dxyError());             
    gentrkdzError_        .push_back(gentrk->dzError());   
    gentrkValidHits_      .push_back(gentrk->numberOfValidHits());                     
    gentrkMissHits_       .push_back(gentrk->numberOfLostHits());  

    
    gentrkPurity_            .push_back(gentrk->reco::TrackBase::qualityByName("highPurity"));              
    ngenTrk_++;
  }  



}

void ggHiNtuplizer::fillSuperClusters(const edm::Event& e, const edm::EventSetup& es, math::XYZPoint& pv)
{ 
  edm::Handle<reco::SuperClusterCollection> scHybridHandle;
  e.getByToken(hybridsc_, scHybridHandle);

 for (reco::SuperClusterCollection::const_iterator schy = scHybridHandle->begin(); schy != scHybridHandle->end(); ++schy) {
      sc_hybrid_E_           .push_back(schy->energy());
      sc_hybrid_Et_          .push_back(schy->energy()/cosh(schy->eta()));
      sc_hybrid_Eta_         .push_back(schy->eta());
      sc_hybrid_Phi_         .push_back(schy->phi());
      sc_hybrid_x_           .push_back(schy->x());
      sc_hybrid_y_           .push_back(schy->y());
      sc_hybrid_z_           .push_back(schy->z());
      sc_hybrid_EtaWidth_    .push_back(schy->etaWidth());
      sc_hybrid_PhiWidth_    .push_back(schy->phiWidth());        
      sc_hybrid_RawE_        .push_back(schy->rawEnergy());         
      sc_hybrid_RawEt_       .push_back(schy->rawEnergy()/cosh(schy->eta()));     
     
    nsc_hybrid_ ++;
   
   }


  edm::Handle<reco::SuperClusterCollection> scMult55Handle;
  e.getByToken(mult55sc_, scMult55Handle);

 for (reco::SuperClusterCollection::const_iterator schy = scMult55Handle->begin(); schy != scMult55Handle->end(); ++schy) {
      sc_mult55_E_           .push_back(schy->energy());
      sc_mult55_Et_          .push_back(schy->energy()/cosh(schy->eta()));
      sc_mult55_Eta_         .push_back(schy->eta());
      sc_mult55_Phi_         .push_back(schy->phi());
      sc_mult55_x_           .push_back(schy->x());
      sc_mult55_y_           .push_back(schy->y());
      sc_mult55_z_           .push_back(schy->z());
      sc_mult55_EtaWidth_    .push_back(schy->etaWidth());
      sc_mult55_PhiWidth_    .push_back(schy->phiWidth());        
      sc_mult55_RawE_        .push_back(schy->rawEnergy());         
      sc_mult55_RawEt_       .push_back(schy->rawEnergy()/cosh(schy->eta()));     
     
    nsc_mult55_ ++;
   
   }

}
     

DEFINE_FWK_MODULE(ggHiNtuplizer);
