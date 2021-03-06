/*
  Based on the jet response analyzer
  Modified by Matt Nguyen, November 2010

*/

#include "HeavyIonsAnalysis/JetAnalysis/interface/HiInclusiveJetAnalyzer.h"

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include <Math/DistFunc.h>
#include "TMath.h"

#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "DataFormats/Common/interface/View.h"

#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "SimDataFormats/HiGenData/interface/GenHIEvent.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMap.h"
#include "L1Trigger/GlobalTrigger/plugins/L1GlobalTrigger.h"

#include "RecoBTag/SecondaryVertex/interface/TrackKinematics.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"

#include "fastjet/contrib/Njettiness.hh"

using namespace std;
using namespace edm;
using namespace reco;

HiInclusiveJetAnalyzer::HiInclusiveJetAnalyzer(const edm::ParameterSet& iConfig) :
  geo(0)
{

  doMatch_ = iConfig.getUntrackedParameter<bool>("matchJets",false);
  jetTagLabel_ = iConfig.getParameter<InputTag>("jetTag");
  jetTag_ = consumes<reco::JetView> (jetTagLabel_);
  jetTagPat_ = consumes<pat::JetCollection> (jetTagLabel_);
  matchTag_ = consumes<reco::JetView> (iConfig.getUntrackedParameter<InputTag>("matchTag"));
  matchTagPat_ = consumes<pat::JetCollection> (iConfig.getUntrackedParameter<InputTag>("matchTag"));
  
  // vtxTag_ = iConfig.getUntrackedParameter<edm::InputTag>("vtxTag",edm::InputTag("hiSelectedVertex"));
  vtxTag_ = consumes<vector<reco::Vertex> >        (iConfig.getUntrackedParameter<edm::InputTag>("vtxTag",edm::InputTag("hiSelectedVertex")));  
  // iConfig.getUntrackedParameter<edm::InputTag>("vtxTag",edm::InputTag("hiSelectedVertex"));
  trackTag_ = consumes<reco::TrackCollection> (iConfig.getParameter<InputTag>("trackTag"));
  useQuality_ = iConfig.getUntrackedParameter<bool>("useQuality",1);
  trackQuality_ = iConfig.getUntrackedParameter<string>("trackQuality","highPurity");

  jetName_ = iConfig.getUntrackedParameter<string>("jetName");
  doGenTaus_ = iConfig.getUntrackedParameter<bool>("doGenTaus",0);
  doSubJets_ = iConfig.getUntrackedParameter<bool>("doSubJets",0);
  doGenSubJets_ = iConfig.getUntrackedParameter<bool>("doGenSubJets",false);
  subjetGenTag_ = consumes<reco::JetView> (iConfig.getUntrackedParameter<InputTag>("subjetGenTag"));

  // useGenTaus = true;
  // if (iConfig.exists("genTau1"))
  //   tokenGenTau1_          = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("genTau1"));
  // else useGenTaus = false;
  // if (iConfig.exists("genTau2"))
  //   tokenGenTau2_          = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("genTau2"));
  // else useGenTaus = false;
  // if (iConfig.exists("genTau3"))
  //   tokenGenTau3_          = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("genTau3"));
  // else useGenTaus = false;
  
  isMC_ = iConfig.getUntrackedParameter<bool>("isMC",false);
  useHepMC_ = iConfig.getUntrackedParameter<bool> ("useHepMC",false);
  fillGenJets_ = iConfig.getUntrackedParameter<bool>("fillGenJets",false);

  doTrigger_ = iConfig.getUntrackedParameter<bool>("doTrigger",false);
  doHiJetID_ = iConfig.getUntrackedParameter<bool>("doHiJetID",false);
  if(doHiJetID_) jetIDweightFile_ = iConfig.getUntrackedParameter<string>("jetIDWeight","weights.xml");
  doStandardJetID_ = iConfig.getUntrackedParameter<bool>("doStandardJetID",false);

  rParam = iConfig.getParameter<double>("rParam");
  hardPtMin_ = iConfig.getUntrackedParameter<double>("hardPtMin",4);
  jetPtMin_ = iConfig.getParameter<double>("jetPtMin");

  if(isMC_){
    genjetTag_ = consumes<vector<reco::GenJet> > (iConfig.getParameter<InputTag>("genjetTag"));
    //genjetTag_ = consumes<edm::View<reco::Jet>>(iConfig.getParameter<InputTag>("genjetTag"));
    if(useHepMC_) eventInfoTag_ = consumes<HepMCProduct> (iConfig.getParameter<InputTag>("eventInfoTag"));
    eventGenInfoTag_ = consumes<GenEventInfoProduct> (iConfig.getParameter<InputTag>("eventInfoTag"));
  }
  verbose_ = iConfig.getUntrackedParameter<bool>("verbose",false);

  useVtx_ = iConfig.getUntrackedParameter<bool>("useVtx",false);
  useJEC_ = iConfig.getUntrackedParameter<bool>("useJEC",true);
  usePat_ = iConfig.getUntrackedParameter<bool>("usePAT",true);

  doLifeTimeTagging_ = iConfig.getUntrackedParameter<bool>("doLifeTimeTagging",false);
  doLifeTimeTaggingExtras_ = iConfig.getUntrackedParameter<bool>("doLifeTimeTaggingExtras",true);
  saveBfragments_  = iConfig.getUntrackedParameter<bool>("saveBfragments",false);
  skipCorrections_  = iConfig.getUntrackedParameter<bool>("skipCorrections",false);

  pfCandidateLabel_ = consumes<reco::PFCandidateCollection> (iConfig.getUntrackedParameter<edm::InputTag>("pfCandidateLabel",edm::InputTag("particleFlowTmp")));

  doTower = iConfig.getUntrackedParameter<bool>("doTower",false);
  if(doTower){
    TowerSrc_ = consumes<CaloTowerCollection>( iConfig.getUntrackedParameter<edm::InputTag>("towersSrc",edm::InputTag("towerMaker")));
  }
  
  // EBSrc_ = iConfig.getUntrackedParameter<edm::InputTag>("EBRecHitSrc",edm::InputTag("ecalRecHit","EcalRecHitsEB"));
  // EESrc_ = iConfig.getUntrackedParameter<edm::InputTag>("EERecHitSrc",edm::InputTag("ecalRecHit","EcalRecHitsEE"));
  // HcalRecHitHFSrc_ = iConfig.getUntrackedParameter<edm::InputTag>("hcalHFRecHitSrc",edm::InputTag("hfreco"));
  // HcalRecHitHBHESrc_ = iConfig.getUntrackedParameter<edm::InputTag>("hcalHBHERecHitSrc",edm::InputTag("hbhereco"));
  doExtraCTagging_ = iConfig.getUntrackedParameter<bool>("doExtraCTagging",false);

  genParticleSrc_ = consumes<reco::GenParticleCollection> (iConfig.getUntrackedParameter<edm::InputTag>("genParticles",edm::InputTag("hiGenParticles")));

  if(doTrigger_){
    L1gtReadout_ = consumes< L1GlobalTriggerReadoutRecord > (iConfig.getParameter<edm::InputTag>("L1gtReadout"));
    hltResName_ = consumes< TriggerResults >(iConfig.getUntrackedParameter<string>("hltTrgResults","TriggerResults::HLT"));


    if (iConfig.exists("hltTrgNames"))
      hltTrgNames_ = iConfig.getUntrackedParameter<vector<string> >("hltTrgNames");

    if (iConfig.exists("hltProcNames"))
      hltProcNames_ = iConfig.getUntrackedParameter<vector<string> >("hltProcNames");
    else {
      hltProcNames_.push_back("FU");
      hltProcNames_.push_back("HLT");
    }
  }
  if(doLifeTimeTagging_){
    bTagJetName_ = iConfig.getUntrackedParameter<string>("bTagJetName");
    ImpactParameterTagInfos_ = consumes<vector<TrackIPTagInfo> > (iConfig.getUntrackedParameter<string>("ImpactParameterTagInfos",(bTagJetName_+"ImpactParameterTagInfos")));
    TrackCountingHighEffBJetTags_ = consumes<JetTagCollection> (iConfig.getUntrackedParameter<string>("TrackCountingHighEffBJetTags",(bTagJetName_+"TrackCountingHighEffBJetTags")));
    TrackCountingHighPurBJetTags_ = consumes<JetTagCollection> (iConfig.getUntrackedParameter<string>("TrackCountingHighPurBJetTags",(bTagJetName_+"TrackCountingHighPurBJetTags")));
    JetProbabilityBJetTags_ = consumes<JetTagCollection> (iConfig.getUntrackedParameter<string>("JetProbabilityBJetTags",(bTagJetName_+"JetProbabilityBJetTags")));
    JetBProbabilityBJetTags_ = consumes<JetTagCollection> (iConfig.getUntrackedParameter<string>("JetBProbabilityBJetTags",(bTagJetName_+"JetBProbabilityBJetTags")));
    SecondaryVertexTagInfos_ = consumes<vector<SecondaryVertexTagInfo> > (iConfig.getUntrackedParameter<string>("SecondaryVertexTagInfos",(bTagJetName_+"SecondaryVertexTagInfos")));
    SecondaryVertexNegativeTagInfos_ = consumes<vector<SecondaryVertexTagInfo> > (iConfig.getUntrackedParameter<string>("SecondaryVertexNegativeTagInfos",(bTagJetName_+"SecondaryVertexNegativeTagInfos")));
    SimpleSecondaryVertexHighEffBJetTags_ = consumes<JetTagCollection> (iConfig.getUntrackedParameter<string>("SimpleSecondaryVertexHighEffBJetTags",(bTagJetName_+"SimpleSecondaryVertexHighEffBJetTags")));
    NegativeSimpleSecondaryVertexHighEffBJetTags_ = consumes<JetTagCollection> (iConfig.getUntrackedParameter<string>("NegativeSimpleSecondaryVertexHighEffBJetTags",(bTagJetName_+"NegativeSimpleSecondaryVertexHighEffBJetTags")));
    SimpleSecondaryVertexHighPurBJetTags_ = consumes<JetTagCollection> (iConfig.getUntrackedParameter<string>("SimpleSecondaryVertexHighPurBJetTags",(bTagJetName_+"SimpleSecondaryVertexHighPurBJetTags")));
    NegativeSimpleSecondaryVertexHighPurBJetTags_ = consumes<JetTagCollection> (iConfig.getUntrackedParameter<string>("NegativeSimpleSecondaryVertexHighPurBJetTags",(bTagJetName_+"NegativeSimpleSecondaryVertexHighPurBJetTags")));
    CombinedSecondaryVertexBJetTags_ = consumes<JetTagCollection> (iConfig.getUntrackedParameter<string>("CombinedSecondaryVertexBJetTags",(bTagJetName_+"CombinedSecondaryVertexBJetTags")));
    NegativeCombinedSecondaryVertexBJetTags_ = consumes<JetTagCollection> (iConfig.getUntrackedParameter<string>("NegativeCombinedSecondaryVertexBJetTags",(bTagJetName_+"NegativeCombinedSecondaryVertexBJetTags")));
    PositiveCombinedSecondaryVertexBJetTags_ = consumes<JetTagCollection> (iConfig.getUntrackedParameter<string>("PositiveCombinedSecondaryVertexBJetTags",(bTagJetName_+"PositiveCombinedSecondaryVertexBJetTags")));
    CombinedSecondaryVertexV2BJetTags_ = consumes<JetTagCollection> (iConfig.getUntrackedParameter<string>("CombinedSecondaryVertexV2BJetTags",(bTagJetName_+"CombinedSecondaryVertexV2BJetTags")));
    NegativeCombinedSecondaryVertexV2BJetTags_ = consumes<JetTagCollection> (iConfig.getUntrackedParameter<string>("NegativeCombinedSecondaryVertexV2BJetTags",(bTagJetName_+"NegativeCombinedSecondaryVertexV2BJetTags")));
    PositiveCombinedSecondaryVertexV2BJetTags_ = consumes<JetTagCollection> (iConfig.getUntrackedParameter<string>("PositiveCombinedSecondaryVertexV2BJetTags",(bTagJetName_+"PositiveCombinedSecondaryVertexV2BJetTags")));
    NegativeSoftPFMuonByPtBJetTags_ = consumes<JetTagCollection> (iConfig.getUntrackedParameter<string>("NegativeSoftPFMuonByPtBJetTags",(bTagJetName_+"NegativeSoftPFMuonByPtBJetTags")));
    PositiveSoftPFMuonByPtBJetTags_ = consumes<JetTagCollection> (iConfig.getUntrackedParameter<string>("PositiveSoftPFMuonByPtBJetTags",(bTagJetName_+"PositiveSoftPFMuonByPtBJetTags")));
  }

  //  cout<<" jet collection : "<<jetTag_<<endl;
  doSubEvent_ = 0;

  if(isMC_){
    //     cout<<" genjet collection : "<<genjetTag_<<endl;
    genPtMin_ = iConfig.getUntrackedParameter<double>("genPtMin",10);
    doSubEvent_ = iConfig.getUntrackedParameter<bool>("doSubEvent",0);
  }

  fastjet::contrib::OnePass_KT_Axes     onepass_kt_axes;
  fastjet::contrib::NormalizedMeasure normalizedMeasure(1.,rParam);
  routine_ = std::auto_ptr<fastjet::contrib::Njettiness> ( new fastjet::contrib::Njettiness( onepass_kt_axes, normalizedMeasure) );

  if(doHiJetID_){
	  string inputArrs[] = { "trackMax/jtpt", "trackHardSum/jtpt", "trackHardN/jtpt", "chargedN/jtpt", "chargedHardSum/jtpt", "chargedHardN/jtpt", "photonN/jtpt", "photonHardSum/jtpt", "photonHardN/jtpt", "neutralN/jtpt", "hcalSum/jtpt", "ecalSum/jtpt", "chargedMax/jtpt", "chargedSum/jtpt", "neutralMax/jtpt", "neutralSum/jtpt", "photonMax/jtpt", "photonSum/jtpt", "eSum/jtpt", "muSum/jtpt" };
	  varAddr = std::unique_ptr<float[]> (new float[sizeof(inputArrs)/sizeof(inputArrs[0])]);

	  reader = std::unique_ptr<TMVA::Reader> (new TMVA::Reader("V:Color:!Silent"));
	  for(unsigned int ivar=0; ivar< (sizeof(inputArrs)/sizeof(inputArrs[0])); ivar++){ reader->AddVariable(inputArrs[ivar].c_str(), &(varAddr.get()[ivar])); }
	  edm::FileInPath fp(jetIDweightFile_.data());
	  std::string transFileName = fp.fullPath(); 
	  reader->BookMVA("BDTG",transFileName.c_str());
   }
}



HiInclusiveJetAnalyzer::~HiInclusiveJetAnalyzer() { }



void
HiInclusiveJetAnalyzer::beginRun(const edm::Run& run,
				 const edm::EventSetup & es) {}

void
HiInclusiveJetAnalyzer::beginJob() {

  //string jetTagName = jetTag_.label()+"_tree";
  string jetTagTitle = jetTagLabel_.label()+" Jet Analysis Tree";
  t = fs1->make<TTree>("t",jetTagTitle.c_str());

  //  TTree* t= new TTree("t","Jet Response Analyzer");
  //t->Branch("run",&jets_.run,"run/I");
  t->Branch("evt",&jets_.evt,"evt/I");
  //t->Branch("lumi",&jets_.lumi,"lumi/I");
  t->Branch("b",&jets_.b,"b/F");
  if (useVtx_) {
    t->Branch("vx",&jets_.vx,"vx/F");
    t->Branch("vy",&jets_.vy,"vy/F");
    t->Branch("vz",&jets_.vz,"vz/F");
  }

  t->Branch("nref",&jets_.nref,"nref/I");
  t->Branch("rawpt",jets_.rawpt,"rawpt[nref]/F");
  if(!skipCorrections_) t->Branch("jtpt",jets_.jtpt,"jtpt[nref]/F");
  t->Branch("jteta",jets_.jteta,"jteta[nref]/F");
  t->Branch("jty",jets_.jty,"jty[nref]/F");
  t->Branch("jtphi",jets_.jtphi,"jtphi[nref]/F");
  t->Branch("jtpu",jets_.jtpu,"jtpu[nref]/F");
  t->Branch("jtm",jets_.jtm,"jtm[nref]/F");
  t->Branch("jtarea",jets_.jtarea,"jtarea[nref]/F");

  t->Branch("jtPfCHF",jets_.jtPfCHF,"jtPfCHF[nref]/F");
  t->Branch("jtPfNHF",jets_.jtPfNHF,"jtPfNHF[nref]/F");
  t->Branch("jtPfCEF",jets_.jtPfCEF,"jtPfCEF[nref]/F");
  t->Branch("jtPfNEF",jets_.jtPfNEF,"jtPfNEF[nref]/F");
  t->Branch("jtPfMUF",jets_.jtPfMUF,"jtPfMUF[nref]/F");

  t->Branch("jtPfCHM",jets_.jtPfCHM,"jtPfCHM[nref]/I");
  t->Branch("jtPfNHM",jets_.jtPfNHM,"jtPfNHM[nref]/I");
  t->Branch("jtPfCEM",jets_.jtPfCEM,"jtPfCEM[nref]/I");
  t->Branch("jtPfNEM",jets_.jtPfNEM,"jtPfNEM[nref]/I");
  t->Branch("jtPfMUM",jets_.jtPfMUM,"jtPfMUM[nref]/I");

  t->Branch("jttau1",jets_.jttau1,"jttau1[nref]/F");
  t->Branch("jttau2",jets_.jttau2,"jttau2[nref]/F");
  t->Branch("jttau3",jets_.jttau3,"jttau3[nref]/F");

  if(doSubJets_) {
    t->Branch("jtSubJetPt",&jets_.jtSubJetPt);
    t->Branch("jtSubJetEta",&jets_.jtSubJetEta);
    t->Branch("jtSubJetPhi",&jets_.jtSubJetPhi);
    t->Branch("jtSubJetM",&jets_.jtSubJetM);
  }
  
  // jet ID information, jet composition
  if(doHiJetID_){
    t->Branch("discr_jetID_cuts", jets_.discr_jetID_cuts,"discr_jetID_cuts[nref]/F");
    t->Branch("discr_jetID_bdt", jets_.discr_jetID_bdt,"discr_jetID_bdt[nref]/F");

    t->Branch("discr_fr01", jets_.discr_fr01,"discr_fr01[nref]/F");

    t->Branch("trackMax", jets_.trackMax,"trackMax[nref]/F");
    t->Branch("trackSum", jets_.trackSum,"trackSum[nref]/F");
    t->Branch("trackN", jets_.trackN,"trackN[nref]/I");
    t->Branch("trackHardSum", jets_.trackHardSum,"trackHardSum[nref]/F");
    t->Branch("trackHardN", jets_.trackHardN,"trackHardN[nref]/I");

    t->Branch("chargedMax", jets_.chargedMax,"chargedMax[nref]/F");
    t->Branch("chargedSum", jets_.chargedSum,"chargedSum[nref]/F");
    t->Branch("chargedN", jets_.chargedN,"chargedN[nref]/I");
    t->Branch("chargedHardSum", jets_.chargedHardSum,"chargedHardSum[nref]/F");
    t->Branch("chargedHardN", jets_.chargedHardN,"chargedHardN[nref]/I");

    t->Branch("photonMax", jets_.photonMax,"photonMax[nref]/F");
    t->Branch("photonSum", jets_.photonSum,"photonSum[nref]/F");
    t->Branch("photonN", jets_.photonN,"photonN[nref]/I");
    t->Branch("photonHardSum", jets_.photonHardSum,"photonHardSum[nref]/F");
    t->Branch("photonHardN", jets_.photonHardN,"photonHardN[nref]/I");

    t->Branch("neutralMax", jets_.neutralMax,"neutralMax[nref]/F");
    t->Branch("neutralSum", jets_.neutralSum,"neutralSum[nref]/F");
    t->Branch("neutralN", jets_.neutralN,"neutralN[nref]/I");

    t->Branch("hcalSum", jets_.hcalSum,"hcalSum[nref]/F");
    t->Branch("ecalSum", jets_.ecalSum,"ecalSum[nref]/F");

    t->Branch("eMax", jets_.eMax,"eMax[nref]/F");
    t->Branch("eSum", jets_.eSum,"eSum[nref]/F");
    t->Branch("eN", jets_.eN,"eN[nref]/I");

    t->Branch("muMax", jets_.muMax,"muMax[nref]/F");
    t->Branch("muSum", jets_.muSum,"muSum[nref]/F");
    t->Branch("muN", jets_.muN,"muN[nref]/I");
  }

  if(doStandardJetID_){
    t->Branch("fHPD",jets_.fHPD,"fHPD[nref]/F");
    t->Branch("fRBX",jets_.fRBX,"fRBX[nref]/F");
    t->Branch("n90",jets_.n90,"n90[nref]/I");
    t->Branch("fSubDet1",jets_.fSubDet1,"fSubDet1[nref]/F");
    t->Branch("fSubDet2",jets_.fSubDet2,"fSubDet2[nref]/F");
    t->Branch("fSubDet3",jets_.fSubDet3,"fSubDet3[nref]/F");
    t->Branch("fSubDet4",jets_.fSubDet4,"fSubDet4[nref]/F");
    t->Branch("restrictedEMF",jets_.restrictedEMF,"restrictedEMF[nref]/F");
    t->Branch("nHCAL",jets_.nHCAL,"nHCAL[nref]/I");
    t->Branch("nECAL",jets_.nECAL,"nECAL[nref]/I");
    t->Branch("apprHPD",jets_.apprHPD,"apprHPD[nref]/F");
    t->Branch("apprRBX",jets_.apprRBX,"apprRBX[nref]/F");

    //  t->Branch("hitsInN90",jets_.n90,"hitsInN90[nref]");
    t->Branch("n2RPC",jets_.n2RPC,"n2RPC[nref]/I");
    t->Branch("n3RPC",jets_.n3RPC,"n3RPC[nref]/I");
    t->Branch("nRPC",jets_.nRPC,"nRPC[nref]/I");

    t->Branch("fEB",jets_.fEB,"fEB[nref]/F");
    t->Branch("fEE",jets_.fEE,"fEE[nref]/F");
    t->Branch("fHB",jets_.fHB,"fHB[nref]/F");
    t->Branch("fHE",jets_.fHE,"fHE[nref]/F");
    t->Branch("fHO",jets_.fHO,"fHO[nref]/F");
    t->Branch("fLong",jets_.fLong,"fLong[nref]/F");
    t->Branch("fShort",jets_.fShort,"fShort[nref]/F");
    t->Branch("fLS",jets_.fLS,"fLS[nref]/F");
    t->Branch("fHFOOT",jets_.fHFOOT,"fHFOOT[nref]/F");
  }

  // Jet ID
  if(doMatch_){
    if(!skipCorrections_) t->Branch("matchedPt", jets_.matchedPt,"matchedPt[nref]/F");
    t->Branch("matchedRawPt", jets_.matchedRawPt,"matchedRawPt[nref]/F");
    t->Branch("matchedPu", jets_.matchedPu,"matchedPu[nref]/F");
    t->Branch("matchedR", jets_.matchedR,"matchedR[nref]/F");
  }

  // b-jet discriminators
  if (doLifeTimeTagging_) {

    t->Branch("discr_ssvHighEff",jets_.discr_ssvHighEff,"discr_ssvHighEff[nref]/F");
    t->Branch("discr_ssvHighPur",jets_.discr_ssvHighPur,"discr_ssvHighPur[nref]/F");

    t->Branch("discr_csvV1",jets_.discr_csvV1,"discr_csvV1[nref]/F");
    t->Branch("discr_csvV2",jets_.discr_csvV2,"discr_csvV2[nref]/F");
    t->Branch("discr_muByIp3",jets_.discr_muByIp3,"discr_muByIp3[nref]/F");
    t->Branch("discr_muByPt",jets_.discr_muByPt,"discr_muByPt[nref]/F");
    t->Branch("discr_prob",jets_.discr_prob,"discr_prob[nref]/F");
    t->Branch("discr_probb",jets_.discr_probb,"discr_probb[nref]/F");
    t->Branch("discr_tcHighEff",jets_.discr_tcHighEff,"discr_tcHighEff[nref]/F");
    t->Branch("discr_tcHighPur",jets_.discr_tcHighPur,"discr_tcHighPur[nref]/F");

    t->Branch("ndiscr_ssvHighEff",jets_.ndiscr_ssvHighEff,"ndiscr_ssvHighEff[nref]/F");
    t->Branch("ndiscr_ssvHighPur",jets_.ndiscr_ssvHighPur,"ndiscr_ssvHighPur[nref]/F");
    t->Branch("ndiscr_csvV1",jets_.ndiscr_csvV1,"ndiscr_csvV1[nref]/F");
    t->Branch("ndiscr_csvV2",jets_.ndiscr_csvV2,"ndiscr_csvV2[nref]/F");
    t->Branch("ndiscr_muByPt",jets_.ndiscr_muByPt,"ndiscr_muByPt[nref]/F");

    t->Branch("pdiscr_csvV1",jets_.pdiscr_csvV1,"pdiscr_csvV1[nref]/F");
    t->Branch("pdiscr_csvV2",jets_.pdiscr_csvV2,"pdiscr_csvV2[nref]/F");

    t->Branch("nsvtx",    jets_.nsvtx,    "nsvtx[nref]/I");
    t->Branch("svtxntrk", jets_.svtxntrk, "svtxntrk[nref]/I");
    t->Branch("svtxdl",   jets_.svtxdl,   "svtxdl[nref]/F");
    t->Branch("svtxdls",  jets_.svtxdls,  "svtxdls[nref]/F");
    t->Branch("svtxdl2d", jets_.svtxdl2d, "svtxdl2d[nref]/F");
    t->Branch("svtxdls2d", jets_.svtxdls2d, "svtxdls2d[nref]/F");
    t->Branch("svtxm",    jets_.svtxm,    "svtxm[nref]/F");
    t->Branch("svtxpt",   jets_.svtxpt,   "svtxpt[nref]/F");
    t->Branch("svtxmcorr", jets_.svtxmcorr, "svtxmcorr[nref]/F");

    t->Branch("nIPtrk",jets_.nIPtrk,"nIPtrk[nref]/I");
    t->Branch("nselIPtrk",jets_.nselIPtrk,"nselIPtrk[nref]/I");

    if (doLifeTimeTaggingExtras_) {
      t->Branch("nIP",&jets_.nIP,"nIP/I");
      t->Branch("ipJetIndex",jets_.ipJetIndex,"ipJetIndex[nIP]/I");
      t->Branch("ipPt",jets_.ipPt,"ipPt[nIP]/F");
      t->Branch("ipProb0",jets_.ipProb0,"ipProb0[nIP]/F");
      t->Branch("ipProb1",jets_.ipProb1,"ipProb1[nIP]/F");
      t->Branch("ip2d",jets_.ip2d,"ip2d[nIP]/F");
      t->Branch("ip2dSig",jets_.ip2dSig,"ip2dSig[nIP]/F");
      t->Branch("ip3d",jets_.ip3d,"ip3d[nIP]/F");
      t->Branch("ip3dSig",jets_.ip3dSig,"ip3dSig[nIP]/F");
      t->Branch("ipDist2Jet",jets_.ipDist2Jet,"ipDist2Jet[nIP]/F");
      t->Branch("ipDist2JetSig",jets_.ipDist2JetSig,"ipDist2JetSig[nIP]/F");
      t->Branch("ipClosest2Jet",jets_.ipClosest2Jet,"ipClosest2Jet[nIP]/F");
      if(doExtraCTagging_){
	      t->Branch("trackPtRel",jets_.trackPtRel,"trackPtRel[nIP]/F");
	      t->Branch("trackPPar",jets_.trackPPar,"trackPPar[nIP]/F");
	      t->Branch("trackPParRatio",jets_.trackPParRatio,"trackPParRatio[nIP]/F");
	      t->Branch("trackDeltaR",jets_.trackDeltaR,"trackDeltaR[nIP]/F");
	      t->Branch("trackPtRatio",jets_.trackPtRatio,"trackPtRatio[nIP]/F");
	      t->Branch("trackSumJetDeltaR",jets_.trackSumJetDeltaR,"trackSumJetDeltaR[nref]/F");
      }

    }

    if (doExtraCTagging_){
	    t->Branch("svtxTrkSumChi2", jets_.svtxTrkSumChi2,"svtxTrkSumChi2[nref]/F");
	    t->Branch("svtxTrkNetCharge",jets_.svtxTrkNetCharge,"svtxTrkNetCharge[nref]/I");
	    t->Branch("svtxNtrkInCone",jets_.svtxNtrkInCone,"svtxNtrkInCone[nref]/I");

	    t->Branch("svJetDeltaR", jets_.svJetDeltaR, "svJetDeltaR[nref]/F");
	    t->Branch("trackSip2dSigAboveCharm",jets_.trackSip2dSigAboveCharm, "trackSip2dSigAboveCharm[nref]/F");
	    t->Branch("trackSip3dSigAboveCharm",jets_.trackSip3dSigAboveCharm, "trackSip3dSigAboveCharm[nref]/F");
	    t->Branch("trackSip2dValAboveCharm",jets_.trackSip2dValAboveCharm, "trackSip2dValAboveCharm[nref]/F");
	    t->Branch("trackSip3dValAboveCharm",jets_.trackSip3dValAboveCharm, "trackSip3dValAboveCharm[nref]/F");

    }

    t->Branch("mue",     jets_.mue,     "mue[nref]/F");
    t->Branch("mupt",    jets_.mupt,    "mupt[nref]/F");
    t->Branch("mueta",   jets_.mueta,   "mueta[nref]/F");
    t->Branch("muphi",   jets_.muphi,   "muphi[nref]/F");
    t->Branch("mudr",    jets_.mudr,    "mudr[nref]/F");
    t->Branch("muptrel", jets_.muptrel, "muptrel[nref]/F");
    t->Branch("muchg",   jets_.muchg,   "muchg[nref]/I");
  }


  if(isMC_){
    t->Branch("beamId1",&jets_.beamId1,"beamId1/I");
    t->Branch("beamId2",&jets_.beamId2,"beamId2/I");

    t->Branch("pthat",&jets_.pthat,"pthat/F");

    // Only matched gen jets
    t->Branch("refpt",jets_.refpt,"refpt[nref]/F");
    t->Branch("refeta",jets_.refeta,"refeta[nref]/F");
    t->Branch("refy",jets_.refy,"refy[nref]/F");
    t->Branch("refphi",jets_.refphi,"refphi[nref]/F");
    t->Branch("refm",jets_.refm,"refm[nref]/F");
    t->Branch("refarea",jets_.refarea,"refarea[nref]/F");
    if(doGenTaus_) {
      t->Branch("reftau1",jets_.reftau1,"reftau1[nref]/F");
      t->Branch("reftau2",jets_.reftau2,"reftau2[nref]/F");
      t->Branch("reftau3",jets_.reftau3,"reftau3[nref]/F");
    }
    t->Branch("refdphijt",jets_.refdphijt,"refdphijt[nref]/F");
    t->Branch("refdrjt",jets_.refdrjt,"refdrjt[nref]/F");
    // matched parton
    t->Branch("refparton_pt",jets_.refparton_pt,"refparton_pt[nref]/F");
    t->Branch("refparton_flavor",jets_.refparton_flavor,"refparton_flavor[nref]/I");
    t->Branch("refparton_flavorForB",jets_.refparton_flavorForB,"refparton_flavorForB[nref]/I");

    if(doGenSubJets_) {
      t->Branch("refptG",jets_.refptG,"refptG[nref]/F");
      t->Branch("refetaG",jets_.refetaG,"refetaG[nref]/F");
      t->Branch("refphiG",jets_.refphiG,"refphiG[nref]/F");
      t->Branch("refmG",jets_.refmG,"refmG[nref]/F");
      t->Branch("refSubJetPt",&jets_.refSubJetPt);
      t->Branch("refSubJetEta",&jets_.refSubJetEta);
      t->Branch("refSubJetPhi",&jets_.refSubJetPhi);
      t->Branch("refSubJetM",&jets_.refSubJetM);
    }
    
    t->Branch("genChargedSum", jets_.genChargedSum,"genChargedSum[nref]/F");
    t->Branch("genHardSum", jets_.genHardSum,"genHardSum[nref]/F");
    t->Branch("signalChargedSum", jets_.signalChargedSum,"signalChargedSum[nref]/F");
    t->Branch("signalHardSum", jets_.signalHardSum,"signalHardSum[nref]/F");

    if(doSubEvent_){
      t->Branch("subid",jets_.subid,"subid[nref]/I");
    }

    if(fillGenJets_){
      // For all gen jets, matched or unmatched
      t->Branch("ngen",&jets_.ngen,"ngen/I");
      t->Branch("genmatchindex",jets_.genmatchindex,"genmatchindex[ngen]/I");
      t->Branch("genpt",jets_.genpt,"genpt[ngen]/F");
      t->Branch("geneta",jets_.geneta,"geneta[ngen]/F");
      t->Branch("geny",jets_.geny,"geny[ngen]/F");
      if(doGenTaus_) {
        t->Branch("gentau1",jets_.gentau1,"gentau1[ngen]/F");
        t->Branch("gentau2",jets_.gentau2,"gentau2[ngen]/F");
        t->Branch("gentau3",jets_.gentau3,"gentau3[ngen]/F");
      }
      t->Branch("genphi",jets_.genphi,"genphi[ngen]/F");
      t->Branch("genm",jets_.genm,"genm[ngen]/F");
      t->Branch("gendphijt",jets_.gendphijt,"gendphijt[ngen]/F");
      t->Branch("gendrjt",jets_.gendrjt,"gendrjt[ngen]/F");

      if(doGenSubJets_) {
        t->Branch("genptG",jets_.genptG,"genptG[ngen]/F");
        t->Branch("genetaG",jets_.genetaG,"genetaG[ngen]/F");
        t->Branch("genphiG",jets_.genphiG,"genphiG[ngen]/F");
        t->Branch("genmG",jets_.genmG,"genmG[ngen]/F");
        t->Branch("genSubJetPt",&jets_.genSubJetPt);
        t->Branch("genSubJetEta",&jets_.genSubJetEta);
        t->Branch("genSubJetPhi",&jets_.genSubJetPhi);
        t->Branch("genSubJetM",&jets_.genSubJetM);
      }

      if(doSubEvent_){
	t->Branch("gensubid",jets_.gensubid,"gensubid[ngen]/I");
      }
    }

    if(saveBfragments_  ) {
      t->Branch("bMult",&jets_.bMult,"bMult/I");
      t->Branch("bJetIndex",jets_.bJetIndex,"bJetIndex[bMult]/I");
      t->Branch("bStatus",jets_.bStatus,"bStatus[bMult]/I");
      t->Branch("bVx",jets_.bVx,"bVx[bMult]/F");
      t->Branch("bVy",jets_.bVy,"bVy[bMult]/F");
      t->Branch("bVz",jets_.bVz,"bVz[bMult]/F");
      t->Branch("bPt",jets_.bPt,"bPt[bMult]/F");
      t->Branch("bEta",jets_.bEta,"bEta[bMult]/F");
      t->Branch("bPhi",jets_.bPhi,"bPhi[bMult]/F");
      t->Branch("bPdg",jets_.bPdg,"bPdg[bMult]/I");
      t->Branch("bChg",jets_.bChg,"bChg[bMult]/I");
    }

  }
  /*
    if(!isMC_){
    t->Branch("nL1TBit",&jets_.nL1TBit,"nL1TBit/I");
    t->Branch("l1TBit",jets_.l1TBit,"l1TBit[nL1TBit]/O");

    t->Branch("nL1ABit",&jets_.nL1ABit,"nL1ABit/I");
    t->Branch("l1ABit",jets_.l1ABit,"l1ABit[nL1ABit]/O");

    t->Branch("nHLTBit",&jets_.nHLTBit,"nHLTBit/I");
    t->Branch("hltBit",jets_.hltBit,"hltBit[nHLTBit]/O");

    }
  */
  TH1D::SetDefaultSumw2();

}


void
HiInclusiveJetAnalyzer::analyze(const Event& iEvent,
				const EventSetup& iSetup) {

  int event = iEvent.id().event();
  int run = iEvent.id().run();
  int lumi = iEvent.id().luminosityBlock();

  jets_.run = run;
  jets_.evt = event;
  jets_.lumi = lumi;

  LogDebug("HiInclusiveJetAnalyzer")<<"START event: "<<event<<" in run "<<run<<endl;

  int bin = -1;
  double hf = 0.;
  double b = 999.;

  if(doHiJetID_ && !geo){
    edm::ESHandle<CaloGeometry> pGeo;
    iSetup.get<CaloGeometryRecord>().get(pGeo);
    geo = pGeo.product();
  }

  // loop the events

  jets_.bin = bin;
  jets_.hf = hf;

  reco::Vertex::Point vtx(0,0,0);
  if (useVtx_) {
    edm::Handle<vector<reco::Vertex> >vertex;
    iEvent.getByToken(vtxTag_, vertex);

    if(vertex->size()>0) {
      jets_.vx=vertex->begin()->x();
      jets_.vy=vertex->begin()->y();
      jets_.vz=vertex->begin()->z();
      vtx = vertex->begin()->position();
    }
  }

  edm::Handle<pat::JetCollection> patjets;
  if(usePat_)iEvent.getByToken(jetTagPat_, patjets);

  edm::Handle<pat::JetCollection> patmatchedjets;
  iEvent.getByToken(matchTagPat_, patmatchedjets);

  edm::Handle<reco::JetView> matchedjets;
  iEvent.getByToken(matchTag_, matchedjets);

  if(doGenSubJets_)
    iEvent.getByToken(subjetGenTag_, gensubjets_);
  
  edm::Handle<reco::JetView> jets;
  iEvent.getByToken(jetTag_, jets);

  edm::Handle<reco::PFCandidateCollection> pfCandidates;
  iEvent.getByToken(pfCandidateLabel_,pfCandidates);

  edm::Handle<reco::TrackCollection> tracks;
  iEvent.getByToken(trackTag_,tracks);

  // edm::Handle<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > > ebHits;
  // edm::Handle<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > > eeHits;

  //edm::Handle<HFRecHitCollection> hfHits;
  //edm::Handle<HBHERecHitCollection> hbheHits;


  edm::Handle<reco::GenParticleCollection> genparts;
  iEvent.getByToken(genParticleSrc_,genparts);

  //Get all the b-tagging handles
  // Track Counting Taggers
  //------------------------------------------------------
  Handle<JetTagCollection> jetTags_TCHighEff;
  Handle<JetTagCollection> jetTags_TCHighPur;

  //------------------------------------------------------
  // Jet Probability tagger
  //------------------------------------------------------
  Handle<vector<TrackIPTagInfo> > tagInfo;
  Handle<JetTagCollection> jetTags_JP;
  Handle<JetTagCollection> jetTags_JB;

  //------------------------------------------------------
  // Secondary vertex taggers
  //------------------------------------------------------
  Handle<vector<SecondaryVertexTagInfo> > tagInfoSVx;
  Handle<vector<SecondaryVertexTagInfo> > tagInfoNegSVx;
  Handle<JetTagCollection> jetTags_SvtxHighEff;
  Handle<JetTagCollection> jetTags_negSvtxHighEff;
  Handle<JetTagCollection> jetTags_SvtxHighPur;
  Handle<JetTagCollection> jetTags_negSvtxHighPur;
  Handle<JetTagCollection> jetTags_CombinedSvtx;
  Handle<JetTagCollection> jetTags_negCombinedSvtx;
  Handle<JetTagCollection> jetTags_posCombinedSvtx;
  Handle<JetTagCollection> jetTags_CombinedSvtxV2;
  Handle<JetTagCollection> jetTags_negCombinedSvtxV2;
  Handle<JetTagCollection> jetTags_posCombinedSvtxV2;

  //------------------------------------------------------
  // Soft muon tagger
  //------------------------------------------------------

  //Handle<SoftLeptonTagInfoCollection> tagInos_softmuon;
  //Handle<JetTagCollection> jetTags_softMu;
  //Handle<JetTagCollection> jetTags_softMuneg;

  if(doLifeTimeTagging_){
    iEvent.getByToken(ImpactParameterTagInfos_, tagInfo);
    iEvent.getByToken(TrackCountingHighEffBJetTags_, jetTags_TCHighEff);
    iEvent.getByToken(TrackCountingHighPurBJetTags_, jetTags_TCHighPur);
    iEvent.getByToken(JetProbabilityBJetTags_, jetTags_JP);
    iEvent.getByToken(JetBProbabilityBJetTags_, jetTags_JB);
    iEvent.getByToken(SecondaryVertexTagInfos_, tagInfoSVx);
    iEvent.getByToken(SecondaryVertexNegativeTagInfos_, tagInfoNegSVx);
    iEvent.getByToken(SimpleSecondaryVertexHighEffBJetTags_, jetTags_SvtxHighEff);
    iEvent.getByToken(NegativeSimpleSecondaryVertexHighEffBJetTags_, jetTags_negSvtxHighEff);
    iEvent.getByToken(SimpleSecondaryVertexHighPurBJetTags_, jetTags_SvtxHighPur);
    iEvent.getByToken(NegativeSimpleSecondaryVertexHighEffBJetTags_, jetTags_negSvtxHighPur);
    iEvent.getByToken(CombinedSecondaryVertexBJetTags_, jetTags_CombinedSvtx);
    iEvent.getByToken(NegativeCombinedSecondaryVertexBJetTags_, jetTags_negCombinedSvtx);
    iEvent.getByToken(PositiveCombinedSecondaryVertexBJetTags_, jetTags_posCombinedSvtx);
    iEvent.getByToken(CombinedSecondaryVertexV2BJetTags_, jetTags_CombinedSvtxV2);
    iEvent.getByToken(NegativeCombinedSecondaryVertexV2BJetTags_, jetTags_negCombinedSvtxV2);
    iEvent.getByToken(PositiveCombinedSecondaryVertexV2BJetTags_, jetTags_posCombinedSvtxV2);
    //iEvent.getByToken(NegativeSoftPFMuonByPtBJetTags_, jetTags_softMuneg);
    //iEvent.getByToken(PositiveSoftPFMuonByPtBJetTags_, jetTags_softMu);
  }


  // get tower information
  edm::Handle<CaloTowerCollection> towers;
  if(doTower){
    iEvent.getByToken(TowerSrc_,towers);
  }

  // FILL JRA TREE
  jets_.b = b;
  jets_.nref = 0;

  if(doTrigger_){
    fillL1Bits(iEvent);
    fillHLTBits(iEvent);
  }

  for(unsigned int j = 0; j < jets->size(); ++j){
    const reco::Jet& jet = (*jets)[j];

    if(jet.pt() < jetPtMin_) continue;
    if (useJEC_ && usePat_){
      jets_.rawpt[jets_.nref]=(*patjets)[j].correctedJet("Uncorrected").pt();
    }

    math::XYZVector jetDir = jet.momentum().Unit();

    if(doLifeTimeTagging_){
      int ith_tagged =    this->TaggedJet(jet,jetTags_SvtxHighEff);
      if(ith_tagged >= 0){
	jets_.discr_ssvHighEff[jets_.nref] = (*jetTags_SvtxHighEff)[ith_tagged].second;
	const SecondaryVertexTagInfo &tagInfoSV = (*tagInfoSVx)[ith_tagged];
	jets_.nsvtx[jets_.nref]     = tagInfoSV.nVertices();

	if ( jets_.nsvtx[jets_.nref] > 0) {

	  jets_.svtxntrk[jets_.nref]  = tagInfoSV.nVertexTracks(0);
	  if(doExtraCTagging_) jets_.svJetDeltaR[jets_.nref] = reco::deltaR(tagInfoSV.flightDirection(0),jetDir);

	  // this is the 3d flight distance, for 2-D use (0,true)
	  Measurement1D m1D = tagInfoSV.flightDistance(0);
	  jets_.svtxdl[jets_.nref]    = m1D.value();
	  jets_.svtxdls[jets_.nref]   = m1D.significance();
	  Measurement1D m2D = tagInfoSV.flightDistance(0,true);
          jets_.svtxdl2d[jets_.nref] = m2D.value();
          jets_.svtxdls2d[jets_.nref] = m2D.significance();
          const Vertex& svtx = tagInfoSV.secondaryVertex(0);
	  
	  double svtxM = svtx.p4().mass();
          double svtxPt = svtx.p4().pt();
	  jets_.svtxm[jets_.nref]    = svtxM; 
	  jets_.svtxpt[jets_.nref]   = svtxPt;
	 
	  if(doExtraCTagging_){
		  double trkSumChi2=0;
		  int trkNetCharge=0;
		  int nTrkInCone=0;
		  int nsvtxtrks=0;
		  TrackRefVector svtxTracks = tagInfoSV.vertexTracks(0);
		  for(unsigned int itrk=0; itrk<svtxTracks.size(); itrk++){
			  nsvtxtrks++;
			  trkSumChi2 += svtxTracks.at(itrk)->chi2();
			  trkNetCharge += svtxTracks.at(itrk)->charge();
			  if(reco::deltaR(svtxTracks.at(itrk)->momentum(),jetDir)<rParam) nTrkInCone++;
		  }
		  jets_.svtxTrkSumChi2[jets_.nref] = trkSumChi2;
		  jets_.svtxTrkNetCharge[jets_.nref] = trkNetCharge;
		  jets_.svtxNtrkInCone[jets_.nref] = nTrkInCone;

		  //try out the corrected mass (http://arxiv.org/pdf/1504.07670v1.pdf)
		  //mCorr=srqt(m^2+p^2sin^2(th)) + p*sin(th)
		  double sinth = svtx.p4().Vect().Unit().Cross(tagInfoSV.flightDirection(0).unit()).Mag2();
		  sinth = sqrt(sinth);
		  jets_.svtxmcorr[jets_.nref] = sqrt(pow(svtxM,2)+(pow(svtxPt,2)*pow(sinth,2)))+svtxPt*sinth;
	  }	

	  if(svtx.ndof()>0)jets_.svtxnormchi2[jets_.nref]  = svtx.chi2()/svtx.ndof();
	}
      }
      ith_tagged    = this->TaggedJet(jet,jetTags_negSvtxHighEff);
      if(ith_tagged >= 0) jets_.ndiscr_ssvHighEff[jets_.nref]   = (*jetTags_negSvtxHighEff)[ith_tagged].second;
      ith_tagged    = this->TaggedJet(jet,jetTags_SvtxHighPur);
      if(ith_tagged >= 0) jets_.discr_ssvHighPur[jets_.nref]  = (*jetTags_SvtxHighPur)[ith_tagged].second;
      ith_tagged    = this->TaggedJet(jet,jetTags_negSvtxHighPur);
      if(ith_tagged >= 0) jets_.ndiscr_ssvHighPur[jets_.nref] = (*jetTags_negSvtxHighPur)[ith_tagged].second;

      ith_tagged          = this->TaggedJet(jet,jetTags_CombinedSvtx);
      if(ith_tagged >= 0) jets_.discr_csvV1[jets_.nref]  = (*jetTags_CombinedSvtx)[ith_tagged].second;
      ith_tagged          = this->TaggedJet(jet,jetTags_negCombinedSvtx);
      if(ith_tagged >= 0) jets_.ndiscr_csvV1[jets_.nref] = (*jetTags_negCombinedSvtx)[ith_tagged].second;
      ith_tagged  = this->TaggedJet(jet,jetTags_posCombinedSvtx);
      if(ith_tagged >= 0) jets_.pdiscr_csvV1[jets_.nref] = (*jetTags_posCombinedSvtx)[ith_tagged].second;

      ith_tagged          = this->TaggedJet(jet,jetTags_CombinedSvtxV2);
      if(ith_tagged >= 0) jets_.discr_csvV2[jets_.nref]  = (*jetTags_CombinedSvtxV2)[ith_tagged].second;
      ith_tagged          = this->TaggedJet(jet,jetTags_negCombinedSvtxV2);
      if(ith_tagged >= 0) jets_.ndiscr_csvV2[jets_.nref] = (*jetTags_negCombinedSvtxV2)[ith_tagged].second;
      ith_tagged  = this->TaggedJet(jet,jetTags_posCombinedSvtxV2);
      if(ith_tagged >= 0) jets_.pdiscr_csvV2[jets_.nref] = (*jetTags_posCombinedSvtxV2)[ith_tagged].second;

      if(ith_tagged >= 0){
	ith_tagged = this->TaggedJet(jet,jetTags_JP);
	jets_.discr_prob[jets_.nref]  = (*jetTags_JP)[ith_tagged].second;

	const TrackIPTagInfo& tagInfoIP= (*tagInfo)[ith_tagged];

	jets_.nIPtrk[jets_.nref] = tagInfoIP.tracks().size();
	jets_.nselIPtrk[jets_.nref] = tagInfoIP.selectedTracks().size();
	if (doLifeTimeTaggingExtras_) {

	  TrackRefVector selTracks=tagInfoIP.selectedTracks();

	  GlobalPoint pv(tagInfoIP.primaryVertex()->position().x(),tagInfoIP.primaryVertex()->position().y(),tagInfoIP.primaryVertex()->position().z());

  	  reco::TrackKinematics allKinematics;	 
 
	  for(int it=0;it<jets_.nselIPtrk[jets_.nref] ;it++)
	    {
	      jets_.ipJetIndex[jets_.nIP + it]= jets_.nref;
	      reco::btag::TrackIPData data = tagInfoIP.impactParameterData()[it];
	      jets_.ipPt[jets_.nIP + it] = selTracks[it]->pt();
	      jets_.ipEta[jets_.nIP + it] = selTracks[it]->eta();
	      jets_.ipDxy[jets_.nIP + it] = selTracks[it]->dxy(tagInfoIP.primaryVertex()->position());
	      jets_.ipDz[jets_.nIP + it] = selTracks[it]->dz(tagInfoIP.primaryVertex()->position());
	      jets_.ipChi2[jets_.nIP + it] = selTracks[it]->normalizedChi2();
	      jets_.ipNHit[jets_.nIP + it] = selTracks[it]->numberOfValidHits();
	      jets_.ipNHitPixel[jets_.nIP + it] = selTracks[it]->hitPattern().numberOfValidPixelHits();
	      jets_.ipNHitStrip[jets_.nIP + it] = selTracks[it]->hitPattern().numberOfValidStripHits();
	      jets_.ipIsHitL1[jets_.nIP + it]  = selTracks[it]->hitPattern().hasValidHitInFirstPixelBarrel();
	      jets_.ipProb0[jets_.nIP + it] = tagInfoIP.probabilities(0)[it];
	      jets_.ip2d[jets_.nIP + it] = data.ip2d.value();
	      jets_.ip2dSig[jets_.nIP + it] = data.ip2d.significance();
	      jets_.ip3d[jets_.nIP + it] = data.ip3d.value();
	      jets_.ip3dSig[jets_.nIP + it] = data.ip3d.significance();
	      jets_.ipDist2Jet[jets_.nIP + it] = data.distanceToJetAxis.value();
	      jets_.ipClosest2Jet[jets_.nIP + it] = (data.closestToJetAxis - pv).mag();  //decay length
	    
	      //additional info for charm tagger dev
	      if(doExtraCTagging_){
		      math::XYZVector trackMom = selTracks[it]->momentum();
		      jets_.trackPtRel[jets_.nIP + it] = ROOT::Math::VectorUtil::Perp(trackMom,jetDir);
		      jets_.trackPPar[jets_.nIP + it] = jetDir.Dot(trackMom);
		      jets_.trackDeltaR[jets_.nIP + it] = ROOT::Math::VectorUtil::DeltaR(trackMom,jetDir); 
		      jets_.trackPtRatio[jets_.nIP + it] = ROOT::Math::VectorUtil::Perp(trackMom, jetDir)/ std::sqrt(trackMom.Mag2());
		      jets_.trackPParRatio[jets_.nIP + it] = jetDir.Dot(trackMom)/ std::sqrt(trackMom.Mag2());
		      const Track track = *(selTracks[it]);
		      allKinematics.add(track);
	      }
	    }

	  if(doExtraCTagging_){
		  jets_.trackSumJetDeltaR[jets_.nref] = ROOT::Math::VectorUtil::DeltaR(allKinematics.vectorSum(),jetDir);

		  //do some sorting to get the impact parameter of all tracks in descending order
		  jets_.trackSip2dSigAboveCharm[jets_.nref] = getAboveCharmThresh(selTracks, tagInfoIP, 1);
		  jets_.trackSip3dSigAboveCharm[jets_.nref] = getAboveCharmThresh(selTracks, tagInfoIP, 2);
		  jets_.trackSip2dValAboveCharm[jets_.nref] = getAboveCharmThresh(selTracks, tagInfoIP, 3);
		  jets_.trackSip3dValAboveCharm[jets_.nref] = getAboveCharmThresh(selTracks, tagInfoIP, 4);
	  }
	  jets_.nIP += jets_.nselIPtrk[jets_.nref];

	}
      }

      ith_tagged = this->TaggedJet(jet,jetTags_JB);
      if(ith_tagged >= 0) jets_.discr_probb[jets_.nref]  = (*jetTags_JB)[ith_tagged].second;

      ith_tagged    = this->TaggedJet(jet,jetTags_TCHighEff);
      if(ith_tagged >= 0) jets_.discr_tcHighEff[jets_.nref]   = (*jetTags_TCHighEff)[ith_tagged].second;
      ith_tagged    = this->TaggedJet(jet,jetTags_TCHighPur);
      if(ith_tagged >= 0) jets_.discr_tcHighPur[jets_.nref]   = (*jetTags_TCHighPur)[ith_tagged].second;

      //ith_tagged = this->TaggedJet(jet,jetTags_softMu);
      //if(ith_tagged >= 0){
      //  if ( (*jetTags_softMu)[ith_tagged].second     > -100000 )
      //    jets_.discr_muByPt[jets_.nref]  = (*jetTags_softMu)[ith_tagged].second;
      //}
      //ith_tagged = this->TaggedJet(jet,jetTags_softMuneg);
      //if(ith_tagged >= 0){
      //  float SoftMN = 0;
      //  if ( (*jetTags_softMuneg)[ith_tagged].second  > -100000 )
      //    SoftMN = ((*jetTags_softMuneg)[ith_tagged].second);
      //  if ( SoftMN > 0 ) SoftMN = -SoftMN;
      //  jets_.ndiscr_muByPt[jets_.nref] = SoftMN;
      //}

      const PFCandidateCollection *pfCandidateColl = &(*pfCandidates);
      int pfMuonIndex = getPFJetMuon(jet, pfCandidateColl);


      if(pfMuonIndex >=0){
	const reco::PFCandidate muon = pfCandidateColl->at(pfMuonIndex);
	jets_.mupt[jets_.nref]    =  muon.pt();
	jets_.mueta[jets_.nref]   =  muon.eta();
	jets_.muphi[jets_.nref]   =  muon.phi();
	jets_.mue[jets_.nref]     =  muon.energy();
	jets_.mudr[jets_.nref]    =  reco::deltaR(jet,muon);
	jets_.muptrel[jets_.nref] =  getPtRel(muon, jet);
	jets_.muchg[jets_.nref]   =  muon.charge();
      }else{
	jets_.mupt[jets_.nref]    =  0.0;
	jets_.mueta[jets_.nref]   =  0.0;
	jets_.muphi[jets_.nref]   =  0.0;
	jets_.mue[jets_.nref]     =  0.0;
	jets_.mudr[jets_.nref]    =  9.9;
	jets_.muptrel[jets_.nref] =  0.0;
	jets_.muchg[jets_.nref]   = 0;
      }
    }

    if(doHiJetID_){
      // Jet ID variables

      jets_.muMax[jets_.nref] = 0;
      jets_.muSum[jets_.nref] = 0;
      jets_.muN[jets_.nref] = 0;

      jets_.eMax[jets_.nref] = 0;
      jets_.eSum[jets_.nref] = 0;
      jets_.eN[jets_.nref] = 0;

      jets_.neutralMax[jets_.nref] = 0;
      jets_.neutralSum[jets_.nref] = 0;
      jets_.neutralN[jets_.nref] = 0;

      jets_.photonMax[jets_.nref] = 0;
      jets_.photonSum[jets_.nref] = 0;
      jets_.photonN[jets_.nref] = 0;
      jets_.photonHardSum[jets_.nref] = 0;
      jets_.photonHardN[jets_.nref] = 0;

      jets_.chargedMax[jets_.nref] = 0;
      jets_.chargedSum[jets_.nref] = 0;
      jets_.chargedN[jets_.nref] = 0;
      jets_.chargedHardSum[jets_.nref] = 0;
      jets_.chargedHardN[jets_.nref] = 0;

      jets_.trackMax[jets_.nref] = 0;
      jets_.trackSum[jets_.nref] = 0;
      jets_.trackN[jets_.nref] = 0;
      jets_.trackHardSum[jets_.nref] = 0;
      jets_.trackHardN[jets_.nref] = 0;

      jets_.hcalSum[jets_.nref] = 0;
      jets_.ecalSum[jets_.nref] = 0;

      jets_.genChargedSum[jets_.nref] = 0;
      jets_.genHardSum[jets_.nref] = 0;

      jets_.signalChargedSum[jets_.nref] = 0;
      jets_.signalHardSum[jets_.nref] = 0;

      jets_.subid[jets_.nref] = -1;

      for(unsigned int icand = 0; icand < tracks->size(); ++icand){
	const reco::Track& track = (*tracks)[icand];
	if(useQuality_ ){
	  bool goodtrack = track.quality(reco::TrackBase::qualityByName(trackQuality_));
	  if(!goodtrack) continue;
	}

	double dr = deltaR(jet,track);
	if(dr < rParam){
	  double ptcand = track.pt();
	  jets_.trackSum[jets_.nref] += ptcand;
	  jets_.trackN[jets_.nref] += 1;

	  if(ptcand > hardPtMin_){
	    jets_.trackHardSum[jets_.nref] += ptcand;
	    jets_.trackHardN[jets_.nref] += 1;

	  }
	  if(ptcand > jets_.trackMax[jets_.nref]) jets_.trackMax[jets_.nref] = ptcand;

	}
      }

      for(unsigned int icand = 0; icand < pfCandidates->size(); ++icand){
        const reco::PFCandidate& track = (*pfCandidates)[icand];
        double dr = deltaR(jet,track);
        if(dr < rParam){
	  double ptcand = track.pt();
	  int pfid = track.particleId();

	  switch(pfid){

	  case 1:
	    jets_.chargedSum[jets_.nref] += ptcand;
	    jets_.chargedN[jets_.nref] += 1;
	    if(ptcand > hardPtMin_){
	      jets_.chargedHardSum[jets_.nref] += ptcand;
	      jets_.chargedHardN[jets_.nref] += 1;
	    }
	    if(ptcand > jets_.chargedMax[jets_.nref]) jets_.chargedMax[jets_.nref] = ptcand;
	    break;

	  case 2:
	    jets_.eSum[jets_.nref] += ptcand;
	    jets_.eN[jets_.nref] += 1;
	    if(ptcand > jets_.eMax[jets_.nref]) jets_.eMax[jets_.nref] = ptcand;
	    break;

	  case 3:
	    jets_.muSum[jets_.nref] += ptcand;
	    jets_.muN[jets_.nref] += 1;
	    if(ptcand > jets_.muMax[jets_.nref]) jets_.muMax[jets_.nref] = ptcand;
	    break;

	  case 4:
	    jets_.photonSum[jets_.nref] += ptcand;
	    jets_.photonN[jets_.nref] += 1;
	    if(ptcand > hardPtMin_){
	      jets_.photonHardSum[jets_.nref] += ptcand;
	      jets_.photonHardN[jets_.nref] += 1;
	    }
	    if(ptcand > jets_.photonMax[jets_.nref]) jets_.photonMax[jets_.nref] = ptcand;
	    break;

	  case 5:
	    jets_.neutralSum[jets_.nref] += ptcand;
	    jets_.neutralN[jets_.nref] += 1;
	    if(ptcand > jets_.neutralMax[jets_.nref]) jets_.neutralMax[jets_.nref] = ptcand;
	    break;

	  default:
	    break;

	  }
	}
      }

      // Calorimeter fractions
      // Jet ID for CaloJets
      if(doTower){
	// changing it to take things from towers
	for(unsigned int i = 0; i < towers->size(); ++i){
	  
	  const CaloTower & hit= (*towers)[i];
	  double dr = deltaR(jet.eta(), jet.phi(), hit.p4(vtx).Eta(), hit.p4(vtx).Phi());
	  if(dr < rParam){
	    jets_.ecalSum[jets_.nref] += hit.emEt(vtx);
	    jets_.hcalSum[jets_.nref] += hit.hadEt(vtx);
	  }
	  
	}
      }

    }

    


    if(doMatch_){

      // Alternative reconstruction matching (PF for calo, calo for PF)

      double drMin = 100;
      for(unsigned int imatch = 0 ; imatch < matchedjets->size(); ++imatch){
	const reco::Jet& mjet = (*matchedjets)[imatch];

	double dr = deltaR(jet,mjet);
	if(dr < drMin){
	  jets_.matchedPt[jets_.nref] = mjet.pt();

	  if(usePat_){
	    const pat::Jet& mpatjet = (*patmatchedjets)[imatch];
	    jets_.matchedRawPt[jets_.nref] = mpatjet.correctedJet("Uncorrected").pt();
	    jets_.matchedPu[jets_.nref] = mpatjet.pileup();
	  }
	  jets_.matchedR[jets_.nref] = dr;
	  drMin = dr;
	}
      }

    }
    //     if(etrk.quality(reco::TrackBase::qualityByName(qualityString_))) pev_.trkQual[pev_.nTrk]=1;

    if(doHiJetID_){

      // JetID Selections for 5 TeV PbPb
      // Cuts by Yen-Jie Lee, BDT by Kurt Jung
      double rawpt = jets_.rawpt[jets_.nref];
      float jtpt = jet.pt();
      int ijet = jets_.nref;
      jets_.discr_jetID_cuts[jets_.nref] = jets_.neutralMax[ijet]/rawpt*0.085 + 
					   jets_.photonMax[ijet]/rawpt*0.337 + 
					   jets_.chargedMax[ijet]/rawpt*0.584 + 
					   jets_.neutralSum[ijet]/rawpt*-0.454 + 
					   jets_.photonSum[ijet]/rawpt*-0.127 + 
					   jets_.chargedSum[ijet]/rawpt*(-0.239) + 
					   jets_.jtpu[ijet]/rawpt*(-0.184) + 0.173;
      //begin bdt - TMVA requires you to load in all input vars and names into separate containers (eyeroll)
      float tempVarAddr[] = {jets_.trackMax[ijet]/jtpt, jets_.trackHardSum[ijet]/jtpt, jets_.trackHardN[ijet]/jtpt, jets_.chargedN[ijet]/jtpt, jets_.chargedHardSum[ijet]/jtpt, jets_.chargedHardN[ijet]/jtpt, jets_.photonN[ijet]/jtpt, jets_.photonHardSum[ijet]/jtpt, jets_.photonHardN[ijet]/jtpt, jets_.neutralN[ijet]/jtpt, jets_.hcalSum[ijet]/jtpt, jets_.ecalSum[ijet]/jtpt, jets_.chargedMax[ijet]/jtpt, jets_.chargedSum[ijet], jets_.neutralMax[ijet]/jtpt, jets_.neutralSum[ijet]/jtpt, jets_.photonMax[ijet]/jtpt, jets_.photonSum[ijet]/jtpt, jets_.eSum[ijet]/jtpt, jets_.muSum[ijet]/jtpt};
      for(unsigned int ivar=0; ivar<(sizeof(tempVarAddr)/sizeof(tempVarAddr[0])); ivar++){ varAddr[ivar] = tempVarAddr[ivar]; }
      jets_.discr_jetID_bdt[jets_.nref] = reader->EvaluateMVA("BDTG");

      /////////////////////////////////////////////////////////////////
      // Jet core pt^2 discriminant for fake jets
      // Edited by Yue Shi Lai <ylai@mit.edu>

      // Initial value is 0
      jets_.discr_fr01[jets_.nref] = 0;
      // Start with no directional adaption, i.e. the fake rejection
      // axis is the jet axis
      float pseudorapidity_adapt = jets_.jteta[jets_.nref];
      float azimuth_adapt = jets_.jtphi[jets_.nref];

      // Unadapted discriminant with adaption search
      for (size_t iteration = 0; iteration < 2; iteration++) {
	float pseudorapidity_adapt_new = pseudorapidity_adapt;
	float azimuth_adapt_new = azimuth_adapt;
	float max_weighted_perp = 0;
	float perp_square_sum = 0;

	for (size_t index_pf_candidate = 0;
	     index_pf_candidate < pfCandidates->size();
	     index_pf_candidate++) {
	  const reco::PFCandidate &p =
	    (*pfCandidates)[index_pf_candidate];

	  switch (p.particleId()) {
	    //case 1:	// Charged hadron
	    //case 3:	// Muon
	  case 4:	// Photon
	    {
	      const float dpseudorapidity =
		p.eta() - pseudorapidity_adapt;
	      const float dazimuth =
		reco::deltaPhi(p.phi(), azimuth_adapt);
	      // The Gaussian scale factor is 0.5 / (0.1 * 0.1)
	      // = 50
	      const float angular_weight =
		exp(-50.0F * (dpseudorapidity * dpseudorapidity +
			      dazimuth * dazimuth));
	      const float weighted_perp =
		angular_weight * p.pt() * p.pt();
	      const float weighted_perp_square =
		weighted_perp * p.pt();

	      perp_square_sum += weighted_perp_square;
	      if (weighted_perp >= max_weighted_perp) {
		pseudorapidity_adapt_new = p.eta();
		azimuth_adapt_new = p.phi();
		max_weighted_perp = weighted_perp;
	      }
	    }
	  default:
	    break;
	  }
	}
	// Update the fake rejection value
	jets_.discr_fr01[jets_.nref] = std::max(
						jets_.discr_fr01[jets_.nref], perp_square_sum);
	// Update the directional adaption
	pseudorapidity_adapt = pseudorapidity_adapt_new;
	azimuth_adapt = azimuth_adapt_new;
      }
    }
    jets_.jtpt[jets_.nref] = jet.pt();
    jets_.jteta[jets_.nref] = jet.eta();
    jets_.jtphi[jets_.nref] = jet.phi();
    jets_.jty[jets_.nref] = jet.eta();
    jets_.jtpu[jets_.nref] = jet.pileup();
    jets_.jtm[jets_.nref] = jet.mass();
    jets_.jtarea[jets_.nref] = jet.jetArea();

    jets_.jttau1[jets_.nref] = -999.;
    jets_.jttau2[jets_.nref] = -999.;
    jets_.jttau3[jets_.nref] = -999.;

    if(doSubJets_) analyzeSubjets(jet);

    if(usePat_){
      if( (*patjets)[j].hasUserFloat(jetName_+"Njettiness:tau1") )
        jets_.jttau1[jets_.nref] = (*patjets)[j].userFloat(jetName_+"Njettiness:tau1");
      if( (*patjets)[j].hasUserFloat(jetName_+"Njettiness:tau2") )
        jets_.jttau2[jets_.nref] = (*patjets)[j].userFloat(jetName_+"Njettiness:tau2");
      if( (*patjets)[j].hasUserFloat(jetName_+"Njettiness:tau3") )
        jets_.jttau3[jets_.nref] = (*patjets)[j].userFloat(jetName_+"Njettiness:tau3");

      if( (*patjets)[j].isPFJet()) {
        jets_.jtPfCHF[jets_.nref] = (*patjets)[j].chargedHadronEnergyFraction();
        jets_.jtPfNHF[jets_.nref] = (*patjets)[j].neutralHadronEnergyFraction();
        jets_.jtPfCEF[jets_.nref] = (*patjets)[j].chargedEmEnergyFraction();
        jets_.jtPfNEF[jets_.nref] = (*patjets)[j].neutralEmEnergyFraction();
        jets_.jtPfMUF[jets_.nref] = (*patjets)[j].muonEnergyFraction();

        jets_.jtPfCHM[jets_.nref] = (*patjets)[j].chargedHadronMultiplicity();
        jets_.jtPfNHM[jets_.nref] = (*patjets)[j].neutralHadronMultiplicity();
        jets_.jtPfCEM[jets_.nref] = (*patjets)[j].electronMultiplicity();
        jets_.jtPfNEM[jets_.nref] = (*patjets)[j].photonMultiplicity();
        jets_.jtPfMUM[jets_.nref] = (*patjets)[j].muonMultiplicity();
      }
        
      if(doStandardJetID_){
	jets_.fHPD[jets_.nref] = (*patjets)[j].jetID().fHPD;
	jets_.fRBX[jets_.nref] = (*patjets)[j].jetID().fRBX;
	jets_.n90[jets_.nref] = (*patjets)[j].n90();

	jets_.fSubDet1[jets_.nref] = (*patjets)[j].jetID().fSubDetector1;
	jets_.fSubDet2[jets_.nref] = (*patjets)[j].jetID().fSubDetector2;
	jets_.fSubDet3[jets_.nref] = (*patjets)[j].jetID().fSubDetector3;
	jets_.fSubDet4[jets_.nref] = (*patjets)[j].jetID().fSubDetector4;
	jets_.restrictedEMF[jets_.nref] = (*patjets)[j].jetID().restrictedEMF;
	jets_.nHCAL[jets_.nref] = (*patjets)[j].jetID().nHCALTowers;
	jets_.nECAL[jets_.nref] = (*patjets)[j].jetID().nECALTowers;
	jets_.apprHPD[jets_.nref] = (*patjets)[j].jetID().approximatefHPD;
	jets_.apprRBX[jets_.nref] = (*patjets)[j].jetID().approximatefRBX;

	//       jets_.n90[jets_.nref] = (*patjets)[j].jetID().hitsInN90;
	jets_.n2RPC[jets_.nref] = (*patjets)[j].jetID().numberOfHits2RPC;
	jets_.n3RPC[jets_.nref] = (*patjets)[j].jetID().numberOfHits3RPC;
	jets_.nRPC[jets_.nref] = (*patjets)[j].jetID().numberOfHitsRPC;

	jets_.fEB[jets_.nref] = (*patjets)[j].jetID().fEB;
	jets_.fEE[jets_.nref] = (*patjets)[j].jetID().fEE;
	jets_.fHB[jets_.nref] = (*patjets)[j].jetID().fHB;
	jets_.fHE[jets_.nref] = (*patjets)[j].jetID().fHE;
	jets_.fHO[jets_.nref] = (*patjets)[j].jetID().fHO;
	jets_.fLong[jets_.nref] = (*patjets)[j].jetID().fLong;
	jets_.fShort[jets_.nref] = (*patjets)[j].jetID().fShort;
	jets_.fLS[jets_.nref] = (*patjets)[j].jetID().fLS;
	jets_.fHFOOT[jets_.nref] = (*patjets)[j].jetID().fHFOOT;
      }

    }

    if(isMC_){

      for(UInt_t i = 0; i < genparts->size(); ++i){
	const reco::GenParticle& p = (*genparts)[i];
	if (p.status()!=1) continue;
	if (p.charge()==0) continue;
	double dr = deltaR(jet,p);
	if(dr < rParam){
	  double ppt = p.pt();
	  jets_.genChargedSum[jets_.nref] += ppt;
	  if(ppt > hardPtMin_) jets_.genHardSum[jets_.nref] += ppt;
	  if(p.collisionId() == 0){
	    jets_.signalChargedSum[jets_.nref] += ppt;
	    if(ppt > hardPtMin_) jets_.signalHardSum[jets_.nref] += ppt;
	  }

	}
      }

    }

    if(isMC_ && usePat_){


      const reco::GenJet * genjet = (*patjets)[j].genJet();

      if(genjet){
	jets_.refpt[jets_.nref] = genjet->pt();
	jets_.refeta[jets_.nref] = genjet->eta();
	jets_.refphi[jets_.nref] = genjet->phi();
        jets_.refm[jets_.nref] = genjet->mass();
        jets_.refarea[jets_.nref] = genjet->jetArea();
        jets_.refy[jets_.nref] = genjet->eta();
	jets_.refdphijt[jets_.nref] = reco::deltaPhi(jet.phi(), genjet->phi());
	jets_.refdrjt[jets_.nref] = reco::deltaR(jet.eta(),jet.phi(),genjet->eta(),genjet->phi());

	if(doSubEvent_){
	  const GenParticle* gencon = genjet->getGenConstituent(0);
	  jets_.subid[jets_.nref] = gencon->collisionId();
	}

        if(doGenSubJets_) analyzeRefSubjets(*genjet);

      }else{
	jets_.refpt[jets_.nref] = -999.;
	jets_.refeta[jets_.nref] = -999.;
	jets_.refphi[jets_.nref] = -999.;
        jets_.refm[jets_.nref] = -999.;
        jets_.refarea[jets_.nref] = -999.;
	jets_.refy[jets_.nref] = -999.;
	jets_.refdphijt[jets_.nref] = -999.;
	jets_.refdrjt[jets_.nref] = -999.;

        if(doGenSubJets_) {
          jets_.refptG[jets_.nref]  = -999.;
          jets_.refetaG[jets_.nref] = -999.;
          jets_.refphiG[jets_.nref] = -999.;
          jets_.refmG[jets_.nref]   = -999.;

          std::vector<float> sjpt;
          std::vector<float> sjeta;
          std::vector<float> sjphi;
          std::vector<float> sjm;
          sjpt.push_back(-999.);
          sjeta.push_back(-999.);
          sjphi.push_back(-999.);
          sjm.push_back(-999.);

          jets_.refSubJetPt.push_back(sjpt);
          jets_.refSubJetEta.push_back(sjeta);
          jets_.refSubJetPhi.push_back(sjphi);
          jets_.refSubJetM.push_back(sjm);
        }
      }
      jets_.reftau1[jets_.nref] = -999.;
      jets_.reftau2[jets_.nref] = -999.;
      jets_.reftau3[jets_.nref] = -999.;
      
      jets_.refparton_flavorForB[jets_.nref] = (*patjets)[j].partonFlavour();
      // matched partons
      const reco::GenParticle & parton = *(*patjets)[j].genParton();

      if((*patjets)[j].genParton()){
	jets_.refparton_pt[jets_.nref] = parton.pt();
	jets_.refparton_flavor[jets_.nref] = parton.pdgId();

	if(saveBfragments_ && abs(jets_.refparton_flavorForB[jets_.nref])==5){

	  usedStringPts.clear();

	  // uncomment this if you want to know the ugly truth about parton matching -matt
	  //if(jet.pt() > 50 &&abs(parton.pdgId())!=5 && parton.pdgId()!=21)
	  // cout<<" Identified as a b, but doesn't match b or gluon, id = "<<parton.pdgId()<<endl;

	  jets_.bJetIndex[jets_.bMult] = jets_.nref;
	  jets_.bStatus[jets_.bMult] = parton.status();
	  jets_.bVx[jets_.bMult] = parton.vx();
	  jets_.bVy[jets_.bMult] = parton.vy();
	  jets_.bVz[jets_.bMult] = parton.vz();
	  jets_.bPt[jets_.bMult] = parton.pt();
	  jets_.bEta[jets_.bMult] = parton.eta();
	  jets_.bPhi[jets_.bMult] = parton.phi();
	  jets_.bPdg[jets_.bMult] = parton.pdgId();
	  jets_.bChg[jets_.bMult] = parton.charge();
	  jets_.bMult++;
	  saveDaughters(parton);
	}
      } else {
	jets_.refparton_pt[jets_.nref] = -999;
	jets_.refparton_flavor[jets_.nref] = -999;
      }
    }
    jets_.nref++;
  }

  if(isMC_){
  
    if(useHepMC_) {
      edm::Handle<HepMCProduct> hepMCProduct;
      iEvent.getByToken(eventInfoTag_,hepMCProduct);
      const HepMC::GenEvent* MCEvt = hepMCProduct->GetEvent();

      std::pair<HepMC::GenParticle*,HepMC::GenParticle*> beamParticles = MCEvt->beam_particles();
      if(beamParticles.first != 0)jets_.beamId1 = beamParticles.first->pdg_id();
      if(beamParticles.second != 0)jets_.beamId2 = beamParticles.second->pdg_id();
    }
    
    edm::Handle<GenEventInfoProduct> hEventInfo;
    iEvent.getByToken(eventGenInfoTag_,hEventInfo);
    //jets_.pthat = hEventInfo->binningValues()[0];

    // binning values and qscale appear to be equivalent, but binning values not always present
    jets_.pthat = hEventInfo->qScale();

    edm::Handle<vector<reco::GenJet> >genjets;
    //edm::Handle<edm::View<reco::Jet>> genjets;
    iEvent.getByToken(genjetTag_, genjets);
    
    //get gen-level n-jettiness
    // edm::Handle<edm::ValueMap<float> > genTau1s;
    // edm::Handle<edm::ValueMap<float> > genTau2s;
    // edm::Handle<edm::ValueMap<float> > genTau3s;
    // if(useGenTaus) {
    //   Printf("get gen taus");
    //   iEvent.getByToken(tokenGenTau1_,genTau1s);
    //   iEvent.getByToken(tokenGenTau2_,genTau2s);
    //   iEvent.getByToken(tokenGenTau3_,genTau3s);
    // }
    
    jets_.ngen = 0;

    //int igen = 0;
    for(unsigned int igen = 0 ; igen < genjets->size(); ++igen){
      //for ( typename edm::View<reco::Jet>::const_iterator genjetIt = genjets->begin() ; genjetIt != genjets->end() ; ++genjetIt ) {
      const reco::GenJet & genjet = (*genjets)[igen];
      //edm::Ptr<reco::Jet> genjetPtr = genjets->ptrAt(genjetIt - genjets->begin());
      //const reco::GenJet genjet = (*dynamic_cast<const reco::GenJet*>(&(*genjetPtr)));
      float genjet_pt = genjet.pt();

      float tau1 =  -999.;
      float tau2 =  -999.;
      float tau3 =  -9999.;
      if(doGenTaus_) {
        tau1 =  getTau(1,genjet);
        tau2 =  getTau(2,genjet);
        tau3 =  getTau(3,genjet);
      }

      // find matching patJet if there is one
      jets_.gendrjt[jets_.ngen] = -1.0;
      jets_.genmatchindex[jets_.ngen] = -1;
      
      for(int ijet = 0 ; ijet < jets_.nref; ++ijet){
        // poor man's matching, someone fix please
        if(fabs(genjet.pt()-jets_.refpt[ijet])<0.00001 &&
           fabs(genjet.eta()-jets_.refeta[ijet])<0.00001){
          if(genjet_pt>genPtMin_) {
            jets_.genmatchindex[jets_.ngen] = (int)ijet;
            jets_.gendphijt[jets_.ngen] = reco::deltaPhi(jets_.refphi[ijet],genjet.phi());
            jets_.gendrjt[jets_.ngen] = sqrt(pow(jets_.gendphijt[jets_.ngen],2)+pow(fabs(genjet.eta()-jets_.refeta[ijet]),2));
          }
          if(doGenTaus_) {
            jets_.reftau1[ijet] = tau1;
            jets_.reftau2[ijet] = tau2;
            jets_.reftau3[ijet] = tau3;
          }
          break;
        }
      }
      
      // threshold to reduce size of output in minbias PbPb
      if(genjet_pt>genPtMin_){
	jets_.genpt[jets_.ngen] = genjet_pt;
	jets_.geneta[jets_.ngen] = genjet.eta();
	jets_.genphi[jets_.ngen] = genjet.phi();
        jets_.genm[jets_.ngen] = genjet.mass();
	jets_.geny[jets_.ngen] = genjet.eta();

        if(doGenTaus_) {
          jets_.gentau1[jets_.ngen] = tau1;
          jets_.gentau2[jets_.ngen] = tau2;
          jets_.gentau3[jets_.ngen] = tau3;
        }

        if(doGenSubJets_)
          analyzeGenSubjets(genjet);

	if(doSubEvent_){
	  const GenParticle* gencon = genjet.getGenConstituent(0);
	  jets_.gensubid[jets_.ngen] = gencon->collisionId();
	}
	jets_.ngen++;
      }
    }

  }


  t->Fill();
  memset(&jets_,0,sizeof jets_);
}


//--------------------------------------------------------------------------------------------------
void HiInclusiveJetAnalyzer::fillL1Bits(const edm::Event &iEvent)
{
  edm::Handle< L1GlobalTriggerReadoutRecord >  L1GlobalTrigger;

  iEvent.getByToken(L1gtReadout_, L1GlobalTrigger);
  const TechnicalTriggerWord&  technicalTriggerWordBeforeMask = L1GlobalTrigger->technicalTriggerWord();

  for (int i=0; i<64;i++)
    {
      jets_.l1TBit[i] = technicalTriggerWordBeforeMask.at(i);
    }
  jets_.nL1TBit = 64;

  int ntrigs = L1GlobalTrigger->decisionWord().size();
  jets_.nL1ABit = ntrigs;

  for (int i=0; i != ntrigs; i++) {
    bool accept = L1GlobalTrigger->decisionWord()[i];
    //jets_.l1ABit[i] = (accept == true)? 1:0;
    if(accept== true){
      jets_.l1ABit[i] = 1;
    }
    else{
      jets_.l1ABit[i] = 0;
    }

  }
}

//--------------------------------------------------------------------------------------------------
void HiInclusiveJetAnalyzer::fillHLTBits(const edm::Event &iEvent)
{
  // Fill HLT trigger bits.
  Handle<TriggerResults> triggerResultsHLT;
  iEvent.getByToken(hltResName_, triggerResultsHLT);

  const TriggerResults *hltResults = triggerResultsHLT.product();
  const TriggerNames & triggerNames = iEvent.triggerNames(*hltResults);

  jets_.nHLTBit = hltTrgNames_.size();

  for(size_t i=0;i<hltTrgNames_.size();i++){

    for(size_t j=0;j<triggerNames.size();++j) {

      if(triggerNames.triggerName(j) == hltTrgNames_[i]){

	//cout <<"hltTrgNames_(i) "<<hltTrgNames_[i]<<endl;
	//cout <<"triggerName(j) "<<triggerNames.triggerName(j)<<endl;
	//cout<<" result "<<triggerResultsHLT->accept(j)<<endl;
	jets_.hltBit[i] = triggerResultsHLT->accept(j);
      }

    }
  }
}

int
HiInclusiveJetAnalyzer::getPFJetMuon(const pat::Jet& pfJet, const reco::PFCandidateCollection *pfCandidateColl)
{

  int pfMuonIndex = -1;
  float ptMax = 0.;


  for(unsigned icand=0;icand<pfCandidateColl->size(); icand++) {
    const reco::PFCandidate pfCandidate = pfCandidateColl->at(icand);

    int id = pfCandidate.particleId();
    if(abs(id) != 3) continue;

    if(reco::deltaR(pfJet,pfCandidate)>0.5) continue;

    double pt =  pfCandidate.pt();
    if(pt>ptMax){
      ptMax = pt;
      pfMuonIndex = (int) icand;
    }
  }

  return pfMuonIndex;

}


double
HiInclusiveJetAnalyzer::getPtRel(const reco::PFCandidate lep, const pat::Jet& jet )
{

  float lj_x = jet.p4().px();
  float lj_y = jet.p4().py();
  float lj_z = jet.p4().pz();

  // absolute values squared
  float lj2  = lj_x*lj_x+lj_y*lj_y+lj_z*lj_z;
  float lep2 = lep.px()*lep.px()+lep.py()*lep.py()+lep.pz()*lep.pz();

  // projection vec(mu) to lepjet axis
  float lepXlj = lep.px()*lj_x+lep.py()*lj_y+lep.pz()*lj_z;

  // absolute value squared and normalized
  float pLrel2 = lepXlj*lepXlj/lj2;

  // lep2 = pTrel2 + pLrel2
  float pTrel2 = lep2-pLrel2;

  return (pTrel2 > 0) ? std::sqrt(pTrel2) : 0.0;
}

// Recursive function, but this version gets called only the first time
void
HiInclusiveJetAnalyzer::saveDaughters(const reco::GenParticle &gen){

  for(unsigned i=0;i<gen.numberOfDaughters();i++){
    const reco::Candidate & daughter = *gen.daughter(i);
    double daughterPt = daughter.pt();
    if(daughterPt<1.) continue;
    double daughterEta = daughter.eta();
    if(fabs(daughterEta)>3.) continue;
    int daughterPdgId = daughter.pdgId();
    int daughterStatus = daughter.status();
    // Special case when b->b+string, both b and string contain all daughters, so only take the string
    if(gen.pdgId()==daughterPdgId && gen.status()==3 && daughterStatus==2) continue;

    // cheesy way of finding strings which were already used
    if(daughter.pdgId()==92){
      for(unsigned ist=0;ist<usedStringPts.size();ist++){
	if(fabs(daughter.pt() - usedStringPts[ist]) < 0.0001) return;
      }
      usedStringPts.push_back(daughter.pt());
    }
    jets_.bJetIndex[jets_.bMult] = jets_.nref;
    jets_.bStatus[jets_.bMult] = daughterStatus;
    jets_.bVx[jets_.bMult] = daughter.vx();
    jets_.bVy[jets_.bMult] = daughter.vy();
    jets_.bVz[jets_.bMult] = daughter.vz();
    jets_.bPt[jets_.bMult] = daughterPt;
    jets_.bEta[jets_.bMult] = daughterEta;
    jets_.bPhi[jets_.bMult] = daughter.phi();
    jets_.bPdg[jets_.bMult] = daughterPdgId;
    jets_.bChg[jets_.bMult] = daughter.charge();
    jets_.bMult++;
    saveDaughters(daughter);
  }
}

// This version called for all subsequent calls
void
HiInclusiveJetAnalyzer::saveDaughters(const reco::Candidate &gen){

  for(unsigned i=0;i<gen.numberOfDaughters();i++){
    const reco::Candidate & daughter = *gen.daughter(i);
    double daughterPt = daughter.pt();
    if(daughterPt<1.) continue;
    double daughterEta = daughter.eta();
    if(fabs(daughterEta)>3.) continue;
    int daughterPdgId = daughter.pdgId();
    int daughterStatus = daughter.status();
    // Special case when b->b+string, both b and string contain all daughters, so only take the string
    if(gen.pdgId()==daughterPdgId && gen.status()==3 && daughterStatus==2) continue;

    // cheesy way of finding strings which were already used
    if(daughter.pdgId()==92){
      for(unsigned ist=0;ist<usedStringPts.size();ist++){
	if(fabs(daughter.pt() - usedStringPts[ist]) < 0.0001) return;
      }
      usedStringPts.push_back(daughter.pt());
    }

    jets_.bJetIndex[jets_.bMult] = jets_.nref;
    jets_.bStatus[jets_.bMult] = daughterStatus;
    jets_.bVx[jets_.bMult] = daughter.vx();
    jets_.bVy[jets_.bMult] = daughter.vy();
    jets_.bVz[jets_.bMult] = daughter.vz();
    jets_.bPt[jets_.bMult] = daughterPt;
    jets_.bEta[jets_.bMult] = daughterEta;
    jets_.bPhi[jets_.bMult] = daughter.phi();
    jets_.bPdg[jets_.bMult] = daughterPdgId;
    jets_.bChg[jets_.bMult] = daughter.charge();
    jets_.bMult++;
    saveDaughters(daughter);
  }
}

double HiInclusiveJetAnalyzer::getEt(math::XYZPoint pos, double energy){
  double et = energy*sin(pos.theta());
  return et;
}

math::XYZPoint HiInclusiveJetAnalyzer::getPosition(const DetId &id, reco::Vertex::Point vtx){
  const GlobalPoint& pos=geo->getPosition(id);
  math::XYZPoint posV(pos.x() - vtx.x(),pos.y() - vtx.y(),pos.z() - vtx.z());
  return posV;
}

int HiInclusiveJetAnalyzer::TaggedJet(Jet calojet, Handle<JetTagCollection > jetTags ) {
  double small = 1.e-5;
  int result = -1;

  for (size_t t = 0; t < jetTags->size(); t++) {
    RefToBase<Jet> jet_p = (*jetTags)[t].first;
    if (jet_p.isNull()) {
      continue;
    }
    if (DeltaR<Candidate>()( calojet, *jet_p ) < small) {
      result = (int) t;
    }
  }
  return result;
}

float HiInclusiveJetAnalyzer::getTau(unsigned num, const reco::GenJet object) const
{
  std::vector<fastjet::PseudoJet> FJparticles;
  if(object.numberOfDaughters()>0) {
    for (unsigned k = 0; k < object.numberOfDaughters(); ++k)
      {
        const reco::Candidate & dp = *object.daughter(k);
        //const reco::CandidatePtr & dp = object.daughterPtr(k);
        if(dp.pt()>0.) FJparticles.push_back( fastjet::PseudoJet( dp.px(), dp.py(), dp.pz(), dp.energy() ) );
      }
  }
  return routine_->getTau(num, FJparticles);
}

void HiInclusiveJetAnalyzer::analyzeSubjets(const reco::Jet jet) {

  std::vector<float> sjpt;
  std::vector<float> sjeta;
  std::vector<float> sjphi;
  std::vector<float> sjm;
  if(jet.numberOfDaughters()>0) {
    for (unsigned k = 0; k < jet.numberOfDaughters(); ++k) {
      const reco::Candidate & dp = *jet.daughter(k);
      sjpt.push_back(dp.pt());
      sjeta.push_back(dp.eta());
      sjphi.push_back(dp.phi());
      sjm.push_back(dp.mass());
    }
  } else {
    sjpt.push_back(-999.);
    sjeta.push_back(-999.);
    sjphi.push_back(-999.);
    sjm.push_back(-999.);
  }
  jets_.jtSubJetPt.push_back(sjpt);
  jets_.jtSubJetEta.push_back(sjeta);
  jets_.jtSubJetPhi.push_back(sjphi);
  jets_.jtSubJetM.push_back(sjm);
  
}

//--------------------------------------------------------------------------------------------------
int HiInclusiveJetAnalyzer::getGroomedGenJetIndex(const reco::GenJet jet) const {
  //Find closest soft-dropped gen jet
  double drMin = 100;
  int imatch = -1;
  for(unsigned int i = 0 ; i < gensubjets_->size(); ++i) {
    const reco::Jet& mjet = (*gensubjets_)[i];

    double dr = deltaR(jet,mjet);
    if(dr < drMin){
      imatch = i;
      drMin = dr;
    }
  }
  return imatch;
}

//--------------------------------------------------------------------------------------------------
void HiInclusiveJetAnalyzer::analyzeRefSubjets(const reco::GenJet jet) {

  //Find closest soft-dropped gen jet
  int imatch = getGroomedGenJetIndex(jet);
  double dr = 999.;
  if(imatch>-1) {
    const reco::Jet& mjet =  (*gensubjets_)[imatch];
    dr = deltaR(jet,mjet);
  }

  std::vector<float> sjpt;
  std::vector<float> sjeta;
  std::vector<float> sjphi;
  std::vector<float> sjm;
  if(imatch>-1 && dr<0.4) {
    const reco::Jet& mjet =  (*gensubjets_)[imatch];
    jets_.refptG[jets_.nref]  = mjet.pt();
    jets_.refetaG[jets_.nref] = mjet.eta();
    jets_.refphiG[jets_.nref] = mjet.phi();
    jets_.refmG[jets_.nref]   = mjet.mass();

    if(mjet.numberOfDaughters()>0) {
      for (unsigned k = 0; k < mjet.numberOfDaughters(); ++k) {
        const reco::Candidate & dp = *mjet.daughter(k);
        sjpt.push_back(dp.pt());
        sjeta.push_back(dp.eta());
        sjphi.push_back(dp.phi());
        sjm.push_back(dp.mass());
      }
    }
  }
  else {
    jets_.refptG[jets_.nref]  = -999.;
    jets_.refetaG[jets_.nref] = -999.;
    jets_.refphiG[jets_.nref] = -999.;
    jets_.refmG[jets_.nref]   = -999.;

    sjpt.push_back(-999.);
    sjeta.push_back(-999.);
    sjphi.push_back(-999.);
    sjm.push_back(-999.);
  }

  jets_.refSubJetPt.push_back(sjpt);
  jets_.refSubJetEta.push_back(sjeta);
  jets_.refSubJetPhi.push_back(sjphi);
  jets_.refSubJetM.push_back(sjm);
}

//--------------------------------------------------------------------------------------------------
void HiInclusiveJetAnalyzer::analyzeGenSubjets(const reco::GenJet jet) {
  //Find closest soft-dropped gen jet
  int imatch = getGroomedGenJetIndex(jet);
  double dr = 999.;
  if(imatch>-1) {
    const reco::Jet& mjet =  (*gensubjets_)[imatch];
    dr = deltaR(jet,mjet);
  }

  std::vector<float> sjpt;
  std::vector<float> sjeta;
  std::vector<float> sjphi;
  std::vector<float> sjm;
  std::vector<float> sjarea;
  if(imatch>-1 && dr<0.4) {
    const reco::Jet& mjet =  (*gensubjets_)[imatch];
    jets_.genptG[jets_.ngen]  = mjet.pt();
    jets_.genetaG[jets_.ngen] = mjet.eta();
    jets_.genphiG[jets_.ngen] = mjet.phi();
    jets_.genmG[jets_.ngen]   = mjet.mass();
    
    if(mjet.numberOfDaughters()>0) {
      for (unsigned k = 0; k < mjet.numberOfDaughters(); ++k) {
        const reco::Candidate & dp = *mjet.daughter(k);
        sjpt.push_back(dp.pt());
        sjeta.push_back(dp.eta());
        sjphi.push_back(dp.phi());
        sjm.push_back(dp.mass());
        //sjarea.push_back(dp.castTo<reco::JetRef>()->jetArea());
      }
    }
  }
  else {
    jets_.genptG[jets_.ngen]  = -999.;
    jets_.genetaG[jets_.ngen] = -999.;
    jets_.genphiG[jets_.ngen] = -999.;
    jets_.genmG[jets_.ngen]   = -999.;
    
    sjpt.push_back(-999.);
    sjeta.push_back(-999.);
    sjphi.push_back(-999.);
    sjm.push_back(-999.);
    sjarea.push_back(-999.);
  }

  jets_.genSubJetPt.push_back(sjpt);
  jets_.genSubJetEta.push_back(sjeta);
  jets_.genSubJetPhi.push_back(sjphi);
  jets_.genSubJetM.push_back(sjm);
  jets_.genSubJetArea.push_back(sjarea);
}

//--------------------------------------------------------------------------------------------------
float HiInclusiveJetAnalyzer::getAboveCharmThresh(TrackRefVector& selTracks, const TrackIPTagInfo& ipData, int sigOrVal)
{

	const double pdgCharmMass = 1.290;
	btag::SortCriteria sc;
	switch(sigOrVal){
		case 1: //2d significance 
			sc = reco::btag::IP2DSig;
		case 2: //3d significance
			sc = reco::btag::IP3DSig;
		case 3:
			sc = reco::btag::IP2DSig;  //values are not sortable!
		case 4:
			sc = reco::btag::IP3DSig;
	}
	std::vector<std::size_t> indices = ipData.sortedIndexes(sc);
	reco::TrackKinematics kin;
	for(unsigned int i=0; i<indices.size(); i++){
		size_t idx = indices[i];
		const Track track = *(selTracks[idx]);
		const btag::TrackIPData &data = ipData.impactParameterData()[idx];
		kin.add(track);
		if(kin.vectorSum().M() > pdgCharmMass){
			switch(sigOrVal){
				case 1:
					return data.ip2d.significance();
				case 2:
					return data.ip3d.significance();
				case 3:
					return data.ip2d.value();
				case 4:
					return data.ip3d.value();
			}
		}
	}
	return 0;

}


DEFINE_FWK_MODULE(HiInclusiveJetAnalyzer);
