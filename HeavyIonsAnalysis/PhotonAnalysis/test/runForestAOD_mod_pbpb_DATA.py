### HiForest Configuration
# Collisions: pp
# Type: Data
# Input: AOD

import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
process = cms.Process('HiForest')
process.options = cms.untracked.PSet()

# setup 'analysis'  options
options = VarParsing.VarParsing ('analysis')

# setup any defaults you want
options.inputFiles = '/store/group/phys_diffraction/diphoton/aug_reco_data_check_for_lumi/HIForward/Aug_reco_pbpb_2015_try2/170824_122451/0001/RECO_FILTER_RAW2DIGI_L1Reco_AOD_double_eg2_trigger_only_1001.root'

#options.outputFile = 'hiforest_aug_reco.root'
options.outputFile = 'hiforest.root'

options.parseArguments()

#####################################################################################
# HiForest labelling info
#####################################################################################

process.load("HeavyIonsAnalysis.JetAnalysis.HiForest_cff")
process.HiForest.inputLines = cms.vstring("HiForest V3",)
import subprocess
version = subprocess.Popen(["(cd $CMSSW_BASE/src && git describe --tags)"], stdout=subprocess.PIPE, shell=True).stdout.read()
if version == '':
    version = 'no git info'
process.HiForest.HiForestVersion = cms.string(version)

#####################################################################################
# Input source
#####################################################################################
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        options.inputFiles
    )
)
process.options = cms.untracked.PSet(
        wantSummary = cms.untracked.bool(True),
        # fileMode = cms.untracked.string('MERGE'),
)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

# Number of events we want to process, -1 = all events
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1))
    #input = cms.untracked.int32(1))


#####################################################################################
# Load Global Tag, Geometry, etc.
#####################################################################################

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
#process.load('FWCore.MessageService.MessageLogger_cfi')

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')
process.GlobalTag = GlobalTag(process.GlobalTag, '75X_dataRun2_PromptHI_v3', '')  #for now track GT manually, since centrality tables updated ex post facto
process.HiForest.GlobalTagLabel = process.GlobalTag.globaltag

from HeavyIonsAnalysis.Configuration.CommonFunctions_cff import overrideJEC_pp5020
process = overrideJEC_pp5020(process)

#####################################################################################
# Define tree output
#####################################################################################

process.TFileService = cms.Service("TFileService",
                                   fileName=cms.string(options.outputFile))

#####################################################################################
# Additional Reconstruction and Analysis: Main Body
#####################################################################################

####################################################################################

#############################
# Jets
#############################

process.highPurityTracks = cms.EDFilter("TrackSelector",
                                        src = cms.InputTag("generalTracks"),
                                        cut = cms.string('quality("highPurity")'))


#####################################################################################

############################
# Event Analysis
###########################

import HLTrigger.HLTcore.hltEventAnalyzerAOD_cfi
process.hltAOD = HLTrigger.HLTcore.hltEventAnalyzerAOD_cfi.hltEventAnalyzerAOD.clone()
process.hltAOD.processName = cms.string("HLT")
process.hltAOD.triggerResults = cms.InputTag("TriggerResults","","HLT")
process.hltAOD.triggerEvent   = cms.InputTag("hltTriggerSummaryAOD","","HLT")


process.load('HeavyIonsAnalysis.EventAnalysis.hievtanalyzer_data_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.hltobject_PbPb_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.hltanalysis_cff')
#triggerL1Name = cms.untracked.string("hltPreHIUPCL1DoubleEG2NotHF2"),
#triggerL1Name = cms.untracked.string("hltL1sL1DoubleEG2NotMinimumBiasHF2AND"),
#triggerL1Name = cms.untracked.string("hltL1sL1EG5NotMinimumBiasHF2AND"),

from HeavyIonsAnalysis.EventAnalysis.dummybranches_cff import addHLTdummybranches
addHLTdummybranches(process)

process.hiEvtAnalyzer.Vertex = cms.InputTag("offlinePrimaryVertices")
process.hiEvtAnalyzer.doCentrality = cms.bool(False)
process.hiEvtAnalyzer.doEvtPlane = cms.bool(False)

process.load("HeavyIonsAnalysis.JetAnalysis.pfcandAnalyzer_cfi")
process.pfcandAnalyzer.skipCharged = False
process.pfcandAnalyzer.pfPtMin = 0
process.load("HeavyIonsAnalysis.JetAnalysis.hcalNoise_cff")

#####################################################################################

#########################
# Track Analyzer
#########################
process.load('HeavyIonsAnalysis.JetAnalysis.ExtraTrackReco_cff')
process.load('HeavyIonsAnalysis.JetAnalysis.TrkAnalyzers_cff')
# process.load("HeavyIonsAnalysis.TrackAnalysis.METAnalyzer_cff")


####################################################################################

#####################
# Photons
#####################
process.load('HeavyIonsAnalysis.PhotonAnalysis.ggHiNtuplizer_cfi')
process.ggHiNtuplizer.doGenParticles = False
process.ggHiNtuplizerGED = process.ggHiNtuplizer.clone(recoPhotonSrc = cms.InputTag('gedPhotons'),
                                                       recoPhotonHiIsolationMap = cms.InputTag('photonIsolationHIProducerppGED')

)


####################################################################################

#####################
# tupel and necessary PAT sequences
#####################

process.load("HeavyIonsAnalysis.VectorBosonAnalysis.tupelSequence_pp_cff")

#####################################################################################

#########################
# Main analysis list
#########################
'''
process.ana_step = cms.Path(process.hltanalysis *
			    process.hltobject *
                            process.hiEvtAnalyzer*
                            #process.jetSequences +
                            process.ggHiNtuplizer +
                            process.ggHiNtuplizerGED +
                            process.pfcandAnalyzer +
                            #process.HiForest +
                            #process.trackSequencesPbPb +
                            process.trackSequencesPP #+
                            #process.hcalNoise #+
                            #process.tupelPatSequence
                            )
'''
#####################################################################################

#########################
# Event Selection
#########################

process.load('HeavyIonsAnalysis.JetAnalysis.EventSelection_cff')
#process.pcollisionEventSelection = cms.Path(process.collisionEventSelectionAOD)
#process.pHBHENoiseFilterResultProducer = cms.Path( process.HBHENoiseFilterResultProducer )
#process.HBHENoiseFilterResult = cms.Path(process.fHBHENoiseFilterResult)
#process.HBHENoiseFilterResultRun1 = cms.Path(process.fHBHENoiseFilterResultRun1)
#process.HBHENoiseFilterResultRun2Loose = cms.Path(process.fHBHENoiseFilterResultRun2Loose)
#process.HBHENoiseFilterResultRun2Tight = cms.Path(process.fHBHENoiseFilterResultRun2Tight)
#process.HBHEIsoNoiseFilterResult = cms.Path(process.fHBHEIsoNoiseFilterResult)
#########process.pprimaryVertexFilter = cms.Path(process.primaryVertexFilter )

process.load('HeavyIonsAnalysis.Configuration.hfCoincFilter_cff')
#process.phfCoincFilter1 = cms.Path(process.hfCoincFilter)
#process.phfCoincFilter2 = cms.Path(process.hfCoincFilter2)
#process.phfCoincFilter3 = cms.Path(process.hfCoincFilter3)
#process.phfCoincFilter4 = cms.Path(process.hfCoincFilter4)
#process.phfCoincFilter5 = cms.Path(process.hfCoincFilter5)

#process.pclusterCompatibilityFilter = cms.Path(process.clusterCompatibilityFilter)

#process.pAna = cms.EndPath(process.skimanalysis)

# Customization
##########################################UE##########################################
from CondCore.DBCommon.CondDBSetup_cfi import *
process.uetable = cms.ESSource("PoolDBESSource",
      DBParameters = cms.PSet(
        messageLevel = cms.untracked.int32(0)
        ),
      timetype = cms.string('runnumber'),
      toGet = cms.VPSet(
          cms.PSet(record = cms.string("JetCorrectionsRecord"),
                   tag = cms.string("UETableCompatibilityFormat_PF_v02_offline"),
                   label = cms.untracked.string("UETable_PF")
          ),
          cms.PSet(record = cms.string("JetCorrectionsRecord"),
                   tag = cms.string("UETableCompatibilityFormat_Calo_v02_offline"),
                   label = cms.untracked.string("UETable_Calo")
          )
      ), 
      connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS")
)
process.es_prefer_uetable = cms.ESPrefer('PoolDBESSource','uetable')
##########################################UE##########################################
import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltfilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltfilter.HLTPaths = ["HLT_HIUPCL1DoubleEG2NotHF2_v1","HLT_HIUPCL1SingleEG5NotHF2_v1"]
#process.hltfilter.HLTPaths = ["HLT_HIUPCL1SingleEG5NotHF2_v1"]
#process.hltfilter.HLTPaths = ["HLT_HIUPCL1DoubleEG2NotHF2_v1"]
#process.hltfilter.HLTPaths = ["HLT_HIUPC*"]
process.hltfilter.throw = False
process.hltfilter.andOr = True
process.hltfilter.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")

process.superFilterPath = cms.Path(process.hltfilter)
process.skimanalysis.superFilters = cms.vstring("superFilterPath")
process.pixel = cms.EDAnalyzer('PixelHitAnalyzer',
                             vertexSrc = cms.vstring('offlinePrimaryVertices'),
                             trackSrc = cms.untracked.InputTag('generalTracks')
)

process.ana_step = cms.Path(process.hltanalysis *
                            process.hltobject *
			    #process.hltAOD *
                            #process.hiEvtAnalyzer*
                            #process.jetSequences +
                            process.ggHiNtuplizer +
                            #process.ggHiNtuplizerGED +
                            #process.pfcandAnalyzer +
                            #process.HiForest +
                            #process.trackSequencesPbPb +
                            #process.trackSequencesPP +
                            process.pixel #+
			    #process.hltfilter
                            #process.hcalNoise #+
                            #process.tupelPatSequence
                            )

##filter all path with the production filter sequence
#for path in process.paths:
#   getattr(process,path)._seq = process.hltfilter * getattr(process,path)._seq
