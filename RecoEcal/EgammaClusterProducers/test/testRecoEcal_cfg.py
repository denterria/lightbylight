import FWCore.ParameterSet.Config as cms

process = cms.Process('testRecoEcal')
process.load('RecoEcal.Configuration.RecoEcal_cff')
process.load("Configuration.StandardSequences.GeometryDB_cff")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)

process.source = cms.Source("PoolSource",
#    debugVerbosity = cms.untracked.uint32(1),
#    debugFlag = cms.untracked.bool(True),
    fileNames = cms.untracked.vstring(
         'file:/afs/cern.ch/work/r/rchudasa/private/CMSSW_7_5_8_patch3/src/Configuration/Generator/python/exclusive_gg/Noptcut_RECO_check_100events.root'
    )
)

process.out = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring(
        'drop *', 
        'keep *_*_*_testRecoEcal'),
    fileName = cms.untracked.string('output_testRecoEcal.root')
)

process.p = cms.Path(process.ecalClusters)
process.outpath = cms.EndPath(process.out)

