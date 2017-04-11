import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        '/store/mc/RunIISummer16DR80Premix/ZToEE_NNPDF30_13TeV-powheg_M_50_120/AODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/004FD0AD-0AC7-E611-99DC-E0071B7AC760.root'
    )
)

#process.content = cms.EDAnalyzer("EventContentAnalyzer")

process.options = cms.untracked.PSet( SkipEvent = cms.untracked.vstring('ProductNotFound') )

process.ele27wptight = cms.EDAnalyzer('HLTTrkAna',
	triggers = cms.InputTag('hltTriggerSummaryAOD::HLT'),
	trigresults = cms.InputTag('TriggerResults::HLT'),
        processname = cms.string('HLT'),
	triggername = cms.string('HLT_Ele27_WPTight_Gsf_v7'),
	filternames0 = cms.vstring('hltEle27WPTightHcalIsoFilter::HLT'),
	filternames1 = cms.vstring('hltEle27WPTightPixelMatchFilter::HLT'),
	filternames2 = cms.vstring('hltEle27WPTightGsfDphiFilter::HLT'),
	genparticles = cms.InputTag('genParticles::HLT'),
	pusummary = cms.InputTag('addPileupInfo::HLT')
)

process.ele23ele12 = process.ele27wptight.clone(
	triggername = cms.string('HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v9'),
	filternames0 = cms.vstring('hltEle23Ele12CaloIdLTrackIdLIsoVLHcalIsoLeg2Filter::HLT'),
	filternames1 = cms.vstring('hltEle23Ele12CaloIdLTrackIdLIsoVLPixelMatchLeg1Filter::HLT','hltEle23Ele12CaloIdLTrackIdLIsoVLPixelMatchLeg2Filter::HLT'),
	filternames2 = cms.vstring('hltEle23Ele12CaloIdLTrackIdLIsoVLDphiLeg1Filter::HLT','hltEle23Ele12CaloIdLTrackIdLIsoVLDphiLeg2Filter::HLT'),
)

process.TFileService = cms.Service("TFileService", 
    #fileName = cms.string("histo_ZToEE.root"),
    fileName = cms.string("histo.root"),
    closeFileFast = cms.untracked.bool(True)
)

process.p = cms.Path(process.ele27wptight + process.ele23ele12)
