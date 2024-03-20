import FWCore.ParameterSet.Config as cms

process = cms.Process("KSHORTS")

# Use the tracks_and_vertices.root file as input.
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring("file://run321167_ZeroBias_AOD.root"))
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1)) #input = cms.untracked.int32(500))

# Suppress messages that are less important than ERRORs.
process.MessageLogger = cms.Service("MessageLogger",
    destinations = cms.untracked.vstring("cout"),
    cout = cms.untracked.PSet(threshold = cms.untracked.string("ERROR")))

# Load part of the CMSSW reconstruction sequence to make vertexing possible.
# We'll need the CMS geometry and magnetic field to follow the true, non-helical
# shapes of tracks through the detector.
process.load("Configuration/StandardSequences/FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag =  GlobalTag(process.GlobalTag, "auto:run2_data")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")

# Copy most of the vertex producer's parameters, but accept tracks with
# progressively more strict quality.
process.load("RecoVertex.V0Producer.generalV0Candidates_cfi")

# loose
process.SecondaryVerticesFromLooseTracks = process.generalV0Candidates.clone(
    trackRecoAlgorithm = cms.InputTag("generalTracks"),
    doKshorts = cms.bool(True),
    doLambdas = cms.bool(True),
    trackQualities = cms.string("loose"),
    innerHitPosCut = cms.double(-1.),
    cosThetaXYCut = cms.double(-1.),
    )

# tight
process.SecondaryVerticesFromTightTracks = process.SecondaryVerticesFromLooseTracks.clone(
    trackQualities = cms.string("tight"),
    )

# highPurity
process.SecondaryVerticesFromHighPurityTracks = process.SecondaryVerticesFromLooseTracks.clone(
    trackQualities = cms.string("highPurity"),
    )

# Run all three versions of the algorithm.
process.path = cms.Path(process.SecondaryVerticesFromLooseTracks *
               		 process.SecondaryVerticesFromTightTracks *
               		 process.SecondaryVerticesFromHighPurityTracks)

# Writer to a new file called output.root.  Save only the new K-shorts and the
# primary vertices (for later exercises).
process.output = cms.OutputModule(
    "PoolOutputModule",
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring("path")),
    outputCommands = cms.untracked.vstring(
   	 "drop *",
   	 "keep *_*_*_KSHORTS",
   	 "keep *_offlineBeamSpot_*_*",
   	 "keep *_offlinePrimaryVertices_*_*",
   	 "keep *_offlinePrimaryVerticesWithBS_*_*",
   	 ),
    fileName = cms.untracked.string("output.root")
    )
process.endpath = cms.EndPath(process.output)
