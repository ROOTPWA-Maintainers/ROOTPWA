import os
import numpy
import pyRootPwa
import pyRootPwa.core
ROOT = pyRootPwa.ROOT

def calcIntegralsOnTheFly(eventFileName, keyFileList, integralFile, binningMap = {} ,maxNmbEvents = -1, startEvent = 0):

	amplitudes      = []
	keyFileContents = []
	hashers         = []
	waveNames       = []
	metadataObject =  pyRootPwa.core.ampIntegralMatrixMetadata()

	outFile = pyRootPwa.ROOT.TFile.Open(integralFile, "CREATE")
	if not outFile: # Do this up here. Without the output file, nothing else makes sense
		pyRootPwa.utils.printErr("could not open outFile. Abort...")
		return False
		
	for keyFile in keyFileList:
		waveDescription = pyRootPwa.core.waveDescription.parseKeyFile(keyFile[0])[keyFile[1]]
		if not metadataObject.addKeyFileContent(waveDescription.keyFileContent()):
			pyRootPwa.utils.printWarn("found same keyFileContent twice")
		(result, amplitude) = waveDescription.constructAmplitude()
		if not result:
			pyRootPwa.utils.printErr('could not construct amplitude for keyfile "' + keyFile[0] + '" (ID '+str(keyFile[1])+'). Abort...')
			return False
		amplitude.init()
		amplitudes.append(amplitude)
		hashers.append(pyRootPwa.core.hashCalculator())
		waveNames.append(pyRootPwa.core.waveDescription.waveNameFromTopology(amplitude.decayTopology()))

	eventFile = pyRootPwa.ROOT.TFile.Open(eventFileName, "READ")


	if not eventFile:
		pyRootPwa.utils.printErr("could not open eventFile. Abort...")
		return False
	eventMeta  = pyRootPwa.core.eventMetadata.readEventFile(eventFile)

	if not metadataObject.setBinningMap(eventMeta.binningMap()):
		pyRootPwa.utils.printErr("could not setBinningMap. Abort...")
		return False
	prodNames  = eventMeta.productionKinematicsParticleNames()
	decayNames = eventMeta.decayKinematicsParticleNames()
	for amplitude in amplitudes:
		topo = amplitude.decayTopology()
		if not topo.initKinematicsData(prodNames, decayNames):
			pyRootPwa.utils.printErr("could not initKinematicsData(). Abort...")
			return False

	eventTree = eventMeta.eventTree()
	nEvents   = eventTree.GetEntries()
	minEvent = startEvent
	maxEvent = nEvents
	if maxNmbEvents	> -1:
		maxEvent = min(maxEvent, startEvent + maxNmbEvents)
	if not metadataObject.addEventMetadata(eventMeta, minEvent, maxEvent):
		pyRootPwa.utils.printErr("could not addEventMetadata. Abort...")
		return False
	if not metadataObject.setBinningMap(binningMap):
		pyRootPwa.utils.printErr("could not setBinningMap. Abort...")
		return False

	prodKinMomenta  = pyRootPwa.ROOT.TClonesArray("TVector3")
	decayKinMomenta = pyRootPwa.ROOT.TClonesArray("TVector3")
	eventTree.SetBranchAddress(pyRootPwa.core.eventMetadata.productionKinematicsMomentaBranchName, prodKinMomenta)
	eventTree.SetBranchAddress(pyRootPwa.core.eventMetadata.decayKinematicsMomentaBranchName, decayKinMomenta)
	ampWaveNameMap = {}
	integralMatrix = pyRootPwa.core.ampIntegralMatrix()
	integralMatrix.setWaveNames(waveNames)

	binningVariables = {}
	for key in binningMap:
		binningVariables[key] = numpy.array(1, dtype = float)
		eventTree.SetBranchAddress(key, binningVariables[key])
	
	for waveName in waveNames:
		ampWaveNameMap[waveName] = 0.+0.j
	pyRootPwa.utils.printInfo("starting event loop")
	skippedEvents = 0
	for evt in range(minEvent, maxEvent):
		if evt%100 == 0:
			pyRootPwa.utils.printInfo("event index = " + str(evt))
		eventTree.GetEvent(evt)
		for key in binningMap:
			if binningVariables[key] < binningMap[key][0] or  binningVariables[key] >= binningMap[key][1]:
				skippedEvents += 1
				continue
		for a, amplitude in enumerate(amplitudes):
			topo = amplitude.decayTopology()
			if not topo.readKinematicsData(prodKinMomenta, decayKinMomenta):
				pyRootPwa.utils.printErr("could not loadKinemaricsData. Abort...")
				return False
			ampl = amplitude()
			hashers[a].Update(ampl)
			ampWaveNameMap[waveNames[a]] = ampl
		if not integralMatrix.addEvent(ampWaveNameMap):
			pyRootPwa.utils.printErr("could not addEvent. Abort...")
			return False
	if not metadataObject.setAmpIntegralMatrix(integralMatrix):
		pyRootPwa.utils.printErr("could not setAmpIntegralMatrix. Abort...")
		return False
	for hasher in hashers:
		if not metadataObject.addAmplitudeHash(hasher.hash()):
			pyRootPwa.utils.printWarn("could not addAmplitudehash")
	if not  metadataObject.writeToFile(outFile):
		pyRootPwa.utils.printErr("could not writeToFile. Abort...")
		return False

	outFile.Close()
	eventFile.Close()
	pyRootPwa.utils.printInfo(str(skippedEvents) + "events rejected due the binning")




