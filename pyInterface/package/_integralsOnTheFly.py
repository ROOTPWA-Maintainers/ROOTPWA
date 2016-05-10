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
		pyRootPwa.utils.printErr("could not open output file. Aborting...")
		return False

	for keyFile in keyFileList:
		waveDescription = pyRootPwa.core.waveDescription.parseKeyFile(keyFile[0])[keyFile[1]]
		if not metadataObject.addKeyFileContent(waveDescription.keyFileContent()):
			pyRootPwa.utils.printWarn("could not add keyfile content.")
		(result, amplitude) = waveDescription.constructAmplitude()
		if not result:
			pyRootPwa.utils.printErr('could not construct amplitude for keyfile "' + keyFile[0] + '" (ID '+str(keyFile[1])+'). Aborting...')
			return False
		amplitude.init()
		amplitudes.append(amplitude)
		hashers.append(pyRootPwa.core.hashCalculator())
		waveNames.append(pyRootPwa.core.waveDescription.waveNameFromTopology(amplitude.decayTopology()))

	eventFile = pyRootPwa.ROOT.TFile.Open(eventFileName, "READ")

	if not eventFile:
		pyRootPwa.utils.printErr("could not open event file. Aborting...")
		return False
	eventMeta  = pyRootPwa.core.eventMetadata.readEventFile(eventFile)

	metadataObject.setBinningMap(eventMeta.binningMap())
	prodNames  = eventMeta.productionKinematicsParticleNames()
	decayNames = eventMeta.decayKinematicsParticleNames()
	for amplitude in amplitudes:
		topo = amplitude.decayTopology()
		if not topo.initKinematicsData(prodNames, decayNames):
			pyRootPwa.utils.printErr("could not initialize the decay topology with the kinematics data. Aborting...")
			return False

	eventTree = eventMeta.eventTree()
	nEvents   = eventTree.GetEntries()
	minEvent = startEvent
	maxEvent = nEvents
	if maxNmbEvents	> -1:
		maxEvent = min(maxEvent, startEvent + maxNmbEvents)
	if not metadataObject.addEventMetadata(eventMeta, minEvent, maxEvent):
		pyRootPwa.utils.printErr("could not add event metadata to integral metadata. Aborting...")
		return False
	metadataObject.setBinningMap(binningMap)

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
	pyRootPwa.utils.printInfo("starting event loop.")
	skippedEvents = 0
	progressBar = pyRootPwa.utils.progressBar(minEvent, maxEvent)
	progressBar.start()
	for evt_i in range(minEvent, maxEvent):
		progressBar.update(evt_i)
		eventTree.GetEvent(evt_i)
		for key in binningMap:
			if binningVariables[key] < binningMap[key][0] or  binningVariables[key] >= binningMap[key][1]:
				skippedEvents += 1
				continue
		for amp_i, amplitude in enumerate(amplitudes):
			topo = amplitude.decayTopology()
			if not topo.readKinematicsData(prodKinMomenta, decayKinMomenta):
				pyRootPwa.utils.printErr("could not load kinematics data. Aborting...")
				return False
			ampl = amplitude()
			hashers[amp_i].Update(ampl)
			ampWaveNameMap[waveNames[amp_i]] = ampl
		if not integralMatrix.addEvent(ampWaveNameMap):
			pyRootPwa.utils.printErr("could not add event to integral matrix. Aborting...")
			return False
	if not metadataObject.setAmpIntegralMatrix(integralMatrix):
		pyRootPwa.utils.printErr("could not add the integral matrix to the metadata object. Aborting...")
		return False
	for hasher in hashers:
		if not metadataObject.addAmplitudeHash(hasher.hash()):
			# !!! why is this not a fatal problem?
			pyRootPwa.utils.printWarn("could not add the amplitude hash.")
	if not metadataObject.writeToFile(outFile):
		pyRootPwa.utils.printErr("could not write integral objects to file. Aborting...")
		return False

	outFile.Close()
	eventFile.Close()
	pyRootPwa.utils.printInfo(str(skippedEvents) + "events rejected because they are outside the binning.")
