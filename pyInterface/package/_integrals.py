import pyRootPwa.core
import pyRootPwa.utils
ROOT = pyRootPwa.utils.ROOT

def calcIntegrals(integralFileName, ampFileDict, maxNmbEvents=0, weightFileName=""):
	integralMatrix = pyRootPwa.core.ampIntegralMatrix()
	ampFiles = []
	ampMetas = []
	outputFile = pyRootPwa.ROOT.TFile.Open(integralFileName, "NEW")
	if not outputFile:
		pyRootPwa.utils.printWarn("cannot open output file '" + integralFileName + "'. Aborting...")
		return False
	for waveName in ampFileDict:
		ampFileName = ampFileDict[waveName]
		ampFile = ROOT.TFile.Open(ampFileName, "READ")
		ampFiles.append(ampFile)
		if not ampFile:
			pyRootPwa.utils.printErr("could not open amplitude file '" + ampFileName + "'. Aborting...")
			return False
		ampMeta = pyRootPwa.core.amplitudeMetadata.readAmplitudeFile(ampFile, waveName)
		ampMetas.append(ampMeta)
		if not ampMeta:
			pyRootPwa.utils.printErr("could not read amplitude file '" + ampFileName + "'.")
	if not integralMatrix.integrate(ampMetas, maxNmbEvents, weightFileName):
		pyRootPwa.utils.printErr("could not run integration. Aborting...")
		del ampFiles
		return False
	integralMetaData = pyRootPwa.core.ampIntegralMatrixMetadata()
	first = True
	for ampMeta in ampMetas:
		integralMetaData.addKeyFileContent(ampMeta.keyfileContent())
		eventMetas = ampMeta.eventMetadata()
		if len(eventMetas) == 0:
			pyRootPwa.utils.printErr("no event metadata found. Aborting...")
			return False
		if len(eventMetas) > 1:
			pyRootPwa.utils.printErr("more than one event file per amplitude file not implemented at the moment. Aborting...")
			return False
		if first:
			eventMeta = eventMetas[0]
			binningMap  = eventMeta.binningMap()
			nmbEvents = ampMeta.amplitudeTree().GetEntries()
			startEvent = 0 # No start Event given at the moment, but you never know what the future brings
			if maxNmbEvents == 0:
				stopEvent = nmbEvents - 1
			else:
				stopEvent  = min(maxNmbEvents, nmbEvents) - 1
			if not 	integralMetaData.addEventMetadata(eventMeta, startEvent, stopEvent):
				pyRootPwa.utils.printErr("could not add event metadata to integral metadata. Aborting...")
				return False
			integralMetaData.setBinningMap(binningMap)
			first = False
		else:
			if not eventMetas[0] == eventMeta:
				pyRootPwa.utils.printErr("amplitude files with non-matching event metadatas. Aborting...")
				return False
		if not integralMetaData.addAmplitudeHash(ampMeta.contentHash()):
			pyRootPwa.utils.printWarn("could not add the amplitude hash.")
			# This error is not fatal, since in special cases the same hash can appear twice:
			# e.g. in freed-isobar analyses with spin zero, the angular dependences are constant
			# and the shape is either 0 or 1. If two such waves accidentally have the same number
			# of events, both will also have the same hash.
	if not integralMetaData.setAmpIntegralMatrix(integralMatrix):
		pyRootPwa.utils.printErr("could not add the integral matrix to the metadata object. Aborting...")
		return False
	del ampFiles
	if not integralMetaData.writeToFile(outputFile):
		pyRootPwa.utils.printErr("could not write integral objects to file. Aborting...")
		return False
	outputFile.Close()
	return True
