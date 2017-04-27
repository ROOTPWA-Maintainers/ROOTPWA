import pyRootPwa.core
import pyRootPwa.utils
ROOT = pyRootPwa.utils.ROOT

def calcIntegrals(integralFileName, eventAndAmpFileDict, multiBin, weightFileName=""):
	outputFile = pyRootPwa.ROOT.TFile.Open(integralFileName, "NEW")
	if not outputFile:
		pyRootPwa.utils.printWarn("cannot open output file '" + integralFileName + "'. Aborting...")
		return False
	integralMetaData = pyRootPwa.core.ampIntegralMatrixMetadata()
	integralMetaData.setBinningMap(multiBin.boundaries)
	integrals = []
	for eventFileName, ampFileDict in eventAndAmpFileDict.iteritems():
		integrals.append(pyRootPwa.core.ampIntegralMatrix())
		ampMetas = []
		(eventFile, eventMeta) = pyRootPwa.utils.openEventFile(eventFileName)
		if not eventFile or not eventMeta:
			pyRootPwa.utils.printErr("could not open event file '" + eventFileName + "'. Aborting...")
			return False
		if not integralMetaData.addEventMetadata(eventMeta):
			pyRootPwa.utils.printErr("could not add event metadata from event file '" + eventFileName + "' to integral metadata. Aborting...")
			return False
		for waveName, ampFileName in ampFileDict.iteritems():
			ampFile = ROOT.TFile.Open(ampFileName, "READ")
			if not ampFile:
				pyRootPwa.utils.printErr("could not open amplitude file '" + ampFileName + "'. Aborting...")
				return False
			ampMeta = pyRootPwa.core.amplitudeMetadata.readAmplitudeFile(ampFile, waveName)
			if not ampMeta:
				pyRootPwa.utils.printErr("could not read metadata from amplitude file '" + ampFileName + "'. Aborting...")
				return False
			ampMetas.append(ampMeta)

			if len(integrals) == 1: # first eventFieldID of the current bin
				if not integralMetaData.addKeyFileContent(ampMeta.keyfileContent()):
					return False
			else:
				if not integralMetaData.hasKeyFileContent(ampMeta.keyfileContent()):
					pyRootPwa.utils.printErr("keyfile content of additional eventFiledID missing in first eventFieldID.")
					return False
					
			if not integralMetaData.addAmplitudeHash(ampMeta.contentHash()):
				# This error is not fatal, since in special cases the same hash can appear twice:
				# e.g. in freed-isobar analyses with spin zero, the angular dependences are constant
				# and the shape is either 0 or 1. If two such waves accidentally have the same number
				# of events, both will also have the same hash.
				pyRootPwa.utils.printWarn("could not add the amplitude hash.")
		if not integrals[-1].integrate(ampMetas, -1, weightFileName, eventMeta, multiBin.boundaries):
			pyRootPwa.utils.printErr("could not run integration. Aborting...")
			return False
	integralMatrix = integrals[0]
	if len(integrals) > 1:
		for integral in integrals[1:]:
			integralMatrix += integral
	if not integralMetaData.setAmpIntegralMatrix(integralMatrix):
		pyRootPwa.utils.printErr("could not add the integral matrix to the metadata object. Aborting...")
		return False
	if not integralMetaData.writeToFile(outputFile):
		pyRootPwa.utils.printErr("could not write integral objects to file. Aborting...")
		return False
	outputFile.Close()
	return True
