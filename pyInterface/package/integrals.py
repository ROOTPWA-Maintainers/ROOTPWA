import os

import pyRootPwa
import pyRootPwa.core
ROOT = pyRootPwa.ROOT

def calcIntegrals(ampFileList, maxNmbEvents=0, weightFileName=""):
	integralMatrix = pyRootPwa.core.ampIntegralMatrix()
	ampFiles = []
	ampMetas = []
	for waveName in ampFileList:
		ampFileName = ampFileList[waveName]
		ampFile = ROOT.TFile.Open(ampFileName, "READ")
		ampFiles.append(ampFile)
		if not ampFile:
			pyRootPwa.utils.printErr("could not open amplitude file '" + ampFileName + "'.")
		ampMeta = pyRootPwa.core.amplitudeMetadata.readAmplitudeFile(ampFile, waveName)
		ampMetas.append(ampMeta)
		if not ampMeta:
			pyRootPwa.utils.printErr("could not read amplitude file '" + ampFileName + "'.")
	if not integralMatrix.integrate(ampMetas, maxNmbEvents, weightFileName):
		pyRootPwa.utils.printErr("could not run integration")
		del ampFiles
		return None
	del ampFiles
	return integralMatrix
