import os

import pyRootPwa
import pyRootPwa.core
ROOT = pyRootPwa.ROOT

def calcIntegrals(ampFileList, maxNmbEvents=0, weightFileName=""):
	treeList = []
	waveNames = []
	integralMatrix = pyRootPwa.core.ampIntegralMatrix()
	ampFiles = []
	for waveName in ampFileList:
		ampFileName = ampFileList[waveName]
		ampFile = ROOT.TFile.Open(ampFileName, "READ")
		ampFiles.append(ampFile)
		if not ampFile:
			pyRootPwa.utils.printErr("could not open amplitude file '" + ampFileName + "'.")
		meta = ampFile.Get(waveName + ".meta")
		if not meta:
			pyRootPwa.utils.printErr("could not get metadata for waveName '" + waveName + "'.")
			del ampFiles
			return False
		waveNames.append(meta.objectBaseName())
		tree = ampFile.Get(waveName + ".amp")
		if not tree:
			pyRootPwa.utils.printErr("could not get amplitude tree for waveName '" + waveName + "'.")
			del ampFiles
			return False
		treeList.append(tree)
	if not integralMatrix.integrate(treeList, waveNames, maxNmbEvents, weightFileName):
		pyRootPwa.utils.printErr("could not run integration")
		del ampFiles
		return None
	del ampFiles
	return integralMatrix
