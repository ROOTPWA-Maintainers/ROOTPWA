import os

import pyRootPwa
import pyRootPwa.core
ROOT = pyRootPwa.ROOT

def calcIntegrals(ampFileList, maxNmbEvents=0, weightFileName=""):
	treeList = []
	waveNames = []
	integralMatrix = pyRootPwa.core.ampIntegralMatrix()
	ampFiles = []
	for ampFileName in ampFileList:
		ampFile = ROOT.TFile.Open(ampFileName, "READ")
		ampFiles.append(ampFile)
		if not ampFile:
			pyRootPwa.utils.printErr("could not open amplitude file '" + ampFileName + "'.")
		foundAmpKey = False
		for key in ampFile.GetListOfKeys():
			if not key:
				pyRootPwa.utils.printWarn("NULL pointer to TKey in file '" + ampFileName + "'.")
				continue
			keyName = key.GetName()
			keyWithoutExt, keyExt = os.path.splitext(keyName)
			if keyExt == ".amp":
				foundAmpKey = True
				tree = ampFile.Get(keyName)
				treeList.append(tree)
				meta = ampFile.Get(keyWithoutExt + ".meta")
				if not meta:
					pyRootPwa.utils.printErr("could not get metadata for waveName '" + keyWithoutExt + "'.")
					del ampFiles
					return False
				waveNames.append(meta.objectBaseName())
		if not foundAmpKey:
			pyRootPwa.utils.printWarn("no TKey in file '" + ampFileName + "'.")
	if not integralMatrix.integrate(treeList, waveNames, maxNmbEvents, weightFileName):
		pyRootPwa.utils.printErr("could not run integration")
		del ampFiles
		return None
	del ampFiles
	return integralMatrix
