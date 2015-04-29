import os

import pyRootPwa
import pyRootPwa.core
ROOT = pyRootPwa.ROOT

def readWaveList(waveListFileName):
	pyRootPwa.utils.printInfo("reading amplitude names and thresholds from wave list file "
	          + "'" + waveListFileName + "'.")
	with open(waveListFileName, 'r') as waveListFile:
# 	if (not waveListFile) {
# 		printErr << "cannot open file '" << waveListFileName << "'. aborting." << endl;
# 		throw;
# 	}
		waveNames = []
		waveThresholds = []
		lineNmb = 0
		for line in waveListFile:
			if (line[0] == '#'):  # comments start with #
				continue
			line = line.replace('\n', '')
			lineArray = line.split(" ")
			if(len(lineArray) >= 1 and len(lineArray) <= 2):
				waveName = lineArray[0]
				if(len(lineArray) == 1):
					threshold = 0
				else:
					threshold = lineArray[1]
				waveNames.append(waveName)
				waveThresholds.append(float(threshold))
			else:
				pyRootPwa.utils.printWarn("cannot parse line '" + line + "' in wave list file "
				          + "'" + waveListFileName + "'.")
#  			if (_debug):
#  				printDebug("reading line " + lineNmb + 1 + ": " + waveName + ", "
#  				           + "threshold = " + threshold + " MeV/c^2")
			lineNmb += 1
	pyRootPwa.utils.printInfo("read " + str(lineNmb) + " lines from wave list file " + "'" + waveListFileName + "'")
	return (waveNames, waveThresholds)


def pwaFit(ampFileList, normIntegralFileName, accIntegralFileName, binningMap, waveListFileName, seed=0, maxNmbEvents=0, startValFileName="", accEventsOverride=0, runHesse=False, rank=1, verbose=False):
	treeDict = {}
	waveNames = []
	ampFiles = []

	normIntFile = ROOT.TFile.Open(normIntegralFileName, "READ")
	if len(normIntFile.GetListOfKeys()) != 1:
		pyRootPwa.utils.printWarn("'" + normIntegralFileName + "' does not contain exactly one TKey.")
		return False
	normIntMatrix = normIntFile.Get(normIntFile.GetListOfKeys()[0].GetName())
	accIntFile = ROOT.TFile.Open(accIntegralFileName, "READ")
	if len(accIntFile.GetListOfKeys()) != 1:
		pyRootPwa.utils.printWarn("'" + normIntegralFileName + "' does not contain exactly one TKey.")
		return False
	accIntMatrix = accIntFile.Get(normIntFile.GetListOfKeys()[0].GetName())

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
				meta = ampFile.Get(keyWithoutExt + ".meta")
				if not meta:
					pyRootPwa.utils.printErr("could not get metadata for waveName '" + keyWithoutExt + "'.")
					del ampFiles
					return False
				waveNames.append(meta.objectBaseName())
				treeDict[meta.objectBaseName()] = tree
		if not foundAmpKey:
			pyRootPwa.utils.printWarn("no TKey in file '" + ampFileName + "'.")
	(waveNames, waveThresholds) = readWaveList(waveListFileName)
	lowerBound = binningMap[binningMap.keys()[0]][0]
	upperBound = binningMap[binningMap.keys()[0]][1]
	fitResult = pyRootPwa.core.pwaFit(
	                                  ampTreesDict = treeDict,
	                                  normMatrix = normIntMatrix,
	                                  accMatrix = accIntMatrix,
	                                  waveNames = waveNames,
	                                  waveThresholds = waveThresholds,
	                                  massBinMin = lowerBound,
	                                  massBinMax = upperBound,
	                                  seed = seed,
	                                  startValFileName = startValFileName,
	                                  accEventsOverride = accEventsOverride,
	                                  rank = rank,
	                                  runHesse = runHesse,
	                                  verbose = verbose
	                                  )
	del ampFiles
	return fitResult
