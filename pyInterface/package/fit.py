import os
import random

import pyRootPwa
import pyRootPwa.core
ROOT = pyRootPwa.ROOT

def readWaveList(waveListFileName, keyFiles):
	pyRootPwa.utils.printInfo("reading amplitude names and thresholds from wave list file "
	          + "'" + waveListFileName + "'.")
	with open(waveListFileName, 'r') as waveListFile:
#	if (not waveListFile) {
#		printErr << "cannot open file '" << waveListFileName << "'. Aborting..." << endl;
#		throw;
#	}
		waveDescThres = []
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
				waveDesc = pyRootPwa.core.waveDescription.parseKeyFile(keyFiles[waveName])
				waveDescThres.append( (waveName, waveDesc, float(threshold)) )
			else:
				pyRootPwa.utils.printWarn("cannot parse line '" + line + "' in wave list file "
				          + "'" + waveListFileName + "'.")
#  			if (_debug):
#  				printDebug("reading line " + lineNmb + 1 + ": " + waveName + ", "
#  				           + "threshold = " + threshold + " MeV/c^2")
			lineNmb += 1
	pyRootPwa.utils.printInfo("read " + str(lineNmb) + " lines from wave list file " + "'" + waveListFileName + "'")
	return waveDescThres


def pwaFit(ampFileList,
           normIntegralFileName,
           accIntegralFileName,
           binningMap,
           waveListFileName,
           keyFiles,
           seed=0,
           cauchy=False,
           cauchyWidth=0.5,
           startValFileName="",
           accEventsOverride=0,
           checkHessian=False,
           saveSpace=False,
           rank=1,
           verbose=False,
           attempts=1
          ):

	waveDescThres = readWaveList(waveListFileName, keyFiles)
	massBinCenter = (binningMap['mass'][1] + binningMap['mass'][0]) / 2. # YOU CAN DO BETTER

	likelihood = pyRootPwa.core.pwaLikelihood()
	likelihood.useNormalizedAmps(True)
	if (not verbose):
		likelihood.setQuiet()
	if cauchy:
		likelihood.setPriorType(pyRootPwa.core.HALF_CAUCHY)
		likelihood.setCauchyWidth(cauchyWidth)
	if (not likelihood.init(waveDescThres,
	                        rank,
	                        massBinCenter)):
		printErr("could not initialize likelihood. Aborting...")
		return [ ]

	normIntFile = ROOT.TFile.Open(normIntegralFileName, "READ")
	if len(normIntFile.GetListOfKeys()) != 1:
		pyRootPwa.utils.printWarn("'" + normIntegralFileName + "' does not contain exactly one TKey.")
		return [ ]
	normIntMatrix = normIntFile.Get(pyRootPwa.core.ampIntegralMatrix.integralObjectName)
	if (not likelihood.addNormIntegral(normIntMatrix)):
		pyRootPwa.utils.printErr("could not add normalization integral. Aborting...")
		return [ ]
	normIntFile.Close()
	accIntFile = ROOT.TFile.Open(accIntegralFileName, "READ")
	if len(accIntFile.GetListOfKeys()) != 1:
		pyRootPwa.utils.printWarn("'" + normIntegralFileName + "' does not contain exactly one TKey.")
		return [ ]
	accIntMatrix = accIntFile.Get(pyRootPwa.core.ampIntegralMatrix.integralObjectName)
	if (not likelihood.addAccIntegral(accIntMatrix, accEventsOverride)):
		pyRootPwa.utils.printErr("could not add acceptance integral. Aborting...")
		return [ ]
	accIntFile.Close()

	for wave in waveDescThres:
		waveName = wave[0]
		ampFileName = ampFileList[waveName]
		ampFile = ROOT.TFile.Open(ampFileName, "READ")
		if not ampFile:
			pyRootPwa.utils.printErr("could not open amplitude file '" + ampFileName + "'.")
			return [ ]
		meta = pyRootPwa.core.amplitudeMetadata.readAmplitudeFile(ampFile, waveName)
		if not meta:
			pyRootPwa.utils.printErr("could not get metadata for waveName '" + waveName + "'.")
			return [ ]
		if (not likelihood.addAmplitude(meta)):
			pyRootPwa.utils.printErr("could not add amplitude '" + waveName + "'. Aborting...")
			return [ ]
	if (not likelihood.finishInit()):
		pyRootPwa.utils.printErr("could not finish initialization of likelihood. Aborting...")
		return [ ]
	lowerBound = binningMap[binningMap.keys()[0]][0]
	upperBound = binningMap[binningMap.keys()[0]][1]

	if attempts == 1:
		seeds = [ seed ]
	elif seed == 0:
		seeds = [ 0 for _ in xrange(attempts) ]
	else:
		random.seed(seed)
		seeds = [ ]
		for _ in xrange(attempts):
			while True:
				randVal = random.randint(1000, 2**32-1)
				if randVal not in seeds:
					break
			seeds.append(randVal)

	fitResults = [ ]
	for fitSeed in seeds:
		fitResult = pyRootPwa.core.pwaFit(likelihood       = likelihood,
		                                  massBinMin       = lowerBound,
		                                  massBinMax       = upperBound,
		                                  seed             = fitSeed,
		                                  startValFileName = startValFileName,
		                                  checkHessian     = checkHessian,
		                                  saveSpace        = saveSpace,
		                                  verbose          = verbose)
		fitResults.append(fitResult)
	return fitResults


def pwaNloptFit(ampFileList,
                normIntegralFileName,
                accIntegralFileName,
                binningMap,
                waveListFileName,
                keyFiles,
                seed=0,
                cauchy=False,
                cauchyWidth=0.5,
                startValFileName="",
                accEventsOverride=0,
                checkHessian=False,
                saveSpace=False,
                rank=1,
                verbose=False,
                attempts=1
               ):

	waveDescThres = readWaveList(waveListFileName, keyFiles)
	massBinCenter = (binningMap['mass'][1] + binningMap['mass'][0]) / 2. # YOU CAN DO BETTER

	likelihood = pyRootPwa.core.pwaLikelihood()
	likelihood.useNormalizedAmps(True)
	if (not verbose):
		likelihood.setQuiet()
	if cauchy:
		likelihood.setPriorType(pyRootPwa.core.HALF_CAUCHY)
		likelihood.setCauchyWidth(cauchyWidth)
	if (not likelihood.init(waveDescThres,
	                        rank,
	                        massBinCenter)):
		printErr("could not initialize likelihood. Aborting...")
		return [ ]

	normIntFile = ROOT.TFile.Open(normIntegralFileName, "READ")
	if len(normIntFile.GetListOfKeys()) != 1:
		pyRootPwa.utils.printWarn("'" + normIntegralFileName + "' does not contain exactly one TKey.")
		return [ ]
	normIntMatrix = normIntFile.Get(pyRootPwa.core.ampIntegralMatrix.integralObjectName)
	if (not likelihood.addNormIntegral(normIntMatrix)):
		pyRootPwa.utils.printErr("could not add normalization integral. Aborting...")
		return [ ]
	normIntFile.Close()
	accIntFile = ROOT.TFile.Open(accIntegralFileName, "READ")
	if len(accIntFile.GetListOfKeys()) != 1:
		pyRootPwa.utils.printWarn("'" + normIntegralFileName + "' does not contain exactly one TKey.")
		return [ ]
	accIntMatrix = accIntFile.Get(pyRootPwa.core.ampIntegralMatrix.integralObjectName)
	if (not likelihood.addAccIntegral(accIntMatrix, accEventsOverride)):
		pyRootPwa.utils.printErr("could not add acceptance integral. Aborting...")
		return [ ]
	accIntFile.Close()

	for wave in waveDescThres:
		waveName = wave[0]
		ampFileName = ampFileList[waveName]
		ampFile = ROOT.TFile.Open(ampFileName, "READ")
		if not ampFile:
			pyRootPwa.utils.printErr("could not open amplitude file '" + ampFileName + "'.")
			return [ ]
		meta = pyRootPwa.core.amplitudeMetadata.readAmplitudeFile(ampFile, waveName)
		if not meta:
			pyRootPwa.utils.printErr("could not get metadata for waveName '" + waveName + "'.")
			return [ ]
		if (not likelihood.addAmplitude(meta)):
			pyRootPwa.utils.printErr("could not add amplitude '" + waveName + "'. Aborting...")
			return [ ]
	if (not likelihood.finishInit()):
		pyRootPwa.utils.printErr("could not finish initialization of likelihood. Aborting...")
		return [ ]
	lowerBound = binningMap[binningMap.keys()[0]][0]
	upperBound = binningMap[binningMap.keys()[0]][1]

	if attempts == 1:
		seeds = [ seed ]
	elif seed == 0:
		seeds = [ 0 for _ in xrange(attempts) ]
	else:
		random.seed(seed)
		seeds = [ ]
		for _ in xrange(attempts):
			while True:
				randVal = random.randint(1000, 2**32-1)
				if randVal not in seeds:
					break
			seeds.append(randVal)

	fitResults = [ ]
	for fitSeed in seeds:
		fitResult = pyRootPwa.core.pwaNloptFit(likelihood       = likelihood,
		                                       massBinMin       = lowerBound,
		                                       massBinMax       = upperBound,
		                                       seed             = fitSeed,
		                                       startValFileName = startValFileName,
		                                       checkHessian     = checkHessian,
		                                       saveSpace        = saveSpace,
		                                       verbose          = verbose)
		fitResults.append(fitResult)
	return fitResults
