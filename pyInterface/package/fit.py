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
				waveDesc = pyRootPwa.core.waveDescription.parseKeyFile(keyFiles[waveName][0])[keyFiles[waveName][1]]
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

	likelihood = pyRootPwa.initLikelihood(waveDescThres = waveDescThres,
	                                      massBinCenter = massBinCenter,
	                                      ampFileList = ampFileList,
	                                      normIntegralFileName = normIntegralFileName,
	                                      accIntegralFileName = accIntegralFileName,
	                                      accEventsOverride = accEventsOverride,
	                                      cauchy = cauchy,
	                                      cauchyWidth = cauchyWidth,
	                                      rank = rank,
	                                      verbose = verbose)
	if not likelihood:
		pyRootPwa.utils.printErr("error while initializing likelihood. Aborting...")
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

	likelihood = pyRootPwa.initLikelihood(waveDescThres = waveDescThres,
	                                      massBinCenter = massBinCenter,
	                                      ampFileList = ampFileList,
	                                      normIntegralFileName = normIntegralFileName,
	                                      accIntegralFileName = accIntegralFileName,
	                                      accEventsOverride = accEventsOverride,
	                                      cauchy = cauchy,
	                                      cauchyWidth = cauchyWidth,
	                                      rank = rank,
	                                      verbose = verbose)
	if not likelihood:
		pyRootPwa.utils.printErr("error while initializing likelihood. Aborting...")
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
