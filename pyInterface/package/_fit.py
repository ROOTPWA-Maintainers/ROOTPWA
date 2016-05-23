import random

import pyRootPwa.core
import pyRootPwa.utils
ROOT = pyRootPwa.utils.ROOT

def pwaFit(eventAndAmpFileDict,
           normIntegralFileName,
           accIntegralFileName,
           multiBin,
           waveListFileName,
           waveDescriptions,
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

	waveDescThres = pyRootPwa.utils.getWaveThresFromWaveList(waveListFileName, waveDescriptions)
	massBinCenter = (multiBin.boundaries['mass'][1] + multiBin.boundaries['mass'][0]) / 2. # YOU CAN DO BETTER

	likelihood = pyRootPwa.initLikelihood(waveDescThres = waveDescThres,
	                                      massBinCenter = massBinCenter,
	                                      eventAndAmpFileDict = eventAndAmpFileDict,
	                                      normIntegralFileName = normIntegralFileName,
	                                      accIntegralFileName = accIntegralFileName,
	                                      multiBin = multiBin,
	                                      accEventsOverride = accEventsOverride,
	                                      cauchy = cauchy,
	                                      cauchyWidth = cauchyWidth,
	                                      rank = rank,
	                                      verbose = verbose)
	if not likelihood:
		pyRootPwa.utils.printErr("error while initializing likelihood. Aborting...")
		return [ ]

	lowerBound = multiBin.boundaries['mass'][0]
	upperBound = multiBin.boundaries['mass'][1]

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

	waveDescThres = pyRootPwa.utils.getWaveDescThresFromWaveList(waveListFileName, keyFiles)
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
