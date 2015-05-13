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


def pwaFit(likelihood, seed=0, massBinLower=0, massBinUpper=0, startValFileName="", checkHessian=False, verbose=False):
	fitResult = pyRootPwa.core.pwaFit(
	                                  likelihood = likelihood,
	                                  seed = seed,
	                                  massBinLower = massBinLower,
	                                  massBinUpper = massBinUpper,
	                                  startValFileName = startValFileName,
	                                  checkHessian = checkHessian,
	                                  verbose = verbose
	                                  )
	return fitResult


def pwaNloptFit(likelihood, seed=0, cauchy=False, massBinLower=0, massBinUpper=0, startValFileName="", checkHessian=False, saveSpace=False, verbose=False):
	fitResult = pyRootPwa.core.pwaNloptFit(
	                                       likelihood = likelihood,
	                                       seed = seed,
	                                       cauchy = cauchy,
	                                       massBinLower = massBinLower,
	                                       massBinUpper = massBinUpper,
	                                       startValFileName = startValFileName,
	                                       checkHessian = checkHessian,
	                                       saveSpace = saveSpace,
	                                       verbose = verbose
	                                       )
	return fitResult
