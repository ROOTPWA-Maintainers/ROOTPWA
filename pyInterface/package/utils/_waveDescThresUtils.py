
import pyRootPwa.core
from _printingUtils import printInfo, printWarn

def getWaveDescThresFromFitResult(fitResult, keyFiles):
	waveDescThres = []
	for waveName in fitResult.waveNames():
		if waveName == "flat":
			continue

		waveDesc = waveDesc = pyRootPwa.core.waveDescription.parseKeyFile(keyFiles[waveName][0])[keyFiles[waveName][1]]

		thresholded = True
		for prodAmpIndex in xrange(fitResult.nmbProdAmps()):
			if fitResult.waveNameForProdAmp(prodAmpIndex) == waveName:
				if fitResult.prodAmp(prodAmpIndex) != 0.:
					thresholded = False

		threshold = 0.
		if thresholded:
			threshold= 1.1 * fitResult.massBinCenter()

		waveDescThres.append( (waveName, waveDesc, threshold) )
	return waveDescThres


def getWaveDescThresFromWaveList(waveListFileName, keyFiles):
	printInfo("reading amplitude names and thresholds from wave list file "
	          + "'" + waveListFileName + "'.")
	waveDescThres = []
	with open(waveListFileName, 'r') as waveListFile:
		lineNmb = 0
		for line in waveListFile:
			if line[0] == '#':  # comments start with #
				continue
			line = line.replace('\n', '')
			lineArray = line.split(" ")
			if len(lineArray) in [1, 2]:
				waveName = lineArray[0]
				if len(lineArray) == 1:
					threshold = 0
				else:
					threshold = lineArray[1]
				waveDesc = pyRootPwa.core.waveDescription.parseKeyFile(keyFiles[waveName][0])[keyFiles[waveName][1]]
				waveDescThres.append( (waveName, waveDesc, float(threshold)) )
			else:
				printWarn("cannot parse line '" + line + "' in wave list file "
				          + "'" + waveListFileName + "'.")
			lineNmb += 1
		printInfo("read " + str(lineNmb) + " lines from wave list file " + "'" + waveListFileName + "'")
	return waveDescThres


def getWaveThresFromWaveList(waveListFileName, waveDescriptions):
	printInfo("reading amplitude names and thresholds from wave list file "
	          + "'" + waveListFileName + "'.")
	waveDescThres = []
	with open(waveListFileName, 'r') as waveListFile:
		lineNmb = 0
		for line in waveListFile:
			if line[0] == '#':  # comments start with #
				continue
			line = line.replace('\n', '')
			lineArray = line.split(" ")
			if len(lineArray) in [1, 2]:
				waveName = lineArray[0]
				if len(lineArray) == 1:
					threshold = 0
				else:
					threshold = lineArray[1]
				waveDesc = waveDescriptions[waveName]
				waveDescThres.append( (waveName, waveDesc, float(threshold)) )
			else:
				printWarn("cannot parse line '" + line + "' in wave list file "
				          + "'" + waveListFileName + "'.")
			lineNmb += 1
		printInfo("read " + str(lineNmb) + " lines from wave list file " + "'" + waveListFileName + "'")
	return waveDescThres
