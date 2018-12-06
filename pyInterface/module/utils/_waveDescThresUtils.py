
import os
import yaml
from _printingUtils import printInfo, printWarn

def getWaveDescThresFromFitResult(fitResult, waveDescriptions):
	waveDescThres = []
	for waveName in fitResult.waveNames():
		if waveName == "flat":
			continue

		waveDesc = waveDescriptions[waveName]

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


def getWaveDescThresFromWaveList(waveListFileName, waveDescriptions):
	printInfo("reading amplitude names and thresholds from wave list file "
	          + "'" + waveListFileName + "'.")
	waveDescThres = []
	with open(waveListFileName, 'r') as waveListFile:
		lineNmb = 0
		for line in waveListFile:
			line = line.strip()
			if line[0] == '#':  # comments start with #
				continue
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


def getWaveDescriptionActiveFromWavelist(waveListFileName, waveDescriptions, multibin):
	'''
	@return [(wavename, wavedescription, isActive)] for all waves.
	isActive means the wave is active in the given multibin:
		- if the range is given, the wave is only active if the bin center lies within the range of the wave
		- if a threshold is given, the wave is only active if the mass bin center >= threshold
	'''
	waveDescriptionActive = []
	if os.path.splitext(waveListFileName)[1] == '.yaml':
		with open(waveListFileName, 'r') as wavelistFile:
			wavelist = yaml.load(wavelistFile)
			for waveData in wavelist:
				isActive = True
				if 'ranges' in waveData:
					isActive = False
					for waveRange in waveData['ranges']:
						if waveRange.inBin(multibin.getBinCenters()):
							isActive = True
							break
				if 'threshold' in waveData and multibin.getBinCenters()['mass'] < waveData['threshold']:
					isActive = False
				waveDescriptionActive.append((waveData['name'], waveDescriptions[waveData['name']], isActive))
	else:
		waveDescriptionThreshold = getWaveDescThresFromWaveList(waveListFileName, waveDescriptions)
		waveDescriptionActive = [(w, d, multibin.getBinCenters()['mass'] >= t) for w,d,t in waveDescriptionThreshold]
	return waveDescriptionActive
