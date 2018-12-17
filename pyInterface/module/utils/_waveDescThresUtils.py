
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
	The yaml wave list can have the following keys:
		- 'wavelist': List of waves with the following keys per entry
			* 'name': wave name
			* 'ranges': (optional) list of multibins in which the wave is active
		- 'referenceWaves': (optional) List of reference waves with the following keys per entry
			* 'name': name of the reference wave
			* 'rank': (optional) Rank in which the wave is the reference wave. If not given, wave is reference in all ranks
			* 'ranges': (optional) list of multibins in which the wave is the reference wave
	@return [(wavename, wavedescription, isActive)] for all waves.
	isActive means the wave is active in the given multibin:
		- if the range is given, the wave is only active if the bin center lies within the range of the wave
		- if a threshold is given, the wave is only active if the mass bin center >= threshold
	'''
	waveDescriptionActive = []
	referenceWaves = None
	if os.path.splitext(waveListFileName)[1] == '.yaml':
		with open(waveListFileName, 'r') as wavelistFile:
			wavelistData = yaml.load(wavelistFile)
			if 'wavelist' not in wavelistData:
				raise ValueError("Cannot find 'wavelist' item in wave-list file '{0}'.".format(waveListFileName))
			for waveData in wavelistData['wavelist']:
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
			if 'referenceWaves' in wavelistData:
				referenceWaves = []
				for ref in wavelistData['referenceWaves']:
					useRef=True
					if 'ranges' in ref:
						useRef=False
						for refRange in ref['ranges']:
							if refRange.inBin(multibin.getBinCenters()):
								useRef=True
					if useRef:
						if not 'name' in ref:
							raise ValueError("Cannot find 'name' of reference wave in '{0}'.".format(waveListFileName))
						if 'rank' in ref:
							ref['rank'] = int(ref['rank'])
						if not ref['name'] in waveDescriptions:
							raise ValueError("Cannot find reference wave '{0}' in wave descriptions.".format(ref['name']))
						ref['description'] = waveDescriptions[ref['name']]
						referenceWaves.append(ref)

	else:
		waveDescriptionThreshold = getWaveDescThresFromWaveList(waveListFileName, waveDescriptions)
		waveDescriptionActive = [(w, d, multibin.getBinCenters()['mass'] >= t) for w,d,t in waveDescriptionThreshold]
	return waveDescriptionActive, referenceWaves
