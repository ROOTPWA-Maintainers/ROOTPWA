
import pyRootPwa.core
import pyRootPwa.utils

def getWaveDescThresFromWaveList(waveListFileName, keyFiles):
	pyRootPwa.utils.printInfo("reading amplitude names and thresholds from wave list file "
	          + "'" + waveListFileName + "'.")
	waveDescThres = []
	with open(waveListFileName, 'r') as waveListFile:
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
			lineNmb += 1
		pyRootPwa.utils.printInfo("read " + str(lineNmb) + " lines from wave list file " + "'" + waveListFileName + "'")
	return waveDescThres
