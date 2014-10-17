#!/usr/bin/env python

import pprint
import random
import sys

import pyRootPwa
import pyRootPwa.core

if __name__ == "__main__":

	# making the output nice for multi-threading
	pyRootPwa.utils.stdoutisatty = sys.stdout.isatty()
	pyRootPwa.utils.stderrisatty = sys.stderr.isatty()

	f = pyRootPwa.fileManager()
	f.dataDirectory = "/home/kbicker/sandbox/rootpwaV3/rpwaV3_data"

	types = [ pyRootPwa.core.eventMetadata.REAL, pyRootPwa.core.eventMetadata.GENERATED, pyRootPwa.core.eventMetadata.ACCEPTED ]
	typeToStringConversion = { types[0]: "", types[1]: ".genPS", types[2]: ".accPS" }
	inputFiles = { types[0]: [], types[1]: [], types[2]: []}
	for eventsType in types:
		for t in range(10):
			tpert = 0.01 if random.random() > 0.9 else 0.
			ltbinb = 0.1 + t*0.1
			htbinb = 0.1 + (t+1.)*0.1 + tpert
			for i in range(40):
				gap = 11. if random.random() > 2.95 else 0.
				pert = 111. if random.random() > 2.9 else 0.
				lmbinb = 1000. + i*60. + gap
				hmbinb = 1000. + (i+1.)*60. + pert
				inputFile = f.InputFile()
				inputFile.dataFileName = str(lmbinb) + "." + str(hmbinb) + "_" + str(ltbinb) + "." + str(htbinb) + typeToStringConversion[eventsType] + ".root"
				inputFile.eventsType = eventsType
				inputFile.binningMap = { "mass": (lmbinb, hmbinb), "tPrime": (ltbinb, htbinb) }
				if random.random() > 0.05:
					inputFiles[eventsType].append(inputFile)

#	inputFiles = f._openInputFiles()

#	pprint.pprint(inputFiles)
	allAxes = []
	for eventsType in inputFiles:
		allAxes.append(f._getBinningAxes(inputFiles[eventsType]))
	pprint.pprint(f._combineAxes(allAxes))
