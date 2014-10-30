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

	f = pyRootPwa.fileManager("../DATA/rpwaV3_data")
	#"/home/kbicker/sandbox/rootpwaV3/rpwaV3_data"

	#print f.getBinFromID(0)
	#print f.getBinFromID(1)
	#print f.getBinFromID(20)
	#pprint.pprint(f.binList)
	pprint.pprint(f.getDataFilePaths())
	print "# of files: " + str(len(f.getDataFilePaths()))
	#print f.getBinID({"mass":3361})
	missingBins = f.getMissingBins()
	pprint.pprint(missingBins)
	for etype in missingBins:
		print "# of bins missing in " + str(etype) + ": " + str(len(missingBins[etype]))

	#print f.getDataFilePath({'mass':3000}, pyRootPwa.core.eventMetadata.eventsTypeEnum.GENERATED))
