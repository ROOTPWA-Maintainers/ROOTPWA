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
	if not f.initialize("rootpwa.config"):
		pyRootPwa.utils.printErr("initializing the file manager failed. Aborting...")
		sys.exit(1)
	#"/home/kbicker/sandbox/rootpwaV3/rpwaV3_data"

	#print f.getBinFromID(0)
	#print f.getBinFromID(1)
	#print f.getBinFromID(20)
	#pprint.pprint(f.binList)
	#pprint.pprint(f.getDataFilePaths())
	#print "# of files: " + str(len(f.getDataFilePaths()))
	#print f.getBinID({"mass":3361})
	#missingBins = f.getMissingBins()
	#pprint.pprint(missingBins)
	#for etype in missingBins:
	#	print "# of bins missing in " + str(etype) + ": " + str(len(missingBins[etype]))
	#print f.getDataFilePath({'mass':3000}, pyRootPwa.core.eventMetadata.eventsTypeEnum.GENERATED))
	#pprint.pprint(f.keyFiles)
	#print f.areDataFilesSynced()
	#print f.areKeyFilesSynced()

	for binID in f.getBinIDList():
		for keyFileID in f.getKeyFileIDList():
			for eventsType in [ pyRootPwa.core.eventMetadata.REAL, pyRootPwa.core.eventMetadata.GENERATED, pyRootPwa.core.eventMetadata.ACCEPTED ]:
				dataFile = f.getDataFile(binID, eventsType)
				keyFile = f.getKeyFile(keyFileID)
				if not dataFile:
					continue
				print str(binID) + "	" + str(keyFileID) + "	" + f.getAmplitudeFilePath(binID, keyFileID, eventsType)
				#print "submit job: calcAmplitudes(" + keyFile.keyFileName + ", " + dataFile.dataFileName + ", -o " + f.getAmplitudeFilePath(binID, keyFileID, eventsType) + ")"
	print ">> found " + str(len(f.getBinIDList())) + " binIDs and " + str(len(f.getKeyFilePaths())) + " keyFiles."
	print ">> submitted " + str(len(f.getBinIDList()) * len(f.getKeyFilePaths()) * len(f.dataFiles)) + " jobs."
