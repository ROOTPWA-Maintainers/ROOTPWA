#!/usr/bin/env python

import argparse
import multiprocessing
import sys

import pyRootPwa
ROOT = pyRootPwa.ROOT

if __name__ == "__main__":

	pyRootPwa.core.printCompilerInfo()
	pyRootPwa.core.printLibraryInfo()
	pyRootPwa.core.printGitHash()

	#initialize the printing functors
	printingCounter = multiprocessing.Array('i', [0]*5)
	pyRootPwa.utils.printErr = pyRootPwa.utils.printErrClass(printingCounter)
	pyRootPwa.utils.printWarn = pyRootPwa.utils.printWarnClass(printingCounter)
	pyRootPwa.utils.printSucc = pyRootPwa.utils.printSuccClass(printingCounter)
	pyRootPwa.utils.printInfo = pyRootPwa.utils.printInfoClass(printingCounter)
	pyRootPwa.utils.printDebug = pyRootPwa.utils.printDebugClass(printingCounter)


	parser = argparse.ArgumentParser(description="print event metadata")
	parser.add_argument("inputFile", type=str, metavar="inputFile", help="input file in ROOTPWA format")
	parser.add_argument("-r", action="store_true", dest="recalculateHash", help="recalculate hash and compare it with the stored one (default: false)")
	args = parser.parse_args()

	inputFile = ROOT.TFile.Open(args.inputFile, "READ")
	if not inputFile:
		pyRootPwa.utils.printErr("could not open input file '" + args.inputFile + "'.")

	eventMeta = pyRootPwa.core.eventMetadata.readEventFile(inputFile, True)
	if not eventMeta:
		pyRootPwa.utils.printErr("could not read event metadata")
		sys.exit(1)

	if args.recalculateHash:
		pyRootPwa.utils.printInfo("recalculating hash...")
		calcHash = eventMeta.recalculateHash(True)
		if calcHash != eventMeta.contentHash():
			pyRootPwa.utils.printErr("hash verification failed, hash from metadata '" +
			                         eventMeta.contentHash() +"' does not match with " +
			                         "calculated hash '" + calcHash + "'.")
		else:
			print("")
			pyRootPwa.utils.printSucc("recalculated hash matches with hash from metadata.")
			print("")

	print(eventMeta)
