#!/usr/bin/env python

import argparse
import os
import sys

import pyRootPwa
ROOT = pyRootPwa.ROOT

if __name__ == "__main__":

	pyRootPwa.core.printCompilerInfo()
	pyRootPwa.core.printLibraryInfo()
	pyRootPwa.core.printGitHash()

	parser = argparse.ArgumentParser(description="print event metadata")
	parser.add_argument("inputFile", type=str, metavar="inputFile", help="input file in ROOTPWA format")
	parser.add_argument("-a", "--amplitudeObjectBaseName", type=str, metavar="objectBaseName", help="amplitude object base name (default: try to guess from input file name)")
	parser.add_argument("-r", action="store_true", dest="recalculateHash", help="recalculate hash and compare it with the stored one (default: %(default)s)")
	args = parser.parse_args()

	inputFile = ROOT.TFile.Open(args.inputFile, "READ")
	if not inputFile:
		pyRootPwa.utils.printErr("could not open input file '" + args.inputFile + "'.")
		sys.exit(1)

	# try to read event metadata
	eventMeta = pyRootPwa.core.eventMetadata.readEventFile(inputFile, True)
	if eventMeta:
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

	# try to read amplitude metadata
	amplitudeObjectBaseName = args.amplitudeObjectBaseName
	if not amplitudeObjectBaseName:
		amplitudeObjectBaseName = os.path.basename(args.inputFile).split("_binID")[0]
	amplitudeMeta = pyRootPwa.core.amplitudeMetadata.readAmplitudeFile(inputFile, amplitudeObjectBaseName, True)
	if amplitudeMeta:
		if args.recalculateHash:
			pyRootPwa.utils.printInfo("recalculating hash...")
			calcHash = amplitudeMeta.recalculateHash(True)
			if calcHash != amplitudeMeta.contentHash():
				pyRootPwa.utils.printErr("hash verification failed, hash from metadata '" +
				                         amplitudeMeta.contentHash() +"' does not match with " +
				                         "calculated hash '" + calcHash + "'.")
			else:
				print("")
				pyRootPwa.utils.printSucc("recalculated hash matches with hash from metadata.")
				print("")

		print(amplitudeMeta)

	# try to read amplitude integral matrix metadata
	ampIntegralMatrixMeta = pyRootPwa.core.ampIntegralMatrixMetadata.readIntegralFile(inputFile, True)
	if ampIntegralMatrixMeta:
		if args.recalculateHash:
			pyRootPwa.utils.printInfo("recalculating hash...")
			calcHash = ampIntegralMatrixMeta.recalculateHash()
			if calcHash != ampIntegralMatrixMeta.contentHash():
				pyRootPwa.utils.printErr("hash verification failed, hash from metadata '" +
				                         ampIntegralMatrixMeta.contentHash() +"' does not match with " +
				                         "calculated hash '" + calcHash + "'.")
			else:
				print("")
				pyRootPwa.utils.printSucc("recalculated hash matches with hash from metadata.")
				print("")

		print(ampIntegralMatrixMeta)

	if not eventMeta and not amplitudeMeta and not ampIntegralMatrixMeta:
		pyRootPwa.utils.printErr("could not read any metadata")
		sys.exit(1)
