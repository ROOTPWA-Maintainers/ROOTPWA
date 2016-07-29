#!/usr/bin/env python

from __future__ import print_function

import argparse
import sys

import pyRootPwa
ROOT = pyRootPwa.ROOT

if __name__ == "__main__":

	pyRootPwa.core.printCompilerInfo()
	pyRootPwa.core.printLibraryInfo()
	pyRootPwa.core.printGitHash()

	parser = argparse.ArgumentParser(description="print event metadata")
	parser.add_argument("inputFile", type=str, metavar="inputFile", help="input file in ROOTPWA format")
	parser.add_argument("-v", "--verify", action="store_true", dest="verifyHash", help="verify hash(es) (default: %(default)s)")
	args = parser.parse_args()

	inputFile = ROOT.TFile.Open(args.inputFile, "READ")
	if not inputFile:
		pyRootPwa.utils.printErr("could not open input file '" + args.inputFile + "'.")
		sys.exit(1)

	# keep track that at least one metadata object was found
	readMeta = False

	# try to read event metadata
	eventMeta = pyRootPwa.core.eventMetadata.readEventFile(inputFile, True)
	if eventMeta:
		readMeta = True

		if args.verifyHash:
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
	for key in inputFile.GetListOfKeys():
		# get the amplitude names from the list of keys
		if key.GetClassName() == "rpwa::amplitudeMetadata":
			if not key.GetName().endswith(".meta"):
				pyRootPwa.utils.printErr("key with name '{}' is an amplitude metadata object but does not follow the naming convention.".format(key.GetName()))
				sys.exit(1)

			amplitudeMeta = pyRootPwa.core.amplitudeMetadata.readAmplitudeFile(inputFile, key.GetName()[:-5], True)

			readMeta = True

			if args.verifyHash:
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
		readMeta = True

		if args.verifyHash:
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

	if not readMeta:
		pyRootPwa.utils.printErr("could not read any metadata")
		sys.exit(1)
