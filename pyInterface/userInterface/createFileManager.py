#!/usr/bin/env python

import argparse
import sys

import pyRootPwa
import pyRootPwa.core

if __name__ == "__main__":

	# parse command line arguments
	parser = argparse.ArgumentParser(
	                                 description="creates a file manager and saves it to the hard drive"
	                                )

	parser.add_argument("-c", type=str, metavar="configFileName", default="rootpwa.config", dest="configFileName", help="path to config file (default: ./rootpwa.config)")
	args = parser.parse_args()

	config = pyRootPwa.rootPwaConfig()
	if not config.initialize(args.configFileName):
		pyRootPwa.utils.printErr("loading config file '" + args.configFileName + "' failed. Aborting...")
		sys.exit(1)
	pyRootPwa.core.particleDataTable.readFile(config.pdgFileName)
	fileManager = pyRootPwa.fileManager()
	if not fileManager.initialize(config):
		pyRootPwa.utils.printErr("loading the file manager failed. Aborting...")
		sys.exit(1)
	if not pyRootPwa.saveFileManager(fileManager, config.fileManagerPath):
		pyRootPwa.utils.printErr("saving the file manager failed. Aborting...")
		sys.exit(1)
	pyRootPwa.utils.printSucc("Saved file manager to '" + config.fileManagerPath + "'.")