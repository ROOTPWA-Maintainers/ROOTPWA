#!/usr/bin/python2.7

import argparse
import sys

import pyRootPwa

if __name__ == "__main__":

	parser = argparse.ArgumentParser(
	                                 description="calculates decay amplitudes "
	                                             "for given wave for events in "
	                                             "input data files and"
	                                             "writes amplitudes to file",
	                                 formatter_class=argparse.RawTextHelpFormatter
	                                )
	parser.add_argument("-k", type=str, metavar="file", required=True, dest="keyFileName", help="path to key file")
	parser.add_argument("-n", type=int, metavar="#", default=-1, dest="maxNmbEvents",  help="maximum number of events to read (default: all)")
	parser.add_argument("-p", type=str, metavar="file", default="./particleDataTable.txt", dest="pdgFileName", help="path to particle data table file (default: ./particleDataTable.txt)")
	parser.add_argument("-o", type=str, metavar="file", default="./out.root", dest="ampFileName", help="path to amplitude file (.amp or .root format; default: ./out.root)")
	parser.add_argument("-m", type=str, metavar="name", default="amplitude", dest="ampLeafName", help="amplitude leaf name (default: 'amplitude')")
	parser.add_argument("-a", action="store_true", dest="asciiOutput", help="write .amp files in ASCII format (default: binary)")
	parser.add_argument("-t", type=str, metavar="name", default="rootPwaEvtTree", dest="inTreeName", help="name of tree in ROOT data files (default: rootPwaEvtTree)")
	parser.add_argument("-l", type=str, metavar="names", default="prodKinParticles;prodKinMomenta;decayKinParticles;decayKinMomenta", dest="leafNames",
	                    help="semicolon separated object/leaf names in input data\n(default: 'prodKinParticles;prodKinMomenta;decayKinParticles;decayKinMomenta')")
	parser.add_argument("-v", action="store_true", dest="debug", help="verbose; print debug output (default: false)")

	parser.add_argument("-f", action="store_true", dest="newKeyFileNameConvention", help="use new keyfile name convention")

	parser.add_argument("input_files", nargs="+", help="input data file(s) (.evt or .root format)")
	arguments = parser.parse_args()

	pyRootPwa.printCompilerInfo()
	pyRootPwa.printLibraryInfo()
	pyRootPwa.printSvnVersion()

	inputFiles = []

	# check input file extensions
	for filename in arguments.input_files:
		if filename.rfind(".") < 0:
			pyRootPwa.printWarn('file extension missing in filename "' + filename + '".')
			continue
		fileExt = filename[filename.rfind(".")+1:]
		if (fileExt != "root") and (fileExt != "evt"):
			pyRootPwa.printWarn('input file "' + filename + '" is neither a .root nor a .evt file. skipping.')
			continue
		inputFiles.append(filename)

	# check output file extension
	if arguments.ampFileName.rfind(".") < 0:
		pyRootPwa.printErr('amplitude file "' + arguments.ampFileName + '" has no file extension. aborting.')
		sys.exit(1)
	fileExt = arguments.ampFileName[arguments.ampFileName.rfind(".")+1:]
	writeRootFile = False
	if fileExt == "root":
		writeRootFile = True
	elif fileExt != "amp":
		pyRootPwa.printErr('amplitude file "' + arguments.ampFileName + '" is neither a .root nor a .amp file. aborting.')
		sys.exit(1)

	# get object and leaf names
	(prodKinPartNamesObjName,
	 prodKinMomentaLeafName,
	 decayKinPartNamesObjName,
	 decayKinMomentaLeafName) = arguments.leafNames.split(";")

#	# open .root and .evt input files
#	inTrees = []
#	inChain = pyRootPwa.ROOT.TChain(arguments.inTreeName)
#
#	for filename in inputFiles:
#		if(filename[filename.rfind(".")+1:] == "root"):
#			pyRootPwa.printInfo('opening ROOT input file "' + filename + '".')
#			if(inChain.Add(filename) < 1):
#				pyRootPwa.printWarn('no events in ROOT input file "' + filename + '".')
#		else:
#			with open(filename, "r") is infile:
#				prodNames = pyRootPwa.ROOT.TClonesArray("TObjString")
#				decayNames = pyRootPwa.ROOT.TClonesArray("TObjString")


	# initialize the particleDataTable
	pyRootPwa.particleDataTable.readFile(arguments.pdgFileName)

	# Parse the keyfile and get the amplitude
	waveDesc = pyRootPwa.waveDescription()
	if not waveDesc.parseKeyFile(arguments.keyFileName):
		pyRootPwa.printErr('problems reading key file "' + arguments.keyFileName + '". aborting')
		sys.exit(1)

	(waveDescConstructionSuccess, amplitude) = waveDesc.constructAmplitude()
	if not waveDescConstructionSuccess:
		pyRootPwa.printErr('problems constructing decay topology from key file. aborting.')
		sys.exit(1)

	printString = 'creating amplitude file "' + arguments.ampFileName + '"'
	if writeRootFile:
		pyRootPwa.printInfo(printString)
		pyRootPwa.ROOT.TFile.Open(arguments.ampFileName, "RECREATE")
		amplitude.decayTopology()
		waveName = waveDesc.waveNameFromTopology(amplitude.decayTopology(), arguments.newKeyFileNameConvention, None)
		waveDesc.Write(waveName)
	else:
		if arguments.asciiOutput:
			printString += "; ASCII mode"
		else:
			printString += "; binary mode"

		pyRootPwa.printInfo(printString)
		pass

