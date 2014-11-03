#!/usr/bin/env python

import argparse
import os
import sys
import math

import pyRootPwa
import pyRootPwa.core

def writeParticleToFile (file, particleName, particleMomentum):
	if pyRootPwa.core.particleDataTable.isInTable(particleName):
		partProperties = pyRootPwa.core.particleDataTable.entry(particleName)
		charge = partProperties.charge
		energy = math.sqrt(particleMomentum.Px()**2 + particleMomentum.Py()**2 + particleMomentum.Pz()**2 + partProperties.mass2)
		file.write(
		                    str(pyRootPwa.core.particleDataTable.geantIdFromParticleName(particleName)) + " " +
		                    str(charge) + " " +
		                    '%.16e' % particleMomentum.Px() + " " +
		                    '%.16e' % particleMomentum.Py() + " " +
		                    '%.16e' % particleMomentum.Pz() + " " +
		                    '%.16e' % energy + "\n"
		                   )
		return True
	else:
		pyRootPwa.utils.printErr("particle '" + particleName + "' not found in particleDataTable.")
		return False

if __name__ == "__main__":

	parser = argparse.ArgumentParser(
	                                 description="Converts ROOTPWA .root file to .evt file."
	                                 )
	parser.add_argument("inputFileName", help="The path to the RootPwa input file")
	parser.add_argument("outputFileName", help="The path to the ASCII evt output file")
	parser.add_argument("-p", "--particleDataTable", help="The path of particleDataTable file (default: '$ROOTPWA/particleData/particleDataTable.txt')",
	                    default='$ROOTPWA/particleData/particleDataTable.txt')

	args = parser.parse_args()

	printWarn = pyRootPwa.utils.printWarn
	printErr = pyRootPwa.utils.printErr
	printSucc = pyRootPwa.utils.printSucc
	ROOT = pyRootPwa.ROOT

	pdtPath = os.path.expandvars(args.particleDataTable)
	if not pyRootPwa.core.particleDataTable.instance.readFile(pdtPath):
		printErr("error loading particleDataTable from '" + pdtPath + "'. Aborting...")
		sys.exit(1)

	inputFile = ROOT.TFile(args.inputFileName, "READ")
	if not inputFile:
		printErr("error opening input file. Aborting...")
		sys.exit(1)
	metaData = pyRootPwa.core.eventMetadata.readEventFile(inputFile)
	if metaData == 0:
		printErr("error reading metaData. Input file is not a RootPWA root file.")
	prodKinPartNames = metaData.productionKinematicsParticleNames()
	decayKinPartNames = metaData.decayKinematicsParticleNames()
	tree = metaData.eventTree()

	with open(args.outputFileName, 'w') as outputEvtFile:
		particleCount = len(prodKinPartNames) + len(decayKinPartNames)
		for event in tree:
			prodKinMomenta  = event.__getattr__(metaData.productionKinematicsMomentaBranchName)
			decayKinMomenta = event.__getattr__(metaData.decayKinematicsMomentaBranchName)
			if particleCount != (prodKinMomenta.GetEntries() + decayKinMomenta.GetEntries()):
				printErr("particle count in metaData does not match particle count in event data.")
				sys.exit(1)
			outputEvtFile.write(str(particleCount) + '\n')

			for particle in range(prodKinMomenta.GetEntries()):
				if not writeParticleToFile(outputEvtFile, prodKinPartNames[particle], prodKinMomenta[particle]):
					printErr("failed writing particle '" + particle + "' to output file.")
					sys.exit(1)

			for particle in range(decayKinMomenta.GetEntries()):
				if not writeParticleToFile(outputEvtFile, decayKinPartNames[particle], decayKinMomenta[particle]):
					printErr("failed writing particle '" + particle + "' to output file.")
					sys.exit(1)

		inputFile.Close()
	printSucc("successfully converted '" + args.inputFileName + "' to '" + args.outputFileName + "'.")
