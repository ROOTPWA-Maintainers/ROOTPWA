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
		outputEvtFile.write(
		                    str(pyRootPwa.core.particleDataTable.geantIdFromParticleName(particleName)) + " " +
		                    str(charge) + " " +
		                    '%.16e' % particleMomentum.Px() + " " +
		                    '%.16e' % particleMomentum.Py() + " " +
		                    '%.16e' % particleMomentum.Pz() + " " +
		                    '%.16e' % energy + "\n"
		                   )
	else:
		printErr("particle (" + particleName + ") not found in particleDataTable. Aborting...")
		sys.exit()

if __name__ == "__main__":

	parser = argparse.ArgumentParser(
	                                 description="Converts ROOTPWA .root file to .evt file."
	                                 )
	parser.add_argument("inputFile", help="The .root file to be read")
	parser.add_argument("-o", "--outputFile", help="The .evt file to be written", default="output.evt")
	parser.add_argument("-p", "--particleDataTable", help="Path of particleDataTable file", default='$ROOTPWA/particleData/particleDataTable.txt')

	args = parser.parse_args()

	printWarn = pyRootPwa.utils.printWarn
	printErr = pyRootPwa.utils.printErr
	printSucc = pyRootPwa.utils.printSucc
	ROOT = pyRootPwa.ROOT

	if not pyRootPwa.core.particleDataTable.instance.readFile(args.particleDataTable):
		printErr("error loading particleDataTable. Aborting...")
		sys.exit(1)

	inputFile = ROOT.TFile(args.inputFile, "READ")
	if not inputFile:
		printErr("error opening input file. Aborting...")
		sys.exit(1)
	metaData = pyRootPwa.core.eventMetadata.readEventFile(inputFile)
	if metaData == 0:
		printErr("error reading metaData. Input file is not a RootPWA root file.")
	prodKinPartNames = metaData.productionKinematicsParticleNames()
	decayKinPartNames = metaData.decayKinematicsParticleNames()
	tree = metaData.eventTree()
	outputFileName = args.outputFile
	if outputFileName == "output.evt":
		binningMap = metaData.binningMap()
		if len(binningMap) > 0:
			for binningVariable in binningMap:
				outputFileName = binningVariable & "(" & binningMap[binningVariable][0] & "-" & binningMap[binningVariable][1] & ")"
			outputFileName += ".evt"

	outputEvtFile = open(outputFileName, 'w')
	particleCount = len(prodKinPartNames) + len(decayKinPartNames)

	for event in tree:
		prodKinMomenta  = event.prodKinMomenta
		decayKinMomenta = event.decayKinMomenta
		if particleCount != prodKinMomenta.GetEntries() + decayKinMomenta.GetEntries():
			printWarn("particle count in metaData does not match particle count in event data.")
		outputEvtFile.write(str(particleCount) + '\n')

		for particle in range(prodKinMomenta.GetEntries()):
			writeParticleToFile(outputEvtFile, prodKinPartNames[particle], prodKinMomenta[particle])

		for particle in range(decayKinMomenta.GetEntries()):
			writeParticleToFile(outputEvtFile, decayKinPartNames[particle], decayKinMomenta[particle])

	inputFile.Close()
	outputEvtFile.close()
	printSucc("successfully converted '" + args.inputFile + "' to '" + args.outputFile + "'.")
