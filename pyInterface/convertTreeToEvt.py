#!/usr/bin/env python

import argparse
import os
import sys

import pyRootPwa
import pyRootPwa.core

if __name__ == "__main__":

	parser = argparse.ArgumentParser(
	                                 description="Converts ROOTPWA .root file to .evt file."
	                                 )
	parser.add_argument("infile", help="The .root file to be read")
	parser.add_argument("-o", "--output", help="The .evt file to be written", default="output.evt")

	args = parser.parse_args()

	printErr = pyRootPwa.utils.printErr
	printWarn = pyRootPwa.utils.printWarn
	printInfo = pyRootPwa.utils.printInfo
	ROOT = pyRootPwa.ROOT

	inputFile = ROOT.TFile(args.infile, "READ")
	metaData = pyRootPwa.core.eventMetadata.readEventFile(inputFile)
	prodKinPartNames = metaData.productionKinematicsParticleNames()
	decayKinPartNames = metaData.decayKinematicsParticleNames()
	tree = metaData.eventTree()
	outputEvtFile = open(args.output, 'w')

	for t in tree:
		prodKinParticles  = t.prodKinMomenta
		decayKinParticles = t.decayKinMomenta
		particleCount = prodKinParticles.GetEntries() + decayKinParticles.GetEntries()
		outputEvtFile.write(str(particleCount) + '\n')

		for i in range(prodKinParticles.GetEntries()):
			outputEvtFile.write(
			                    str(pyRootPwa.core.particleDataTable.geantIdFromParticleName(prodKinPartNames[i])) + " 0 " + # TODO: ADD CHARGE
			                    str(prodKinParticles[i].Px()) + " " +
			                    str(prodKinParticles[i].Py()) + " " +
			                    str(prodKinParticles[i].Pz()) + " " +
			                    #str(prodKinParticles[i].M()) +
			                    "\n"
			                   )

		for i in range(decayKinParticles.GetEntries()):
			outputEvtFile.write(
			                    str(pyRootPwa.core.particleDataTable.geantIdFromParticleName(decayKinPartNames[i])) + " 0 " + # TODO: ADD CHARGE 
			                    str(decayKinParticles[i].Px()) + " " +
			                    str(decayKinParticles[i].Py()) + " " +
			                    str(decayKinParticles[i].Pz()) + " " +
			                    #str(decayKinParticles[i].M()) +
			                    "\n"
			                   )
	outputEvtFile.close()
	inputFile.Close()
