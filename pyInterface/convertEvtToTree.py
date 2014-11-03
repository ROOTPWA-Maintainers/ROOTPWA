#!/usr/bin/env python

import argparse
import sys

import pyRootPwa
import pyRootPwa.core

class EventFile:

	evtFile = None

	def __init__(self, inputFile):
		self.evtFile = inputFile


	def getEvent(self):
		linesToRead = self.evtFile.readline()[:-1]
		if linesToRead == "":
			return None
		lines = [linesToRead]
		for i in range(0, int(linesToRead)):
			lines.append(self.evtFile.readline()[:-1])
			if lines[-1] == "":
				raise Exception("Unexpected end of file")
		return Event(lines)


	def writeEvent(self, event):
		for line in event.lines:
			self.evtFile.write(line + "\n")

class Event:

	lines = []

	def __init__(self, lines):
		self.lines = lines

	def reorder(self, targetProdNames, targetDecayNames):
		currentProdNames  = self.getProductionKinematicsParticleNames()
		currentDecayNames = self.getDecayKinematicsParticleNames()
		if currentProdNames[0] != targetProdNames[0]:
			printErr("production kinematics particles do not match over all events. Aborting...")
			sys.exit(1)
		newlines = []
		newlines.append(self.lines[0])
		newlines.append(self.lines[1])

		for i in range(len(targetDecayNames)):
			for j in range(len(currentDecayNames)):
				if currentDecayNames[j] == targetDecayNames[i]:
					newlines.append(self.lines[j+2])
					currentDecayNames[j]=""
					break
		self.lines = newlines

	def getProductionKinematicsParticleNames(self):
		retVar = []
		splitUp = self.lines[1].split()
		try:
			particleName = pyRootPwa.core.particleDataTable.particleNameFromGeantId(int(splitUp[0]))
			retVar.append(particleName)
		except:
			printErr("invalid particle ID (" + str(splitUp[0]) + ").")
			sys.exit(1)
		return retVar

	def getDecayKinematicsParticleNames(self):
		retVar = []
		for line in self.lines[2:]:
			splitUp = line.split()
			try:
				particleName = pyRootPwa.core.particleDataTable.particleNameFromGeantId(int(splitUp[0]))
				retVar.append(particleName)
			except:
				printErr("invalid particle ID (" + splitUp[0] + ").")
				sys.exit(1)
		return retVar

	def getProductionKinematicsMomenta(self):
		retVar = []
		splitUp = self.lines[1].split()
		retVar.append(pyRootPwa.ROOT.TVector3(float(splitUp[2]), float(splitUp[3]), float(splitUp[4])))
		return retVar

	def getDecayKinematicsMomenta(self):
		retVar = []
		for line in self.lines[2:]:
			splitUp = line.split()
			retVar.append(pyRootPwa.ROOT.TVector3(float(splitUp[2]), float(splitUp[3]), float(splitUp[4])))
		return retVar

	def __str__(self):
		retval = ""
		for line in self.lines:
			retval += line + '\n'
		return retval[:-1]

if __name__ == "__main__":

	parser = argparse.ArgumentParser(
	                                  description="Converts .evt file to .root file."
	                                )
	parser.add_argument("inputFileName", help="The .evt file to be read")
	parser.add_argument("outputFileName", help="The .root file to be written")
	parser.add_argument("-u", "--userstring", help="User string", default="")
	parser.add_argument("-t", "--type", dest="eventsTypeString", help="type of data (can be 'real', 'generated' or 'accepted', default: 'other')", default="other")
	parser.add_argument("-b", "--binning", action='append', help="declare current bin in the form 'binningVariable;lowerBound;upperBound' (e.g. 'mass;1000;1100')."+
	                                                             "You can use the argument multiple times for multiple binning variables")

	args = parser.parse_args()

	printErr = pyRootPwa.utils.printErr
	printWarn = pyRootPwa.utils.printWarn
	printInfo = pyRootPwa.utils.printInfo
	printSucc = pyRootPwa.utils.printSucc

	eventsType = pyRootPwa.core.eventMetadata.OTHER
	if args.eventsTypeString == "real":
		eventsType = pyRootPwa.core.eventMetadata.REAL
	elif args.eventsTypeString == "generated":
		eventsType = pyRootPwa.core.eventMetadata.GENERATED
	elif args.eventsTypeString == "accepted":
		eventsType = pyRootPwa.core.eventMetadata.ACCEPTED
	elif args.eventsTypeString == "other":
		pass
		# do nothing
	else:
		printErr("type '" + args.eventsTypeString + "' is invalid as an event data type.")
	printInfo("set eventsType to '" + str(eventsType) + "'.")
	printInfo("set userString to '" + args.userstring + "'.")

	outputFile = pyRootPwa.ROOT.TFile.Open(args.outputFileName, "NEW")
	if not outputFile:
		printErr("could not open output file '" + args.outputFileName + "'. Aborting...")
		sys.exit(1)

	binningMap = {}
	if args.binning:
		for bin in args.binning:
			try:
				splitUp = bin.split(";")
				if len(splitUp)==3:
					binningMap[splitUp[0]] = (float(splitUp[1]), float(splitUp[2]))
					printInfo("adding to binning map: " + splitUp[0] + " -> (" + splitUp[1] + "," + splitUp[2] + ")")
				else:
					printErr("did not get the right amount of semicolon seperated values for " + splitUp[0] + "-bin. Aborting...")
					sys.exit(1)
			except ValueError:
				printErr("could not convert binning map boundaries of " + splitUp[0] + "-bin to float. Aborting...")
				sys.exit(1)

	with open(args.inputFileName, 'r') as inputFile:
		printInfo("Opened input file '" + args.inputFileName + "'.")
		eventFile = EventFile(inputFile)

		event = eventFile.getEvent()

		fileWriter = pyRootPwa.core.eventFileWriter()
		initialProdNames  = event.getProductionKinematicsParticleNames()
		initialDecayNames = event.getDecayKinematicsParticleNames()
		success = fileWriter.initialize(outputFile,
		                                args.userstring,
		                                eventsType,
		                                initialProdNames,
		                                initialDecayNames,
		                                binningMap,
		                                [])

		while event is not None and success:
			event.reorder(initialProdNames, initialDecayNames)
			fileWriter.addEvent(event.getProductionKinematicsMomenta(), event.getDecayKinematicsMomenta())
			event = eventFile.getEvent()

	fileWriter.finalize()
	printSucc("successfully converted '" + args.inputFileName + "' to '" + args.outputFileName + "'.")
