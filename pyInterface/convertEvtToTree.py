#!/usr/bin/env python

import argparse
import os
import sys

import pyRootPwa
import pyRootPwa.core

class EventFile:

    evtfile = None

    def __init__(self, infile):
        self.evtfile = infile


    def get_event(self):
        n_lines_to_read = self.evtfile.readline()[:-1]
        if n_lines_to_read == "":
            return None
        lines = [n_lines_to_read]
        for i in range(0, int(n_lines_to_read)):
            lines.append(self.evtfile.readline()[:-1])
            if lines[-1] == "":
                raise Exception("Unexpected end of file")
        return Event(lines)


    def write_event(self, event):
        for line in event.lines:
            self.evtfile.write(line + "\n")

class Event:

    lines = []

    def __init__(self, lines):
        self.lines = lines

    def sort(self):
        new_lines = self.lines[2:]
        new_lines = sorted(new_lines, key=lambda entry: int(entry.split()[0]))
        self.lines = self.lines[0:2] + new_lines

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
    parser.add_argument("infile", help="The .evt file to be read")
    parser.add_argument("-o", "--output", help="The .root file to be written", default="output.root")
    parser.add_argument("-u", "--userstring", help="User string", default="")
    parser.add_argument("-t", "--type", dest="eventsTypeString", help="type of data (can be 'real', 'generated' or 'accepted', default: 'other')", default="other")

    args = parser.parse_args()

    printErr = pyRootPwa.utils.printErr
    printWarn = pyRootPwa.utils.printWarn
    printInfo = pyRootPwa.utils.printInfo

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
        printErr("type '" + eventsTypeString + "' is invalid as an event data type.")

    inputFile = args.infile

    outputFile = pyRootPwa.ROOT.TFile.Open(args.output, "NEW")
    if not outputFile:
        printErr("could not open output file. Aborting...")
        sys.exit(1)

    with open(inputFile, 'r') as infile:
        printInfo("Opened input file " + inputFile)
        in_event_file = EventFile(infile)

        event = in_event_file.get_event()

        fileWriter = pyRootPwa.core.eventFileWriter()
        success = fileWriter.initialize(outputFile,
                                        args.userstring,
                                        eventsType,
                                        event.getProductionKinematicsParticleNames(),
                                        event.getDecayKinematicsParticleNames(),
                                        {},
                                        [])

        while event is not None and success:
            event.sort()
            fileWriter.addEvent(event.getProductionKinematicsMomenta(), event.getDecayKinematicsMomenta())
            event = in_event_file.get_event()

    fileWriter.finalize()
