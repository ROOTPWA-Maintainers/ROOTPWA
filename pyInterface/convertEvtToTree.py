#!/usr/bin/env python

import argparse
import os
import sys

import pyRootPwa
import pyRootPwa.core

_geantParticleCodes = {}
_geantParticleCodes[0] = "unknown"
_geantParticleCodes[1] = "gamma"
_geantParticleCodes[2] = "e"
_geantParticleCodes[3] = "e"
_geantParticleCodes[7] = "pi0"
_geantParticleCodes[8] = "pi"
_geantParticleCodes[9] = "pi"
_geantParticleCodes[11] = "K"
_geantParticleCodes[12] = "K"
_geantParticleCodes[13] = "n"
_geantParticleCodes[14] = "p"
_geantParticleCodes[15] = "pbar"
_geantParticleCodes[16] = "K0"
_geantParticleCodes[17] = "eta"
_geantParticleCodes[18] = "lambda"
_geantParticleCodes[57] = "rho(770)"
_geantParticleCodes[58] = "rho(770)"
_geantParticleCodes[59] = "rho(770)"
_geantParticleCodes[60] = "omega(782)"
_geantParticleCodes[61] = "eta'(958)"
_geantParticleCodes[62] = "phi(1020)"
_geantParticleCodes[45] = "d"

def getParticleNameFromGeantParticleID(id):
    try:
        return _geantParticleCodes[int(id)]
    except (KeyError,ValueError):
        return False

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
        particleName = getParticleNameFromGeantParticleID(splitUp[0])
        if not particleName == False:
            retVar.append(particleName)
        else:
            printErr("invalid particle ID (" + splitUp[0] + ").")
            sys.exit(1)
        return retVar

    def getDecayKinematicsParticleNames(self):
        retVar = []
        for line in self.lines[2:]:
            splitUp = line.split()
            particleName = getParticleNameFromGeantParticleID(splitUp[0])
            if not particleName == False:
                retVar.append(particleName)
            else:
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

    outputFile = pyRootPwa.ROOT.TFile(args.output, "NEW")
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
