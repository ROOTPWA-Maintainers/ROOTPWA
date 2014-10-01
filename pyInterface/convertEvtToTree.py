#!/usr/bin/env python

import argparse
import os

import pyRootPwa

def getParticleNameFromGeantParticleID(number):
    if number==7:
        return "PION 0"
    elif number==8:
        return "PION +"
    elif number==9:
        return "PION -"
    else:
        return "ERROR"

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
        retVar.append(getParticleNameFromGeantParticleID(splitUp[0]))
        return retVar

    def getDecayKinematicsParticleNames(self):
        retVar = []
        for line in self.lines[2:]:
            splitUp = self.lines[1].split()
            retVar.append(getParticleNameFromGeantParticleID(splitUp[0]))
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
    parser.add_argument("-u", "--userstring", help="User string")

    args = parser.parse_args()

    input_file = args.infile
    output_file = os.path.abspath(args.output)
    output_file = pyRootPwa.ROOT.TFile(output_file, "NEW")

    with open(input_file, 'r') as infile:
        print("Opened input file " + input_file)
        in_event_file = EventFile(infile)

        event = in_event_file.get_event()

        fileWriter = pyRootPwa.core.eventFileWriter()
        success = fileWriter.initialize(output_file, args.userstring, event.getProductionKinematicsParticleNames(), event.getDecayKinematicsParticleNames(), {}, [])

        while event is not None and success:
            event.sort()
            fileWriter.addEvent(event.getProductionKinematicsMomenta(), event.getDecayKinematicsMomenta())
            event = in_event_file.get_event()

    fileWriter.finalize()
