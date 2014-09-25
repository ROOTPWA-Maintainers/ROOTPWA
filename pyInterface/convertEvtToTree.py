#!/usr/bin/env python

import argparse
import os
import tempfile

import pyRootPwa

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


    def __str__(self):
        retval = ""
        for line in self.lines:
            retval += line + '\n'
        return retval[:-1]

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
                                      description="Order particles in .evt files.",
                                      epilog="BE AWARE: If the output parameter is not given, "
                                             "the original file will be overwritten with the ordered version!"
                                    )
    parser.add_argument("infile", help="The .evt file to be read", nargs='?', default='/nfs/mnemosyne/user/odrotleff/private/proj/rootpwaMaster/ROOTPWA/DATA/1120.1190/1120.1190.evt')
    parser.add_argument("-o", "--output", help="The .evt file to be written", default='./test.root')

    args = parser.parse_args()

    input_file = args.infile
    output_file = args.output
    
    print(input_file + " " + output_file)

    if output_file is None:
        output_file = tempfile.mkstemp()[1]
        change_evt_file = Truep

    print("Starting to order events...")
    
    fileWriter = pyRootPwa.core.eventFileWriter()
    fileWriter.initialize(output_file, "test", ["test"], ["test"])
    
    with open(input_file, 'r') as infile:
        print("Opened input file " + input_file)
        in_event_file = EventFile(infile)
        
        event = in_event_file.get_event()
        while event is not None:
            event.sort()
            out_event_file.write_event(event)
            event = in_event_file.get_event()
