#!/usr/bin/env python
##########################################################################
#
#    Copyright 2012
#
#    This file is part of rootpwa
#
#    rootpwa is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    rootpwa is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with rootpwa. If not, see <http://www.gnu.org/licenses/>.
#
##########################################################################
#-------------------------------------------------------------------------
#
# Description:
#      converts all .amp files in the given mass bin directories into ROOT files
#
#
# Author List:
#      Karl Bicker          TUM            (original author)
#
#
#-------------------------------------------------------------------------


import argparse
import os
import tempfile

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
	parser.add_argument("infile", help="The .evt file to be read")
	parser.add_argument("-o", "--output", help="The .evt file to be written")

	args = parser.parse_args()

	input_file = args.infile
	output_file = args.output

	change_evt_file = False
	if output_file is None:
		output_file = tempfile.mkstemp()[1]
		change_evt_file = True

	print("Starting to order events...")
	with open(input_file, 'r') as infile:
		print("Opened input file " + input_file)
		in_event_file = EventFile(infile)
		with open(output_file, 'w') as outfile:
			print("Opened output file " + output_file)
			out_event_file = EventFile(outfile)
			event = in_event_file.get_event()
			while event is not None:
				event.sort()
				out_event_file.write_event(event)
				event = in_event_file.get_event()
	if change_evt_file:
		print("Moving " + output_file + " to " + input_file)
		os.system('mv -f ' + output_file + ' ' + input_file)
