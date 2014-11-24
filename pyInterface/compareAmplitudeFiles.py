#!/usr/bin/env python

import argparse
import sys

import pyRootPwa
import pyRootPwa.core

if __name__ == "__main__":

	parser = argparse.ArgumentParser(
	                                 description="compare amplitude files"
	                                )

	parser.add_argument("file1", type=str, metavar="<filename>", help="path to first amplitude file")
	parser.add_argument("file2", type=str, metavar="<filename>", help="path to second amplitude file")
	args = parser.parse_args()

	pyRootPwa.utils.stdoutisatty = sys.stdout.isatty()
	pyRootPwa.utils.stderrisatty = sys.stderr.isatty()
	ROOT = pyRootPwa.ROOT

	file1 = ROOT.TFile(args.file1, "READ")
	file2 = ROOT.TFile(args.file2, "READ")
	tree1 = None
	tree2 = None
	for key in file1.GetListOfKeys():
		if key.GetName()[-3:] == "amp": tree1 = file1.Get(key.GetName())
	for key in file2.GetListOfKeys():
		if key.GetName()[-3:] == "amp": tree2 = file2.Get(key.GetName())

	if not tree1.GetEntries() == tree2.GetEntries():
		print "both trees have to have the same number of entries. Aborting..."
		sys.exit(1)

	real1 = []
	imag1 = []

	for event in tree1:
		leaf1  = event.__getattr__(pyRootPwa.core.amplitudeMetadata.amplitudeLeafName)
		real1.append(leaf1.amp().real())
		imag1.append(leaf1.amp().imag())

	real2 = []
	imag2 = []

	for event in tree2:
		leaf2  = event.__getattr__(pyRootPwa.core.amplitudeMetadata.amplitudeLeafName)
		real2.append(leaf2.amp().real())
		imag2.append(leaf2.amp().imag())

	diff = []
	diffSum = 0

	for i in range(len(real1)):
		diff.append((real1[i] - real2[i])**2 + (imag1[i] - imag2[i])**2)
		diffSum += diff[i]
	print "Sum of differences squared: " + str(diffSum)
		