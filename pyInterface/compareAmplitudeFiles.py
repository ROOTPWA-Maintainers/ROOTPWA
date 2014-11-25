#!/usr/bin/env python

import argparse
import sys
import numpy

import pyRootPwa
import pyRootPwa.core
ROOT = pyRootPwa.ROOT

def getDiffListAbs(list1, list2):
	diff = []
	for i in range(len(list1)):
		diff.append(list1[i] - list2[i])
	return diff

def getHistogram(name, title, lst):
	minHist = minimum(lst) * 1.1
	maxHist = maximum(lst) * 1.1
	if abs(maxHist) > abs(minHist):
		minHist = - abs(maxHist)
	else:
		maxHist = abs(minHist)
	if maxHist == 0: maxHist = 10**(-11)
	if minHist == 0: minHist = -10**(-11)
	minHist = -1e-14
	maxHist = 1e-14
	hist = ROOT.TH1D(name, title, 150, minHist, maxHist)
	for val in lst:
		hist.Fill(val)
	return hist

def saveHistogram(path, hist):
	rootFile = ROOT.TFile.Open(path, "RECREATE")
	hist.Write()
	rootFile.Close()

def saveTree(path, name, title, lst):
	rootFile = ROOT.TFile.Open(path, "RECREATE")
	tree = ROOT.TTree(name, title)
	arr = numpy.zeros(1)
	tree.Branch("diff", arr, "diff/D")
	for element in lst:
		arr[0] = element
		tree.Fill()
	hist = ROOT.TH1D(name + "_hist", title, 301, -1e-11, 1e-11)
	tree.Draw("diff>>" + name + "_hist")
	rootFile.Write()
	rootFile.Close()

def minimum(lst):
	currMin = lst[0]
	for val in lst[1:]:
		if val < currMin: currMin = val
	return currMin

def maximum(lst):
	currMax = lst[0]
	for val in lst[1:]:
		if val > currMax: currMax = val
	return currMax

def extractRealImagListsFromAmpFile(fileName):
	ampFile = ROOT.TFile(fileName, "READ")
	tree = None
	for currKey in ampFile.GetListOfKeys():
		if currKey.GetName()[-3:] == "amp": tree = ampFile.Get(currKey.GetName())

	real = []
	imag = []
	for currEvent in tree:
		leaf  = currEvent.__getattr__(pyRootPwa.core.amplitudeMetadata.amplitudeLeafName)
		real.append(leaf.amp().real())
		imag.append(leaf.amp().imag())
	return (real, imag)

if __name__ == "__main__":

	parser = argparse.ArgumentParser(
	                                 description="compare amplitude files"
	                                )

	parser.add_argument("file1", type=str, metavar="<filename>", help="path to first amplitude file")
	parser.add_argument("file2", type=str, metavar="<filename>", help="path to second amplitude file")
	args = parser.parse_args()

	pyRootPwa.utils.stdoutisatty = sys.stdout.isatty()
	pyRootPwa.utils.stderrisatty = sys.stderr.isatty()

	real1, imag1 = extractRealImagListsFromAmpFile(args.file1)
	real2, imag2 = extractRealImagListsFromAmpFile(args.file2)
	if not len(real1) == len(real2):
		pyRootPwa.utils.printErr("both trees have to have the same number of entries. Aborting...")
		sys.exit(1)

	diffAbsReal = getDiffListAbs(real1, real2)
	diffAbsImag = getDiffListAbs(imag1, imag2)

	histReal = getHistogram("deltaReal", "delta real", diffAbsReal)
	histImag = getHistogram("deltaImag", "delta imag", diffAbsImag)
	saveHistogram("histReal.root", histReal)
	saveHistogram("histImag.root", histImag)
	pyRootPwa.utils.printSucc("comparing the files '" + args.file1 + "' and '" + args.file2 + "' done.")
