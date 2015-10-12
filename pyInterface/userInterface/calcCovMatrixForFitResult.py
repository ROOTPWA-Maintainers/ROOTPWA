#!/usr/bin/env python

import argparse
import sys

import pyRootPwa
import pyRootPwa.core
ROOT = pyRootPwa.ROOT


def getFitResultFromFile(fitResultFileName):
	fitResultFile = ROOT.TFile.Open(fitResultFileName, "READ")
	if not fitResultFile:
		pyRootPwa.utils.printErr("Could not open generated fit result file '" + fitResultFileName + "'.")
		return None
	fitResultTree = fitResultFile.Get("pwa")
	result = pyRootPwa.core.fitResult()
	result.setBranchAddress(fitResultTree, "fitResult_v2")
	if fitResultTree.GetEntries() != 1:
		pyRootPwa.utils.printErr("More than one fit result in TTree, somebody should probably implement this properly...")
		fitResultFile.Close()
		return None
	fitResultTree.GetEntry(0)
	if not result.converged():
		pyRootPwa.utils.printErr("Fit not converged for fit result file '" + fitResultFileName + "'.")
		fitResultFile.Close()
		return None
	return result


def getWaveDescThresFromFitResult(fitResult, keyFiles):
	waveDescThres = []
	for waveName in fitResult.waveNames():
		if waveName == "flat":
			continue

		waveDesc = waveDesc = pyRootPwa.core.waveDescription.parseKeyFile(keyFiles[waveName][0])[keyFiles[waveName][1]]

		thresholded = True
		for prodAmpIndex in xrange(fitResult.nmbProdAmps()):
			if fitResult.waveNameForProdAmp(prodAmpIndex) == waveName:
				if fitResult.prodAmp(prodAmpIndex) != 0.:
					tresholded = False

		threshold = 0.
		if thresholded:
			treshold = 1.1 * fitResult.massBinCenter()

		waveDescThres.append( (waveName, waveDesc, threshold) )
	return waveDescThres


if __name__ == "__main__":

	parser = argparse.ArgumentParser(
	                                 description="pwa fit executable"
	                                )

	parser.add_argument("inputFileName", type=str, metavar="inputFile", help="path to fit result")
	parser.add_argument("outputFileName", type=str, metavar="outputFile", help="path to output file")
	parser.add_argument("-c", type=str, metavar="configFileName", dest="configFileName", default="./rootpwa.config", help="path to config file (default: './rootpwa.config')")
	parser.add_argument("-b", type=int, metavar="#", dest="binID", default=0, help="bin ID of fit (default: 0)")
	parser.add_argument("-C", "--cauchyPriors", help="use half-Cauchy priors (default: false)", action="store_true")
	parser.add_argument("-P", "--cauchyPriorWidth", type=float, metavar ="WIDTH", default=0.5, help="width of half-Cauchy prior (default: 0.5)")
	parser.add_argument("-A", type=int, metavar="#", dest="accEventsOverride", default=0, help="number of input events to normalize acceptance to (default: use number of events from normalization integral file)")
	parser.add_argument("-v", "--verbose", help="verbose; print debug output (default: false)", action="store_true")
	args = parser.parse_args()

	config = pyRootPwa.rootPwaConfig()
	if not config.initialize(args.configFileName):
		pyRootPwa.utils.printErr("loading config file '" + args.configFileName + "' failed. Aborting...")
		sys.exit(1)
	pyRootPwa.core.particleDataTable.readFile(config.pdgFileName)
	fileManager = pyRootPwa.loadFileManager(config.fileManagerPath)
	if not fileManager:
		pyRootPwa.utils.printErr("loading the file manager failed. Aborting...")
		sys.exit(1)

	ampFileList = fileManager.getAmplitudeFilePaths(args.binID, pyRootPwa.core.eventMetadata.REAL)
	if not ampFileList:
		pyRootPwa.utils.printErr("could not retrieve valid amplitude file list. Aborting...")
		sys.exit(1)
	binningMap = fileManager.getBinFromID(args.binID)

	psIntegralPath  = fileManager.getIntegralFilePath(args.binID, pyRootPwa.core.eventMetadata.GENERATED)
	accIntegralPath = fileManager.getIntegralFilePath(args.binID, pyRootPwa.core.eventMetadata.ACCEPTED)

	result = getFitResultFromFile(args.inputFileName)
	if not result:
		pyRootPwa.utils.printErr("could not get fit result from file '" + args.inputFileName + "'. Aborting...")
		sys.exit(1)

	waveDescThres = getWaveDescThresFromFitResult(result, fileManager.getKeyFiles())
	if not waveDescThres:
		pyRootPwa.utils.printErr("error while getting wave names, descriptions and thresholds. Aborting...")
		sys.exit(1)

	likelihood = pyRootPwa.initLikelihood(waveDescThres = waveDescThres,
	                                      massBinCenter = result.massBinCenter(),
	                                      ampFileList = ampFileList,
	                                      normIntegralFileName = psIntegralPath,
	                                      accIntegralFileName = accIntegralPath,
	                                      accEventsOverride = args.accEventsOverride,
	                                      cauchy = args.cauchyPriors,
	                                      cauchyWidth = args.cauchyPriorWidth,
	                                      rank = result.rank(),
	                                      verbose = args.verbose)

	pars = []
	for i in range(likelihood.nmbPars()):
		parName = likelihood.parName(i);
		pars.append(result.fitParameter(parName))

	covMatrix = likelihood.CovarianceMatrix(pars)
	if args.verbose:
		covMatrix.Print()

	newResult = pyRootPwa.core.fitResult()
	newResult.fill(result.nmbEvents(),
	               result.normNmbEvents(),
	               result.massBinCenter(),
	               result.logLikelihood(),
	               result.rank(),
	               result.prodAmps(),
	               result.prodAmpNames(),
	               covMatrix,
	               result.fitParCovIndices(),
	               result.normIntegralMatrix(),
	               result.acceptedNormIntegralMatrix(),
	               result.phaseSpaceIntegralVector(),
	               result.converged(),
	               True)
	valTreeName   = "pwa"
	valBranchName = "fitResult_v2"
	outputFile = pyRootPwa.ROOT.TFile.Open(args.outputFileName, "NEW")
	if ((not outputFile) or outputFile.IsZombie()):
		pyRootPwa.utils.printErr("cannot open output file '" + args.outputFileName + "'. Aborting...")
		sys.exit(1)
	pyRootPwa.utils.printInfo("file '" + args.outputFileName + "' is empty. "
	        + "creating new tree '" + valTreeName + "' for PWA result.")
	tree = pyRootPwa.ROOT.TTree(valTreeName, valTreeName)
	if not newResult.branch(tree, valBranchName):
		pyRootPwa.utils.printErr("failed to create new branch '" + valBranchName + "' in file '" + args.outputFileName + "'.")
		sys.exit(1)
	tree.Fill()
	nmbBytes = tree.Write()
	outputFile.Close()
	if nmbBytes == 0:
		pyRootPwa.utils.printErr("problems writing fit result to TKey 'fitResult' "
		       + "in file '" + args.outputFileName + "'")
		sys.exit(1)
	else:
		pyRootPwa.utils.printSucc("wrote fit result to TKey 'fitResult' "
		        + "in file '" + args.outputFileName + "'")
