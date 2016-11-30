from _printingUtils import printErr, printWarn

import pyRootPwa.core

import _root
ROOT = _root.ROOT


def getFitResultFromFile(fitResultFileName,
                         fitResultTreeName = "pwa",
                         fitResultBranchName = "fitResult_v2"
                        ):
	fitResultFile = ROOT.TFile.Open(fitResultFileName, "READ")
	if not fitResultFile:
		printErr("Could not open generated fit result file '" + fitResultFileName + "'.")
		return None

	fitResultTree = fitResultFile.Get(fitResultTreeName)
	if not fitResultTree:
		printErr("could not find fit result tree '" + fitResultTreeName +
		                         "' in file '" + fitResultFileName + "'. Aborting...")
		return None

	result = pyRootPwa.core.fitResult()
	result.setBranchAddress(fitResultTree, fitResultBranchName)

	if fitResultTree.GetEntries() != 1:
		printErr("More than one fit result in TTree, somebody should probably implement this properly...")
		fitResultFile.Close()
		return None
	fitResultTree.GetEntry(0)

	if not result.converged():
		printErr("Fit not converged for fit result file '" + fitResultFileName + "'.")
		return None

	fitResultFile.Close()
	return result


def getBestFitResultsFromFile(fitResultFileName,
                              fitResultTreeName = "pwa",
                              fitResultBranchName = "fitResult_v2"
                             ):
	fitResultFile = ROOT.TFile.Open(fitResultFileName, "READ")
	if not fitResultFile:
		printErr("Could not open generated fit result file '" + fitResultFileName + "'.")
		return None

	fitResultTree = fitResultFile.Get(fitResultTreeName)
	if not fitResultTree:
		printErr("could not find fit result tree '" + fitResultTreeName +
		                         "' in file '" + fitResultFileName + "'. Aborting...")
		return None

	result = pyRootPwa.core.fitResult()
	result.setBranchAddress(fitResultTree, fitResultBranchName)

	bestResults = { }

	for i in xrange(fitResultTree.GetEntries()):
		fitResultTree.GetEntry(i)

		# skip fit results that did no converge
		if not result.converged():
			continue

		if not result.massBinCenter() in bestResults:
			bestResults[result.massBinCenter()] = pyRootPwa.core.fitResult(result)
		else:
			if result.logLikelihood() < bestResults[result.massBinCenter()].logLikelihood():
				bestResults[result.massBinCenter()] = pyRootPwa.core.fitResult(result)

	fitResultFile.Close()

	return bestResults


def getBestFitResultFromFile(fitResultFileName,
                             massBinCenter = None,
                             fitResultTreeName = "pwa",
                             fitResultBranchName = "fitResult_v2"
                            ):
	fitResultFile = ROOT.TFile.Open(fitResultFileName, "READ")
	if not fitResultFile:
		printErr("Could not open generated fit result file '" + fitResultFileName + "'.")
		return None

	fitResultTree = fitResultFile.Get(fitResultTreeName)
	if not fitResultTree:
		printErr("could not find fit result tree '" + fitResultTreeName +
		         "' in file '" + fitResultFileName + "'. Aborting...")
		return None

	result = pyRootPwa.core.fitResult()
	result.setBranchAddress(fitResultTree, fitResultBranchName)

	bestResult = None

	for i in xrange(fitResultTree.GetEntries()):
		fitResultTree.GetEntry(i)

		# skip fit results that did no converge
		if not result.converged():
			continue

		if bestResult is None:
			bestResult = pyRootPwa.core.fitResult(result)
		elif massBinCenter is None:
			if bestResult.massBinCenter() != result.massBinCenter():
				printWarn("fit result file '" + fitResultFileName + "' " +
				          "contains more than one mass bin, return the " +
				          "fit result with the best likelihood.")
			if result.logLikelihood() < bestResult.logLikelihood():
				bestResult = pyRootPwa.core.fitResult(result)
		else:
			if abs(massBinCenter - result.massBinCenter()) < abs(massBinCenter - bestResult.massBinCenter()):
				bestResult = pyRootPwa.core.fitResult(result)
			elif abs(massBinCenter - result.massBinCenter()) == abs(massBinCenter - bestResult.massBinCenter()):
				if result.logLikelihood() < bestResult.logLikelihood():
					bestResult = pyRootPwa.core.fitResult(result)

	fitResultFile.Close()

	return bestResult
