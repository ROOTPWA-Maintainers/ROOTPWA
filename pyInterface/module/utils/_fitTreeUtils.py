import pyRootPwa.core

from _printingUtils import printErr, printWarn
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


def getFitResultsFromFiles(fitResultFileNames,
                           fitResultTreeName = "pwa",
                           fitResultBranchName = "fitResult_v2",
                           stripMatricesFromFurtherAttempts = False,
                           onlyConvergedResults = False
                          ):
	'''
	Loads all fit results from the fit-result file, grouped by their multibin.
	@param stripMatricesFromFurtherAttempts: if True, the integral and covariance matrices are striped
	                                         from each fit result except the best one and the best converged one in each multibin
	@param onlyConvergedResults: if True, only the converged fit results are loaded and considered for the best fit result

	@return: { <multiBin>: [ <fit-result> ] }
	         [ <fit-result> ] is ordered by logLikelihood from biggest to smallest
	'''

	return pyRootPwa.core.getFitResultsFromFilesInMultibins(fitResultFileNames, fitResultTreeName, fitResultBranchName,
	                                                        False, stripMatricesFromFurtherAttempts, onlyConvergedResults)


def getFitResultsFromFile(fitResultFileName,
                          fitResultTreeName = "pwa",
                          fitResultBranchName = "fitResult_v2",
                          stripMatricesFromFurtherAttempts = False,
                          onlyConvergedResults = False
                         ):
	'''
	Loads all fit results from the fit-result file, grouped by their multibin.
	@param stripMatricesFromFurtherAttempts: if True, the integral and covariance matrices are striped
	                                         from each fit result except the best one and the best converged one in each multibin
	@param onlyConvergedResults: if True, only the converged fit results are loaded and considered for the best fit result
	@return: { <multiBin>: [ <fit-result> ] }
	         [ <fit-result> ] is ordered by logLikelihood from biggest to smallest
	'''

	return getFitResultsFromFiles([fitResultFileName], fitResultTreeName, fitResultBranchName, stripMatricesFromFurtherAttempts, onlyConvergedResults)


def getBestFitResultsFromFile(fitResultFileName,
                              fitResultTreeName = "pwa",
                              fitResultBranchName = "fitResult_v2"
                             ):

	bestResultsInMultibin = pyRootPwa.core.getFitResultsFromFilesInMultibins([fitResultFileName], fitResultTreeName, fitResultBranchName, True, False, True)
	multiBins = set([ b.getSubMultiBin(exclude="mass") for b in bestResultsInMultibin.keys()])
	if len(multiBins) != 1:
		printWarn("Fit result file '{0}' contains more than one bin in a variable other than mass: {1}".format(fitResultFileName, list(multiBins)))
	return {b.getBinCenters()["mass"]: results[0] for b,results in bestResultsInMultibin.iteritems()}


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

	# only print the warning about multiple mass-bin centers in one file
	# once
	printedMassBinCenterWarning = False

	bestResult = None

	for i in xrange(fitResultTree.GetEntries()):
		fitResultTree.GetEntry(i)

		# skip fit results that did no converge
		if not result.converged():
			continue

		if bestResult is None:
			bestResult = pyRootPwa.core.fitResult(result)
		elif massBinCenter is None:
			if bestResult.massBinCenter() != result.massBinCenter() and not printedMassBinCenterWarning:
				printWarn("fit result file '" + fitResultFileName + "' " +
				          "contains more than one mass bin, return the " +
				          "fit result with the best likelihood.")
				printedMassBinCenterWarning = True
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
