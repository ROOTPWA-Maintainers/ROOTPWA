

import pyRootPwa
ROOT = pyRootPwa.ROOT


def getFitResultDict(fitResultTree):
	fitRes = pyRootPwa.core.fitResult()
	fitRes.setBranchAddress(fitResultTree, "fitResult_v2")
	retval = {}
	for i in xrange(fitResultTree.GetEntries()):
		fitResultTree.GetEntry(i)
		massBinCenter = fitRes.massBinCenter()
		resCopy = pyRootPwa.core.fitResult(fitRes)
		if massBinCenter not in retval:
			retval[massBinCenter] = [ ]
		retval[massBinCenter].append(resCopy)
	for key in retval:
		retval[key] = sorted(retval[key], key=lambda t: t.logLikelihood())
	return retval

def getBestFitResults(fitResultTree):
	retval = { }
	fitResultDict = getFitResultDict(fitResultTree)
	for massBinCenter, results in fitResultDict.iteritems():
		retval[massBinCenter] = results[0]
	return retval

