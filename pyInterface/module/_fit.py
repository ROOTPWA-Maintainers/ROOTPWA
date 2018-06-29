import random

import pyRootPwa.core
import pyRootPwa.utils
ROOT = pyRootPwa.utils.ROOT

def pwaFit(eventAndAmpFileDict,
           normIntegralFileName,
           accIntegralFileName,
           multiBin,
           waveListFileName,
           waveDescriptions,
           seed=0,
           cauchy=False,
           cauchyWidth=0.5,
           startValFileName="",
           accEventsOverride=0,
           useNormalizedAmps=True,
           checkHessian=False,
           saveSpace=False,
           rank=1,
           verbose=False,
           attempts=1
          ):

	waveDescThres = pyRootPwa.utils.getWaveDescThresFromWaveList(waveListFileName, waveDescriptions)
	massBinCenter = (multiBin.boundaries['mass'][1] + multiBin.boundaries['mass'][0]) / 2. # YOU CAN DO BETTER

	likelihood = pyRootPwa.initLikelihood(waveDescThres = waveDescThres,
	                                      massBinCenter = massBinCenter,
	                                      eventAndAmpFileDict = eventAndAmpFileDict,
	                                      normIntegralFileName = normIntegralFileName,
	                                      accIntegralFileName = accIntegralFileName,
	                                      multiBin = multiBin,
	                                      accEventsOverride = accEventsOverride,
	                                      useNormalizedAmps = useNormalizedAmps,
	                                      cauchy = cauchy,
	                                      cauchyWidth = cauchyWidth,
	                                      rank = rank,
	                                      verbose = verbose)
	if not likelihood:
		pyRootPwa.utils.printErr("error while initializing likelihood. Aborting...")
		return [ ]

	if attempts == 1:
		seeds = [ seed ]
	elif seed == 0:
		seeds = [ 0 for _ in xrange(attempts) ]
	else:
		random.seed(seed)
		seeds = [ ]
		for _ in xrange(attempts):
			while True:
				randVal = random.randint(1000, 2**32-1)
				if randVal not in seeds:
					break
			seeds.append(randVal)

	fitResults = [ ]
	for fitSeed in seeds:
		fitResult = pyRootPwa.core.pwaFit(likelihood         = likelihood,
		                                  multibinBoundaries = multiBin.boundaries,
		                                  seed               = fitSeed,
		                                  startValFileName   = startValFileName,
		                                  checkHessian       = checkHessian,
		                                  saveSpace          = saveSpace,
		                                  verbose            = verbose)
		fitResults.append(fitResult)
	return fitResults


def pwaNloptFit(eventAndAmpFileDict,
                normIntegralFileName,
                accIntegralFileName,
                multiBin,
                waveListFileName,
                waveDescriptions,
                seed=0,
                cauchy=False,
                cauchyWidth=0.5,
                startValFileName="",
                accEventsOverride=0,
                useNormalizedAmps=True,
                checkHessian=False,
                saveSpace=False,
                rank=1,
                verbose=False,
                attempts=1,
                keepMatricesOnlyOfBest= False
               ):

	if keepMatricesOnlyOfBest:
		saveSpace = True

	waveDescThres = pyRootPwa.utils.getWaveDescThresFromWaveList(waveListFileName, waveDescriptions)
	massBinCenter = (multiBin.boundaries['mass'][1] + multiBin.boundaries['mass'][0]) / 2. # YOU CAN DO BETTER

	likelihood = pyRootPwa.initLikelihood(waveDescThres = waveDescThres,
	                                      massBinCenter = massBinCenter,
	                                      eventAndAmpFileDict = eventAndAmpFileDict,
	                                      normIntegralFileName = normIntegralFileName,
	                                      accIntegralFileName = accIntegralFileName,
	                                      multiBin = multiBin,
	                                      accEventsOverride = accEventsOverride,
	                                      useNormalizedAmps = useNormalizedAmps,
	                                      cauchy = cauchy,
	                                      cauchyWidth = cauchyWidth,
	                                      rank = rank,
	                                      verbose = verbose)
	if not likelihood:
		pyRootPwa.utils.printErr("error while initializing likelihood. Aborting...")
		return [ ]

	if attempts == 1:
		seeds = [ seed ]
	elif seed == 0:
		seeds = [ 0 for _ in xrange(attempts) ]
	else:
		random.seed(seed)
		seeds = [ ]
		for _ in xrange(attempts):
			while True:
				randVal = random.randint(1000, 2**32-1)
				if randVal not in seeds:
					break
			seeds.append(randVal)

	fitResults = [ ]
	iBest = -1
	iBestValid = -1
	negLogLikeBest = 1e300
	negLogLikeBestValid = 1e300
	for fitSeed in seeds:
		fitResult = pyRootPwa.core.pwaNloptFit(likelihood         = likelihood,
		                                       multibinBoundaries = multiBin.boundaries,
		                                       seed               = fitSeed,
		                                       startValFileName   = startValFileName,
		                                       checkHessian       = checkHessian,
		                                       saveSpace          = saveSpace,
		                                       verbose            = verbose)
		fitResults.append(fitResult)
		iFitResult = len(fitResults)-1
		if fitResult.logLikelihood() < negLogLikeBest:
			negLogLikeBest = fitResult.logLikelihood()
			iBest = iFitResult
		if fitResult.converged() and fitResult.logLikelihood() < negLogLikeBestValid:
			negLogLikeBestValid = fitResult.logLikelihood()
			iBestValid = iFitResult

	if keepMatricesOnlyOfBest: # calculate covariance matrix only of best result and best valid result
		# get matrices
		norm, acc, psVector = likelihood.integralMatrices(True)
		fitResults[iBestValid] = _addMatrices(fitResults[iBestValid], likelihood, norm, acc, psVector)
		if iBest != iBestValid:
			fitResults[iBest] = _addMatrices(fitResults[iBest], likelihood, norm, acc, psVector)

	return fitResults


def addCovarianceMatrix(result, likelihood, verbose = False):
	'''
	Calculate parameter covariance matrix and adds it to a new fit result
	@return: New fit result containing the covariance matrix
	'''
	return _addMatrices(result, likelihood, verbose = verbose)


def _addMatrices(result, likelihood, normIntegralMatrix = None, acceptedNormIntegralMatrix = None, phaseSpaceIntegralVector = None, verbose = False):
	'''
	Calculate parameter covariance matrix and adds it and the given integral matrices to a new fit result
	@return: New fit result containing the covariance and integral matrices
	'''
	pars = []
	for parameter in likelihood.parameters():
		pars.append(result.fitParameter(parameter.parName()))

	# analytically calculate Hessian
	hessian = likelihood.Hessian(pars)
	# calculate and check eigenvalues
	eigenVectors = likelihood.HessianEigenVectors(hessian)
	if verbose:
		pyRootPwa.utils.printInfo("eigenvalues of (analytic) Hessian:")
	for i in xrange(len(eigenVectors)):
		if verbose:
			pyRootPwa.utils.printInfo("    {: .15e}".format(eigenVectors[i][1]))
		if eigenVectors[i][1] <= 0:
			pyRootPwa.utils.printWarn("eigenvalue {:d} of Hessian is not positive ({: .15e}).".format(i, eigenVectors[i][1]))

	covMatrix = likelihood.CovarianceMatrix(hessian)
	if verbose:
		pyRootPwa.utils.printInfo("(analytic) covariance matrix:")
		covMatrix.Print()

	oldCovMatrix = result.fitParCovMatrix()
	if oldCovMatrix.GetNcols() > 0 and oldCovMatrix.GetNrows() > 0:
		pyRootPwa.utils.printWarn("fit result from input already has a covariance matrix. it will be overwritten.")

	normIntegralMatrix = result.normIntegralMatrix() if normIntegralMatrix is None else normIntegralMatrix
	acceptedNormIntegralMatrix = result.acceptedNormIntegralMatrix() if acceptedNormIntegralMatrix is None else acceptedNormIntegralMatrix
	phaseSpaceIntegralVector = result.phaseSpaceIntegralVector() if phaseSpaceIntegralVector is None else phaseSpaceIntegralVector

	newResult = pyRootPwa.core.fitResult()
	newResult.fill(result,
	               covMatrix,
	               None,
	               normIntegralMatrix,
	               acceptedNormIntegralMatrix,
	               phaseSpaceIntegralVector)
	return newResult
