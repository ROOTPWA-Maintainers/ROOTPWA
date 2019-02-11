'''
@author: F. Kaspar, S. Wallner
'''

import timeit
import sys
import os
# pylint: disable=E0611,E0401
import distutils.version
# pylint: enable=E0611,E0401
import nlopt
import numpy as np

import pyRootPwa


class Fitter(object):
	'''
	classdocs
	'''


	def __init__(self, model, checkLevel=1, storageLevel=1, startParameterGenerator=None):
		'''
		@param checkLevel: check level for hessian matrix checks 0 no checks, 1: check only best result, 2: check all results
		@param storageLevel: Storage level of matrices 0: store nothing, 1: store hessian and integral matrix only of best result, 2: store all matrices of all results
		'''
		self.model = model
		self.checkLevel = checkLevel
		self.storageLevel = storageLevel
		self.startParameterGenerator = startParameterGenerator


	def fit(self, nmbAttempts = 1, verbosity = 1):
		'''
		@param verbosity: -1: Nothing, 0: Only summary 1: Progress bar, 2: Information from each attempt
		'''

		results = []

		progress = None
		if verbosity == 1:
			progress = pyRootPwa.utils.progressBar(maximum=nmbAttempts)
			progress.start()

		for iAttempt in xrange(nmbAttempts):
			startParameters = self.startParameterGenerator()

			tStart = timeit.default_timer()
			results.append(self.optimize(startParameters))
			tDuration = timeit.default_timer() - tStart

			if verbosity > 1:
				if results[iAttempt]['success']:
					pyRootPwa.utils.printSucc("Minimization successful (-log(L) = {0:.2f})! ".format(results[iAttempt]["negLlhd"]) + self.getStatus(results[iAttempt]['status']))
				else:
					pyRootPwa.utils.printWarn("Minimization NOT successfull! " + self.getStatus(results[iAttempt]['status']))
				pyRootPwa.utils.printInfo("Minimization took "+str(tDuration)+" seconds. And "+str(results[iAttempt]['nmbEvals'])+" calls to objective (+gradient) function.")

			if progress:
				progress.update(iAttempt)

		iBestAttempt         = np.argmin([r['negLlhd'] for r in results])
		iBestConvegedAttempt = np.argmin([r['negLlhd'] if r['success'] else 1e300 for r in results])

		# check results
		if self.storageLevel > 2 or self.checkLevel > 2:
			attemptsToCeck = xrange(nmbAttempts)
		elif self.storageLevel == 1 or self.checkLevel == 1:
			attemptsToCeck = sorted(list(set([iBestAttempt, iBestConvegedAttempt])))
		else:
			attemptsToCeck = []
		for iAttempt in attemptsToCeck:
			hessian = self.checkHessian(results[iAttempt], verbosity)
			if self.storageLevel >= 2 or (self.storageLevel == 1 and (iAttempt == iBestAttempt or iAttempt == iBestConvegedAttempt)):
				results[iAttempt]['hessian'] = hessian

		if verbosity >= 0:
			nmbSuccess = len([r for r in results if r["success"]])
			pyRootPwa.utils.printInfo("Finished {0} fit attempts.".format(nmbAttempts))
			if nmbSuccess > 0:
				pyRootPwa.utils.printSucc("{0} fit attempts successful.".format(nmbSuccess))
			if nmbAttempts - nmbSuccess > 0:
				pyRootPwa.utils.printWarn("{0} fit attempts not successful.".format(nmbAttempts - nmbSuccess))
			nmbValidHessian = len([r for r in results if r["hessianValid"] is True])
			nmbInvalidHessian = len([r for r in results if r["hessianValid"] is False])
			if nmbValidHessian > 0:
				pyRootPwa.utils.printSucc("{0} fit attempts with valid Hessian matrix.".format(nmbValidHessian))
			if nmbInvalidHessian > 0:
				pyRootPwa.utils.printWarn("{0} fit attempts with invalid Hessian matrix.".format(nmbInvalidHessian))

		return results

	def checkHessian(self, result, verbosity):
		tStart = timeit.default_timer()
		hessian = self.model.likelihood.hessianMatrixFitter(result['parameters'])
		tDuration = timeit.default_timer() - tStart
		if verbosity > 1:
			pyRootPwa.utils.printInfo("Calculation of Hessian took "+str(tDuration)+" seconds.")

		eigenvalues, _ = np.linalg.eig(hessian)
		result['hessianValid'] = True
		for val in eigenvalues:
			if val < 0:
				result['hessianValid'] = False
		return hessian

	def optimize(self, startParameters):
		'''
		@return: result dictionary. Keys are 'parameters', 'startParameters', 'negLlhd', 'success', 'status', 'nmbEvals', 'hessian', 'hessianValid'
		'''
		raise NotImplementedError("This method must be implemented in the derived classes")

	def getStatus(self, statuscode):
		'''
		Get fit status information based on the status code of the fitter
		'''
		raise NotImplementedError("This method must be implemented in the derived classes")


def writeResultsRpwa(model, results, outputFileName, integralsStorageLevel, valTreeName = "pwa", valBranchName = "fitResult_v2"):
	'''
	@param integralsStorageLevel: 0=do not store integrals, 1=store integrals of best and best converged result, 2=store integrals for all results
	'''

	if os.path.exists(outputFileName) or os.path.exists(outputFileName+'.root'):
		pyRootPwa.utils.printErr("Output file already exists! Aborting ...")
		sys.exit(1)
	outputFile = pyRootPwa.ROOT.TFile.Open(outputFileName, "UPDATE")
	if (not outputFile) or outputFile.IsZombie():
		pyRootPwa.utils.printErr("cannot open output file '" + outputFileName + "'. Aborting...")
		sys.exit(1)
	fitResult = pyRootPwa.core.fitResult()
	tree = outputFile.Get(valTreeName)
	if not tree:
		pyRootPwa.utils.printInfo("file '" + outputFileName + "' is empty. "
				+ "creating new tree '" + valTreeName + "' for PWA result.")
		tree = pyRootPwa.ROOT.TTree(valTreeName, valTreeName)
		if not fitResult.branch(tree, valBranchName):
			pyRootPwa.utils.printErr("failed to create new branch '" + valBranchName + "' in file '" + outputFileName + "'.")
			sys.exit(1)
	else:
		fitResult.setBranchAddress(tree, valBranchName)


	normIntegralMatrix = buildIntegralMatrixFromSubmatrices(model.getNormSubmatrices())
	accIntegralMatrix = buildIntegralMatrixFromSubmatrices(model.getAccSubmatrices())
	normIntegrals = list(np.hstack(model.getNormIntegrals()))

	if model.waveNamesPosRefl and model.waveNamesNegRefl and model.rankPosRefl != model.rankNegRefl:
		pyRootPwa.utils.printWarn("Cannot store different ranks for positve and negative reflectivity in fit result. Using rank of positive reflectivity sector")
	if model.waveNamesPosRefl:
		rank = model.rankPosRefl
	else:
		rank = model.rankNegRefl

	iBestAttempt         = np.argmin([r['negLlhd'] for r in results])
	iBestConvegedAttempt = np.argmin([r['negLlhd'] if r['success'] else 1e300 for r in results])

	for iResult, result in enumerate(results):

		if result['hessian'] is not None:
			cov = np.linalg.inv(result['hessian'])
			fitparcovMatrix = pyRootPwa.ROOT.TMatrixD(cov.shape[0], cov.shape[0])
			for i in xrange(cov.shape[0]):
				for j in xrange(cov.shape[1]):
					fitparcovMatrix[i][j] = cov[i][j]
			hasHessian = True
		else:
			fitparcovMatrix  = None
			hasHessian   = False

		if integralsStorageLevel >= 2 or (integralsStorageLevel == 1 and (iResult == iBestAttempt or iResult == iBestConvegedAttempt)):
			normIntegralMatrixResult = normIntegralMatrix
			accIntegralMatrixResult  = accIntegralMatrix
			normIntegralsResult = normIntegrals
		else:
			normIntegralMatrixResult = None
			accIntegralMatrixResult  = None
			normIntegralsResult = None

		fitResult.fill(
						model.likelihood.nmbEvents,
						1,
						model.multibin.boundaries,
						result['negLlhd'],
						rank,
						list(model.parameterMapping.paraFitter2AmpsForRpwaFitresult(result['parameters'])),
						model.amplitudeNames(),
						fitparcovMatrix,
						model.parameterMapping.paraFitterCovMatrixIndicesForRpwaFitresult(),
						normIntegralMatrixResult,
						accIntegralMatrixResult,
						normIntegralsResult,
						result['success'],
						hasHessian
			)
		tree.Fill()
	nmbBytes = tree.Write()
	outputFile.Close()
	if nmbBytes == 0:
		pyRootPwa.utils.printErr("problems writing fit result to TKey 'fitResult' "
				+ "in file '" + outputFileName + "'")
		sys.exit(1)
	else:
		pyRootPwa.utils.printSucc("wrote fit result to TKey 'fitResult' "
				+ "in file '" + outputFileName + "'")


class NLoptFitter(Fitter):

	def __init__(self, model, checkLevel=1, storageLevel=1, startValueGenerator=None,
				algorithm = nlopt.LD_LBFGS,
				xtolRel = 1e-4,
				ftolAbs = 1e-6,
				maxeval = 50000,
				vectorStorage = None):
		'''
		@param xtolRel: Relative tolerance in the parameter space to stop the minimization
		@param ftolAbs: Absolute tolerance of the log-likelihood to assume to stop the minimization
		@param maxeval: Maximal number of evaluations (-1 means unlimited)
		@param vectorStorage: Number of gradients to keep for the approximation of the hessian matrix
		'''
		Fitter.__init__(self, model, checkLevel=checkLevel, storageLevel=storageLevel, startParameterGenerator=startValueGenerator)

		self.opt = nlopt.opt(algorithm, model.parameterMapping.nmbParameters)
		self.opt.set_maxeval(maxeval)
		self.opt.set_xtol_rel(xtolRel)
		self.opt.set_ftol_abs(ftolAbs)
		if vectorStorage is not None:
			self.opt.set_vector_storage(vectorStorage)
		self.opt.set_min_objective(self.model.likelihood.f)


	def optimize(self, startParameters):
		optimizedParameters = self.opt.optimize(startParameters)
		result = {}
		result['status'] = self.opt.last_optimize_result()
		result['success'] = self.opt.last_optimize_result() > 0
		result['parameters'] = optimizedParameters
		result['startParameters'] = startParameters
		result['negLlhd']  = self.opt.last_optimum_value()
		result['hessian'] = None
		result['hessianValid'] = None
# pylint: disable=E1101
		if distutils.version.LooseVersion(nlopt.__version__) >= distutils.version.LooseVersion("2.5.0"):
# pylint: enable=E1101
			result['nmbEvals'] = self.opt.get_numevals()
		else:
			result['nmbEvals'] = None
		return result


	def getStatus(self, statuscode):
		if statuscode == 1:
			return ":)"
		elif statuscode == 2:
			return "Stop-value reached."
		elif statuscode == 3:
			return "Relative/absolute function tolerance reached."
		elif statuscode == 4:
			return "Relative/absolute parameter tolerance reached."
		elif statuscode == 5:
			return "Maximal number of evaluations reached."
		elif statuscode == 6:
			return "Maximal evaluation time reached."
		elif statuscode == -1:
			return "Generic failure."
		elif statuscode == -2:
			return "Invalid arguments."
		elif statuscode == -3:
			return "Ran out of memory."
		elif statuscode == -4:
			return "Roundoff errors limited progress."
		elif statuscode == -5:
			return "Forced termination."
		else:
			raise ValueError("Status code '{0}' not implemented!".format(statuscode))



def buildIntegralMatrixFromSubmatrices(submatrices):
	'''
	Build one `complexMatrix` matrix for all given submatrices.
	The block off-diagonal elements will be set to 0.
	@param submatrices: List of submatrices, each matrix is a 2D complex-valued numpy array
	'''
	nmbWaves = np.sum([m.shape[0] for m in submatrices])
	totMatrix = pyRootPwa.core.complexMatrix(nmbWaves, nmbWaves)
	for i in xrange(nmbWaves):
		for j in xrange(nmbWaves):
			totMatrix.set(i,j,0.0)

	offset = 0
	for submatrix in submatrices:
		for i in xrange(submatrix.shape[0]):
			for j in xrange(submatrix.shape[1]):
				totMatrix.set(i+offset, j+offset, submatrix[i,j])
		offset += submatrix.shape[0]
	return totMatrix
