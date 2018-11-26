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


	def fit(self, nmbAttempts = 1):

		results = []

		for iAttempt in xrange(nmbAttempts):
			startParameters = self.startParameterGenerator()

			tStart = timeit.default_timer()
			results.append(self.optimize(startParameters))
			tDuration = timeit.default_timer() - tStart

			pyRootPwa.utils.printInfo("Minimization took "+str(tDuration)+" seconds. And "+str(results[iAttempt]['nmbEvals'])+" calls to objective (+gradient) function.")
			if results[iAttempt]['success']:
				pyRootPwa.utils.printSucc("Minimization successfull!")
			else:
				pyRootPwa.utils.printWarn("Minimization NOT successfull!")

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
			hessian = self.checkHessian(results[iAttempt])
			if self.storageLevel > 2 or (self.storageLevel == 1 and (iAttempt == iBestAttempt or iAttempt == iBestConvegedAttempt)):
				results[iAttempt]['hessian'] = hessian

		return results

	def checkHessian(self, result):
		tStart = timeit.default_timer()
		hessian = self.model.likelihood.hessianMatrixFitter(result['parameters'])
		tDuration = timeit.default_timer() - tStart
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


	def writeResultsRpwa(self, results, outputFileName, valTreeName = "pwa", valBranchName = "fitResult_v2"):

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


		normIntegralMatrix = buildIntegralMatrixFromSubmatrices(self.model.getNormSubmatrices())
		accIntegralMatrix = buildIntegralMatrixFromSubmatrices(self.model.getAccSubmatrices())
		normIntegrals = list(np.hstack(self.model.getNormIntegrals()))

		if self.model.waveNamesPosRefl and self.model.waveNamesNegRefl and self.model.rankPosRefl != self.model.rankNegRefl:
			pyRootPwa.utils.printWarn("Cannot store different ranks for positve and negative reflectivity in fit result. Using rank of positive reflectivity sector")
		if self.model.waveNamesPosRefl:
			rank = self.model.rankPosRefl
		else:
			rank = self.model.rankNegRefl

		for result in results:

			if result['hessian'] is not None:
				cov = np.linalg.inv(result['hessian'])
				fitparcovMatrix = pyRootPwa.ROOT.TMatrixD(cov.shape[0], cov.shape[0])
				for i in xrange(cov.shape[0]):
					for j in xrange(cov.shape[1]):
						fitparcovMatrix[i][j] = cov[i][j]
				hasHessian = True
				normIntegralMatrixResult = normIntegralMatrix
				accIntegralMatrixResult  = accIntegralMatrix
				normIntegralsResult = normIntegrals
			else:
				fitparcovMatrix  = None
				hasHessian   = False
				normIntegralMatrixResult = None
				accIntegralMatrixResult  = None
				normIntegralsResult = None
			fitResult.fill(
			               self.model.likelihood.nmbEvents,
			               1,
			               self.model.multibin.boundaries,
			               result['negLlhd'],
			               rank,
			               list(self.model.parameterMapping.paraFitter2AmpsForRpwaFitresult(result['parameters'])),
			               self.model.amplitudeNames(),
			               fitparcovMatrix,
			               self.model.parameterMapping.paraFitterCovMatrixIndicesForRpwaFitresult(),
			               normIntegralMatrixResult,
			               accIntegralMatrixResult,
			               normIntegralsResult,
			               result['hessianValid'],
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
				xtolRel = 1e-6,
				ftolAbs = 1e-8,
				maxeval = -1,
				vectorStorage = 100):
		'''
		@param xtolRel: Relative tolerance in the parameter space to stop the minimization
		@param ftolAbs: Absolute tolerance of the log-likelihood to assume to stop the minimization
		@param maxeval: Maximal number of evaluations (-1 means unlimited)
		@param vectorStorage: Number of gradients to keep for the approximation of the hessian matrix
		'''
		Fitter.__init__(self, model, checkLevel=checkLevel, storageLevel=storageLevel, startParameterGenerator=startValueGenerator)

		self.opt = nlopt.opt(algorithm, model.parameterMapping.nmbParameters)
		self.opt.set_xtol_rel(xtolRel)
		self.opt.set_ftol_abs(ftolAbs)
		self.opt.set_maxeval(maxeval)
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
