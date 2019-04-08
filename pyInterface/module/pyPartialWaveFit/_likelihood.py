'''
@author: F. Kaspar, S. Wallner
'''
# pylint: disable=E1120,E1101,W0221

import sys
import inspect

import autograd
import autograd.numpy as np
import pyRootPwa

class Likelihood(object):
	'''
	classdocs
	'''


	def __init__(self, decayAmplitudes, accMatrices, normMatrices, normIntegrals, parameterMapping):
		'''
		@param decayAmplitudes: List of lists of complex-valued matrix of normalized decay amplitudes, one for each data-set and incoherent sector
								The first index is the wave number, the second index is the event number
								decayAmplitudes[<dataset>][<sector>][<wave>,<event>]
		@param accMatrices: List of lists of integral matrices with acceptance effects, one for each data-set and incoherent sector
		                    accMatrices[<dataset>][<sector>][<waveA>,<waveB>]
		@param normMatrices: List of lists of integral matrices without acceptance effects, one for each data-set and for each incoherent sector
		                     normMatrices[<dataset>][<sector>][<waveA>,<waveB>]
		'''

		if not decayAmplitudes:
			raise Exception("No decay amplitudes given")

		if len(decayAmplitudes) != len(accMatrices) or len(decayAmplitudes) != len(normMatrices):
			raise Exception('number of dataset not matching!')

		self.nmbDatasets = len(decayAmplitudes)
		for iDataset in range(self.nmbDatasets):
			if len(decayAmplitudes[iDataset]) != len(accMatrices[iDataset]) or len(decayAmplitudes[iDataset]) != len(normMatrices[iDataset]):
				raise Exception('Number of sectors not matching!')


		self.nmbSectors = len(decayAmplitudes[0])
		self.nmbEvents = np.array([ amps[0].shape[1] for amps in decayAmplitudes ])
		self.nmbWavesInSectors = []
		for iDataset in range(self.nmbDatasets):
			for i, amplitudes in enumerate(decayAmplitudes[iDataset]):
				if self.nmbEvents[iDataset] != amplitudes.shape[1] and amplitudes.shape[1] != 1:
					raise Exception('Events not matching')

				if amplitudes.shape[0] !=  accMatrices[iDataset][i].shape[0] or amplitudes.shape[0] !=  normMatrices[iDataset][i].shape[0] \
				or amplitudes.shape[0] !=  accMatrices[iDataset][i].shape[1] or amplitudes.shape[0] !=  normMatrices[iDataset][i].shape[1]:
					raise Exception('Matrices not matching')

				if iDataset == 0:
					self.nmbWavesInSectors.append(amplitudes.shape[0])
				else:
					if self.nmbWavesInSectors[i] != amplitudes.shape[0]:
						raise Exception('Data sets have different number of waves')

		self.decayAmplitudes = decayAmplitudes
		self.accMatrices = accMatrices
		self.normMatrices = normMatrices
		self.normIntegrals = normIntegrals


		self.parameterMapping = parameterMapping

		self._valueAndGrad = autograd.value_and_grad(self._negLlhd)
		self._grad = autograd.grad(self._negLlhd)
		self._hessian = autograd.jacobian(self._gradientForHessian)
		self.countFdF = 0 # number of calls to f with gradient
		self.countF = 0 # number of calls to f without gradient


	def setParameters(self):
		pass # no parameters for this likelihood class


	def negLlhd(self, transitionAmps, datasetRatioParameters):
		'''
		@return Negative log-likelihood for the given set of transition amplitudes
		'''
		datasetRatios, datasetRatiosLog = self.parameterMapping.datasetRatioParameters2DatasetRatios(datasetRatioParameters)
		negLlhd = 0.0
		for iDataset in xrange(self.nmbDatasets):
			nBar = 0.0
			dataTerm = 0.0
			for iSector in xrange(self.nmbSectors):
				nBar = nBar + datasetRatios[iDataset]*np.real(np.dot( np.dot( self.accMatrices[iDataset][iSector],np.conj(transitionAmps[iSector]) ), transitionAmps[iSector]))
				sumOfAmps = np.dot(transitionAmps[iSector],self.decayAmplitudes[iDataset][iSector])
				dataTerm = dataTerm + np.real(sumOfAmps * np.conj(sumOfAmps))

			# likelihood calculation
			negLlhd += -np.sum( np.log( dataTerm ) ) + nBar - self.nmbEvents[iDataset]*datasetRatiosLog[iDataset]
		return negLlhd


# pylint: disable=C0103,E1102
	def f(self, paraFitter, gradFitter):
		paraLlhd = self.parameterMapping.paraFitter2Llhd(paraFitter)
		if gradFitter.size > 0:
			self.countFdF += 1
			negLL, gradLlhd = self._valueAndGrad(paraLlhd)
			self.parameterMapping.gradLlhd2Fitter(gradLlhd, gradFitter)
		else:
			self.countF += 1
			negLL = self._negLlhd(paraLlhd)
		return negLL

	def _negLlhd(self, paraLlhd):
		'''
		Wrapper function to separate transition amplitudes from other fit parameters
		'''
		return self.negLlhd(*self.parameterMapping.paraLlhd2negLlhd(paraLlhd))


	def _gradientForHessian(self, paraLlhd):
		'''
		Split complex-valued gradient by real and imaginary parts
		in order to put this into the jacobian for the hessian calculation
		'''
		gradient = self._grad(paraLlhd)
		return np.ravel(np.column_stack((np.real(gradient), - np.imag(gradient))))

	def hessian(self, paraLlhd):
		'''
		@return: Hessian matrix of the complex-valued likelihood parameters
				 2*nmbLlhdParameter x 2*nmbLlhdParameters matrix of re_i, imag_i as components
		'''
		hessianMatrix = np.conj(self._hessian(paraLlhd)).view(np.float64)
		return hessianMatrix

	def hessianMatrixFitter(self, paraFitter):
		'''
		@return: Hessian matrix of the real-valued fitter parameters
				 nmbParameter x nmbParameter matrix
		'''
		hessianMatrix = self.hessian(self.parameterMapping.paraFitter2Llhd(paraFitter))
		hessianMatrixFitterParameter = np.empty((self.parameterMapping.nmbParameters, self.parameterMapping.nmbParameters))
		self.parameterMapping.hessianLlhd2Fitter(hessianMatrix, hessianMatrixFitterParameter)
		return hessianMatrixFitterParameter


def getLikelihoodClassNames():
	likelihoods = []
	thisModule = sys.modules[__name__]
	for name, _ in inspect.getmembers(thisModule, inspect.isclass):
		if name.startswith('Likelihood'):
			likelihoods.append(name)
	return likelihoods



class LikelihoodCauchy(Likelihood):
	'''
	Likelihood with cauchy regularization term
	'''
	def __init__(self, decayAmplitudes, accMatrices, normMatrices, normIntegrals, parameterMapping):
		Likelihood.__init__(self, decayAmplitudes, accMatrices, normMatrices, normIntegrals, parameterMapping)
		self.sectors = range(len(decayAmplitudes))[:-1] # apply regularization to all but the last one, which corresponds to the flat wave usually
		self.width = [0.5]*len(self.decayAmplitudes)


	def setParameters(self, width, correctAcceptance=False, sectors= None):
		if isinstance(width, float):
			if not correctAcceptance:
				self.width = [width] * len(self.decayAmplitudes)
			else:
				self.width = [width / np.sqrt(np.abs(accMatrix.diagonal())) for accMatrix in self.accMatrices[0]]
		else:
			self.width = width
			if correctAcceptance:
				pyRootPwa.utils.printErr("'widthTimesAcceptance' can only be used with scalar width!")
				raise Exception()

		if sectors is not None:
			if True in [ i < 0 or i >= self.nmbSectors for i in sectors]:
				pyRootPwa.utils.printErr("One of the sectors for the Cauchy regularization is out of range")
				raise Exception()
			self.sectors = sectors


	def negLlhd(self, transitionAmps, datasetRatios):
		negLlhd = Likelihood.negLlhd(self, transitionAmps, datasetRatios)
		for i in self.sectors:
			absT = np.abs(transitionAmps[i])
			negLlhd = negLlhd - np.sum(np.log(1.0/(1.0+absT**2/self.width[i]**2)))
		return negLlhd


class LikelihoodConnected(object):


	def __init__(self, likelihoods, binWidths):

		self.likelihoods = likelihoods
		self.binWidthsNormalization = np.empty((len(likelihoods), likelihoods[0].parameterMapping.nmbLlhdParameters-likelihoods[0].parameterMapping.nmbDatasets))
		for i in range(self.binWidthsNormalization.shape[0]):
			self.binWidthsNormalization[i,:] = np.sqrt(1.0/binWidths[i]*0.02)

		self.countFdF = 0  # number of calls to f with gradient
		self.countF = 0  # number of calls to f without gradient

		# this is a bad hack :(
		self.nmbEvents = np.sum([likelihood.nmbEvents for likelihood in self.likelihoods])
		self.nmbParameters = np.sum([likelihood.parameterMapping.nmbParameters for likelihood in self.likelihoods])

		self.parameterMapping = None

		self.strength = 0.08

		self._connectionValueAndGrad = autograd.value_and_grad(self._connection)
		self._connectionGrad = autograd.grad(self._connection)
		self._connectionHessian = autograd.jacobian(self._gradientForConnectionHessian)


	# make this customizable
	def _connection(self, paraLlhd):
		paras = np.dstack([np.hstack(self.likelihoods[i].parameterMapping.paraLlhd2negLlhd(self.parameterMapping.paraLlhdOfBin(paraLlhd,i))[0]) for i in xrange(len(self.likelihoods))])
		paras = paras.reshape((paras.shape[1], paras.shape[2])).transpose()
		# normalize to 20 MeV bins
		parasNormed = paras*self.binWidthsNormalization
		return np.sum(self.strength * (np.abs(parasNormed[1:] - parasNormed[:-1]) ** 2))
# 		return np.sum(self.strength * (np.abs(paras[1:-1] - paras[:-2]) ** 2 + np.abs(paras[1:-1] - paras[2:]) ** 2))
		# return self.strength * np.log(1. + np.sum( np.abs(paras[1:-1]-paras[:-2])**2 + np.abs(paras[1:-1]-paras[2:])**2) )
		# return self.strength * np.sum( np.sqrt(0.01+np.abs(paras[1:-1]-paras[:-2])**2) + np.sqrt(0.01+np.abs(paras[1:-1]-paras[2:])**2) )


	def setParameters(self, strength=0.08, scaleStrengthByEvents=False, **kwargs):
		'''
		@param strength: Strength parameter of the connection term
		All other keyword arguments are passed to the individual likelihoods
		'''
		self.strength = strength
		if scaleStrengthByEvents:
			if not isinstance(strength, float):
				raise Exception()
			strengths = np.empty((len(self.likelihoods),np.sum(self.likelihoods[0].nmbWavesInSectors)))
			for iLikelihood, likelihood in enumerate(self.likelihoods):
				nmbEventAccCorr = likelihood.nmbEvents/np.abs(likelihood.accMatrices[0][-1][0,0])
				strengths[iLikelihood,:] = strength/nmbEventAccCorr
			self.strength =0.5*(strengths[1:,:]+strengths[0:-1,:])
		else:
			self.strength = strength

		if kwargs:
			for likelihood in self.likelihoods:
				likelihood.setParameters(**kwargs)


	def _gradientForConnectionHessian(self, paraLlhd):
		'''
		Split complex-valued gradient by real and imaginary parts
		in order to put this into the jacobian for the hessian calculation
		'''
		gradient = self._connectionGrad(paraLlhd)
		return np.ravel(np.column_stack((np.real(gradient), -np.imag(gradient))))


	def connectionHessian(self, paraLlhd):
		'''
		@return: Hessian matrix of the complex-valued likelihood parameters
				 2*nmbLlhdParameter x 2*nmbLlhdParameters matrix of re_i, imag_i as components
		'''
		hessianMatrix = np.conj(self._connectionHessian(paraLlhd)).view(np.float64)
		return hessianMatrix


	def connectionHessianMatrixFitter(self, paraFitter):
		'''
		@return: Hessian matrix of the real-valued fitter parameters
				 nmbParameter x nmbParameter matrix
		'''
		hessianMatrix = self.connectionHessian(self.parameterMapping.paraFitter2Llhd(paraFitter))
		hessianMatrixFitterParameter = np.empty((self.parameterMapping.nmbParameters, self.parameterMapping.nmbParameters))
		self.parameterMapping.hessianLlhd2Fitter(hessianMatrix, hessianMatrixFitterParameter)
		return hessianMatrixFitterParameter


# pylint: disable=C0103,E1102
	def f(self, paraFitter, gradFitter):

		negLL = 0

		paraLlhd = self.parameterMapping.paraFitter2Llhd(paraFitter)

		if gradFitter.size > 0:
			self.countFdF += 1
			for iLikelihood, likelihood in enumerate(self.likelihoods):
				paras = self.parameterMapping.paraFitterOfBin(paraFitter, iLikelihood)
				negLL += likelihood.f(paras, gradFitter[self.parameterMapping.offsetsFitterForBins[iLikelihood]:self.parameterMapping.offsetsFitterForBins[iLikelihood+1]])

			negPrior, gradL = self._connectionValueAndGrad(paraLlhd)
			negLL += negPrior
			gradFitterPrior = np.zeros(len(gradFitter), dtype=np.float64)
			self.parameterMapping.gradLlhd2Fitter(gradL, gradFitterPrior)

			gradFitter += gradFitterPrior

		else:
			self.countF += 1

			for iLikelihood, likelihood in enumerate(self.likelihoods):
				paras = self.parameterMapping.paraFitterOfBin(paraFitter, iLikelihood)
				negLL += likelihood.f(paras, gradFitter)

			negLL += self._connection(paraLlhd)

		return negLL


	def hessianMatrixFitter(self, paraFitter):

		hessianMatrixFitterParameter = np.zeros((self.nmbParameters, self.nmbParameters))
		for iLikelihood, likelihood in enumerate(self.likelihoods):
			paras = self.parameterMapping.paraFitterOfBin(paraFitter, iLikelihood)
			hessianMatrixFitterParameter[self.parameterMapping.offsetsFitterForBins[iLikelihood]:self.parameterMapping.offsetsFitterForBins[iLikelihood+1],
			                             self.parameterMapping.offsetsFitterForBins[iLikelihood]:self.parameterMapping.offsetsFitterForBins[iLikelihood+1]]\
			  = likelihood.hessianMatrixFitter(paras)

		hessianMatrixFitterParameter += self.connectionHessianMatrixFitter(paraFitter)
		return hessianMatrixFitterParameter
