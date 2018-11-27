'''
@author: F. Kaspar, S. Wallner
'''
# pylint: disable=E1120,E1101

import autograd
import autograd.numpy as np

class Likelihood(object):
	'''
	classdocs
	'''


	def __init__(self, decayAmplitudes, accMatrices, normMatrices, normIntegrals, parameterMapping):
		'''
		@param decayAmplitudes: List of complex-valued matrix of normalized decay amplitudes, one for each incoherent sector
                                The first index is the wave number, the second index is the event number
		@param accMatrices: List of integral matrices with acceptance effects, one for each incoherent sector
		@param normMatrices: List of integral matrices without acceptance effects, one for each incoherent sector
		'''

		if not decayAmplitudes:
			raise Exception("No decay amplitudes given")

		if len(decayAmplitudes) != len(accMatrices) or len(decayAmplitudes) != len(normMatrices):
			raise Exception('Number of sectors not matching!')

		self.nmbSectors = len(decayAmplitudes)
		self.nmbEvents = decayAmplitudes[0].shape[1]
		self.nmbWavesInSectors = []
		for i, amplitudes in enumerate(decayAmplitudes):
			if self.nmbEvents != amplitudes.shape[1] and amplitudes.shape[1] != 1:
				raise Exception('Events not matching')

			if amplitudes.shape[0] !=  accMatrices[i].shape[0] or amplitudes.shape[0] !=  normMatrices[i].shape[0] \
			   or amplitudes.shape[0] !=  accMatrices[i].shape[1] or amplitudes.shape[0] !=  normMatrices[i].shape[1]:
				raise Exception('Matrices not matching')

			self.nmbWavesInSectors.append(amplitudes.shape[0])

		self.decayAmplitudes = decayAmplitudes
		self.accMatrices = accMatrices
		self.normMatrices = normMatrices
		self.normIntegrals = normIntegrals

		self.parameterMapping = parameterMapping

		self._valueAndGrad = autograd.value_and_grad(self._negLlhd)
		self._grad = autograd.grad(self._negLlhd)
		self._hessian = autograd.jacobian(self._gradientForHessian)


	def negLlhd(self, transitionAmps):
		'''
		@return Negative log-likelihood for the given set of transition amplitudes
		'''

		nBar = 0.0
		dataTerm = 0.0

		for i in xrange(self.nmbSectors):
			nBar = nBar + np.real(np.dot( np.dot( self.accMatrices[i],np.conj(transitionAmps[i]) ), transitionAmps[i]))
			sumOfAmps = np.dot(transitionAmps[i],self.decayAmplitudes[i])
			dataTerm = dataTerm + np.real(sumOfAmps * np.conj(sumOfAmps))

		# likelihood calculation
		return -np.sum( np.log( dataTerm ) ) + nBar


# pylint: disable=C0103,E1102
	def f(self, paraFitter, gradFitter):
		paraLlhd = self.parameterMapping.paraFitter2Llhd(paraFitter)
		if gradFitter.size > 0:
			negLL, gradLlhd = self._valueAndGrad(paraLlhd)
			self.parameterMapping.gradLlhd2Fitter(gradLlhd, gradFitter)
		else:
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
		return np.ravel(np.column_stack((np.real(gradient), -np.imag(gradient))))

	def hessian(self, paraLlhd):
		'''
		@return: Hessian matrix of the complex-valued likelihood parameters
		         2*nmbLlhdParameter x 2*nmbLlhdParameters matrix of re_i, imag_i as components
		'''
		hessianMatrix = self._hessian(paraLlhd).view(np.float64)
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
