'''
@author: F. Kaspar, S. Wallner
'''
# pylint: disable=E1101

import inspect
import autograd.numpy as np


class ParameterMapping(object):
	'''
	The parameters transforms between three parameter spaces
		- The fitter parameters: array of real-valued parameters
		- The likelihood parameters: array of complex-valued parameters
			- The space in which autograd works
		- The parameters of the likelihood function (paraNegLlhd): List of parameter arrays, specific to the likelihood function
	'''

	def __init__(self, model):

		self.nmbParameters = None
		self.model = model

	def paraFitter2Llhd(self, paraFitter):
		'''
		@return: The paraLlhd form the paraFitter
		'''
		raise NotImplementedError("Needs to be implemented in specialized class!")

	def paraLlhd2Fitter(self, paraLlhd):
		'''
		@return: The fitter parameter from the likelihood parameters
		'''
		raise NotImplementedError("Needs to be implemented in specialized class!")

	def paraLlhd2negLlhd(self, paraLlhd):
		'''
		@return: Parameter list for the likelihood function
		'''
		raise NotImplementedError("Needs to be implemented in specialized class!")

	def paraNegLlhd2Llhd(self, paraNegLlhd):
		'''
		@return: The paraLlhd form the negative log-likelihood function parameters
		'''
		raise NotImplementedError("Needs to be implemented in specialized class!")

	def gradLlhd2Fitter(self, gradLlhd, gradFitter):
		'''
		@param gradLlhd: complex-valued array of derivatives in the likelihood parameter space from autgrad.grad
		@param gradLlhd: real-valued array of derivatives in the fitter parameter space
		'''
		raise NotImplementedError("Needs to be implemented in specialized class!")

	def hessianLlhd2Fitter(self, hessianLlhd, hessianFitter):
		'''
		@param hessianLlhd: nmbLlhdParameters*2 x nmbLlhdParameters*2 real-valued matrix
		@param hessianFitter: nmbParameters x nmbParameters real-valued matrix
		'''
		raise NotImplementedError("Needs to be implemented in specialized class!")




class ParameterMappingRpwa(ParameterMapping):
	'''
	classdocs
	'''

	def __init__(self, model):
		'''
		The fitter parameters are mapped onto the likelihood parameters in the following wave
		In the fitter parameters (paraFitter):
			The fitter parameters are real-valued.
			The first fitter parameters are the real- and imaginary parts of the transition amplitudes
				- first ordered by incoherent sector
				- second ordered by wave (first real, second imaginary part)
				- the anchor wave has only one parameter, its real part
			The last parameters are the real valued additional parameters of the likelihood
		In the likelihood parameters (paraLlhd):
			paraLlhd is one list of complex-valued parameters.
			The imaginary parts of the anchor waves are fixed to 0
			The additional parameters are real-valued, but in paraLlhd, they are complex-valued of an imaginary part of zero
		'''
		ParameterMapping.__init__(self, model)

		self.wavesInSectors = model.wavesInSectors
		self.nmbWavesInSectors = [ len(wavesInSector) for wavesInSector in self.wavesInSectors ]
		self.nmbSectors = len(self.nmbWavesInSectors)
		self.nmbAmplitudeParameters = np.sum([ 2*n-1 for n in self.nmbWavesInSectors])

		negLlhdFunction = model.clsLikelihood.negLlhd
		argsNegLlhdFunction = inspect.getargspec(negLlhdFunction).args
		self.nmbAdditionalParameters = len(argsNegLlhdFunction) -2 # the first one is self, the second one are the transition amplitudes
		self.nmbParameters = self.nmbAmplitudeParameters + self.nmbAdditionalParameters


		# build indices of paraLlhd where a new sector starts
		self.paraLlhdStartSectors = [0]
		for i, nWavesInSector in enumerate(self.nmbWavesInSectors):
			self.paraLlhdStartSectors.append(nWavesInSector + self.paraLlhdStartSectors[i])

		# build indices of paraFitter where a new sector starts
		self.paraFitterStartSectors = [0]
		for i, nWavesInSector in enumerate(self.nmbWavesInSectors):
			self.paraFitterStartSectors.append(2*nWavesInSector-1 + self.paraFitterStartSectors[i])

		self.nmbLlhdParameters = self.paraLlhdStartSectors[-1] + self.nmbAdditionalParameters


		# build index mapping between paraLlhd and paraFitter
		self.indicesRealFitterPara = [] # mapping paraLlhd index -> real value fitter parameter index
		self.indicesImagFitterPara = [] # mapping paraLlhd index -> imag value fitter parameter index
		for iSector, nWavesInSector in enumerate(self.nmbWavesInSectors):
			self.indicesRealFitterPara.append(2*0+self.paraFitterStartSectors[iSector])
			self.indicesImagFitterPara.append(None)
			for iWave in xrange(1, nWavesInSector):
				# -1 because the first have has only the real part
				self.indicesRealFitterPara.append(2*iWave-1+0+self.paraFitterStartSectors[iSector])
				self.indicesImagFitterPara.append(2*iWave-1+1+self.paraFitterStartSectors[iSector])
		for iAdditionalParameter in xrange(self.nmbAdditionalParameters):
			self.indicesRealFitterPara.append(iAdditionalParameter + self.nmbAmplitudeParameters)
			self.indicesImagFitterPara.append(None)
		# build list of indices of imaginary parts which are not None, i.e. which are real parameters
		self.indicesImagFitterParaNonNone = [i for i in self.indicesImagFitterPara if i is not None]
		self.indicesLlhdParaNonNoneImag = [j for j,i in enumerate(self.indicesImagFitterPara) if i is not None]

		# build parameter names
		self.paraNamesFitter = []
		for iSector, wavesInSector in enumerate(self.wavesInSectors):
			self.paraNamesFitter.append("V{s}_{w}_re".format(s=iSector, w=wavesInSector[0]))
			for wave in wavesInSector[1:]:
				self.paraNamesFitter.append("V{s}_{w}_re".format(s=iSector, w=wave))
				self.paraNamesFitter.append("V{s}_{w}_im".format(s=iSector, w=wave))
		for iAdditionalParameter in xrange(self.nmbAdditionalParameters):
			self.paraNamesFitter.append(argsNegLlhdFunction[iAdditionalParameter+2])


	def indexParaLlhd(self, iSector, iWave):
		'''
		@return: index of the paraLlhd for the given secor and wave
		'''
		return self.paraLlhdStartSectors[iSector] + iWave


	def paraFitter2Llhd(self, paraFitter):
		'''
		@return: The paraLlhd form the paraFitter
		'''
		paraLlhdImag = np.zeros(self.nmbLlhdParameters)
		paraLlhdImag[self.indicesLlhdParaNonNoneImag] = paraFitter[self.indicesImagFitterParaNonNone]
		paraLlhd = paraFitter[self.indicesRealFitterPara] + 1j*paraLlhdImag
		return paraLlhd

	def paraLlhd2Fitter(self, paraLlhd):
		paraFitter = np.empty(self.nmbParameters)
		paraFitter[self.indicesRealFitterPara] = np.real(paraLlhd)
		paraFitter[self.indicesImagFitterParaNonNone] = np.imag(paraLlhd[self.indicesLlhdParaNonNoneImag])
		return paraFitter


	def paraLlhd2negLlhd(self, paraLlhd):
		'''
		Split paraLlhd (one long parameter list) into [ [<list of T-vectors in sectors>], <additional parameters>]) for the `negLlhd` function.
		The complex-valued additional parameters in paraLlhd are mapped to real ones.
		'''
		amplitudes = tuple(paraLlhd[self.paraLlhdStartSectors[i]:self.paraLlhdStartSectors[i+1]] for i in range(self.nmbSectors))
		return [ amplitudes ] + list(np.real(paraLlhd[self.paraLlhdStartSectors[-1]:]))

	def paraNegLlhd2Llhd(self, paraNegLlhd):
		return np.hstack([np.hstack(paraNegLlhd[0])] + paraNegLlhd[1:])


	def gradLlhd2Fitter(self, gradLlhd, gradFitter):
		'''
		@param gradLlhd: complex-valued array of derivatives in the likelihood parameter space from autgrad.grad
		@param gradLlhd: real-valued array of derivatives in the fitter parameter space
		'''
		# complex conj: autograd requires this for complex valued parameters (complex valued parameters, but real optimization ;) )
		gradFitter[self.indicesRealFitterPara] = np.real(gradLlhd)
		gradFitter[self.indicesImagFitterParaNonNone] = -np.imag(gradLlhd[self.indicesLlhdParaNonNoneImag])


	def hessianLlhd2Fitter(self, hessianLlhd, hessianFitter):
		'''
		@param hessianLlhd: nmbLlhdParameters*2 x nmbLlhdParameters*2 real-valued matrix
		@param hessianFitter: nmbParameters x nmbParameters real-valued matrix
		'''
		indices = []
		for i,indexImag in enumerate(self.indicesImagFitterPara):
			indices.append(2*i)
			if indexImag is not None: # only take those imaginary parts that belong to a actual parameter
				indices.append(2*i+1)
		hessianFitter[:,:] = hessianLlhd[:,indices][indices]


	def paraFitter2AmpsForRpwaFitresult(self, paraFitter):
		return self.paraFitter2Llhd(paraFitter)[:self.nmbLlhdParameters - self.nmbAdditionalParameters]

	def paraFitterCovMatrixIndicesForRpwaFitresult(self):
		indices = []
		for i in xrange(self.nmbLlhdParameters - self.nmbAdditionalParameters):
			indices.append((self.indicesRealFitterPara[i], self.indicesImagFitterPara[i] if self.indicesImagFitterPara[i] is not None else -1))
		return indices
