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

	def __init__(self, model, realWaves = None, zeroWaves = None):
		'''
		@param realWaves: List lists of waves whose amplitude should be fixed to real, i.e. the imaginary part is fixed to zero. One list for each incoherent sector.
		@param zeroWaves: List of waves whose amplitude should be fixed to zero. Fixing of the reference wave is ignored.
		'''
		ParameterMapping.__init__(self, model)

		self.wavesInSectors = model.wavesInSectors
		self.nmbWavesInSectors = [ len(wavesInSector) for wavesInSector in self.wavesInSectors ]
		self.nmbSectors = len(self.nmbWavesInSectors)
		self.nmbDatasets = model.nmbDatasets

		if realWaves is None:
			realWaves = [ [] ] * self.nmbSectors

		if len(realWaves) != self.nmbSectors:
			raise ValueError("List of real-valued wave does not have n sectors.")

		if zeroWaves is None:
			zeroWaves = []

		negLlhdFunction = model.clsLikelihood.negLlhd
		argsNegLlhdFunction = inspect.getargspec(negLlhdFunction).args
		self.nmbAdditionalParameters = len(argsNegLlhdFunction) -3 # the first one is self, the second one are the transition amplitudes, the third one are the data-set ratios
		self.nmbAmplitudeParameters = np.sum([ 2*n-len(realWaves[i]) for i,n in enumerate(self.nmbWavesInSectors)]) - 2*len(zeroWaves)
		self.nmbDatasetRatioParameters = self.nmbDatasets-1
		self.nmbParameters = self.nmbAmplitudeParameters + self.nmbDatasetRatioParameters + self.nmbAdditionalParameters


		# build indices of paraLlhd where a new sector starts
		self.paraLlhdStartSectors = [0]
		for i, nWavesInSector in enumerate(self.nmbWavesInSectors):
			self.paraLlhdStartSectors.append(nWavesInSector + self.paraLlhdStartSectors[i])
		self.paraLlhdStartDatasetRatios = self.paraLlhdStartSectors[-1]
		self.paraLlhdStartStartAdditionParameters = self.paraLlhdStartDatasetRatios + self.nmbDatasets

		self.nmbLlhdParameters = self.paraLlhdStartSectors[-1] + self.nmbDatasets + self.nmbAdditionalParameters

		# build index mapping between paraLlhd and paraFitter
		iFitterPara = 0
		self.indicesRealFitterPara = [] # mapping paraLlhd index -> real value fitter parameter index
		self.indicesImagFitterPara = [] # mapping paraLlhd index -> imag value fitter parameter index
		for iSector, nWavesInSector in enumerate(self.nmbWavesInSectors):
			for iWave in xrange( nWavesInSector):
				if self.wavesInSectors[iSector][iWave] not in zeroWaves:
					self.indicesRealFitterPara.append(iFitterPara)
					iFitterPara += 1
					if self.wavesInSectors[iSector][iWave] not in realWaves[iSector]:
						self.indicesImagFitterPara.append(iFitterPara)
						iFitterPara += 1
					else:
						self.indicesImagFitterPara.append(None)
				else:
					self.indicesRealFitterPara.append(None)
					self.indicesImagFitterPara.append(None)
		for iDataset in xrange(self.nmbDatasets):
			if iDataset == 0:
				self.indicesRealFitterPara.append(None)
			else:
				self.indicesRealFitterPara.append(iFitterPara)
				iFitterPara += 1
			self.indicesImagFitterPara.append(None)
		for iAdditionalParameter in xrange(self.nmbAdditionalParameters):
			self.indicesRealFitterPara.append(iFitterPara)
			iFitterPara += 1
			self.indicesImagFitterPara.append(None)
		# build list of indices which are not None, i.e. which are real parameters
		self.indicesRealFitterParaNonNone = np.array([i for i in self.indicesRealFitterPara if i is not None], dtype=np.int64)
		self.indicesLlhdParaNonNoneReal   = np.array([j for j,i in enumerate(self.indicesRealFitterPara) if i is not None], dtype=np.int64)
		self.indicesImagFitterParaNonNone = np.array([i for i in self.indicesImagFitterPara if i is not None], dtype=np.int64)
		self.indicesLlhdParaNonNoneImag   = np.array([j for j,i in enumerate(self.indicesImagFitterPara) if i is not None], dtype=np.int64)

		assert iFitterPara == self.nmbParameters

		# build parameter names
		self.paraNamesFitter = []
		for iSector, wavesInSector in enumerate(self.wavesInSectors):
			for wave in wavesInSector:
				if wave not in zeroWaves:
					self.paraNamesFitter.append("V{s}_{w}_re".format(s=iSector, w=wave))
					if wave not in realWaves[iSector]:
						self.paraNamesFitter.append("V{s}_{w}_im".format(s=iSector, w=wave))
		for iDatasetRatioParameter in xrange(self.nmbDatasetRatioParameters):
			self.paraNamesFitter.append("r_{0}".format(iDatasetRatioParameter+1))
		for iAdditionalParameter in xrange(self.nmbAdditionalParameters):
			self.paraNamesFitter.append(argsNegLlhdFunction[iAdditionalParameter+2])


	def indexParaLlhd(self, iSector, iWave):
		'''
		@return: index of the paraLlhd for the given secor and wave
		'''
		return self.paraLlhdStartSectors[iSector] + iWave

	def indexParaLlhdDatasetRatio(self, iDataset):
		return self.paraLlhdStartDatasetRatios + iDataset


	def paraFitter2Llhd(self, paraFitter):
		'''
		@return: The paraLlhd form the paraFitter
		'''
		paraLlhdReal = np.zeros(self.nmbLlhdParameters)
		paraLlhdReal[self.indicesLlhdParaNonNoneReal] = paraFitter[self.indicesRealFitterParaNonNone]
		paraLlhdImag = np.zeros(self.nmbLlhdParameters)
		paraLlhdImag[self.indicesLlhdParaNonNoneImag] = paraFitter[self.indicesImagFitterParaNonNone]
		# set data-set ratio reference
		paraLlhdReal[self.indexParaLlhdDatasetRatio(0)] = 1.0
		paraLlhd = paraLlhdReal + 1j*paraLlhdImag
		return paraLlhd

	def paraLlhd2Fitter(self, paraLlhd):
		paraFitter = np.empty(self.nmbParameters)
		paraFitter[self.indicesRealFitterParaNonNone] = np.real(paraLlhd[self.indicesLlhdParaNonNoneReal])
		paraFitter[self.indicesImagFitterParaNonNone] = np.imag(paraLlhd[self.indicesLlhdParaNonNoneImag])
		return paraFitter


	def paraLlhd2negLlhd(self, paraLlhd):
		'''
		Split paraLlhd (one long parameter list) into [ [<list of T-vectors in sectors>], <list of likelihood ratios>, <additional parameters>]) for the `negLlhd` function.
		The complex-valued additional parameters in paraLlhd are mapped to real ones.
		'''
		amplitudes = tuple(paraLlhd[self.paraLlhdStartSectors[i]:self.paraLlhdStartSectors[i+1]] for i in range(self.nmbSectors))
		datasetRatios = np.real(paraLlhd[self.paraLlhdStartDatasetRatios:self.paraLlhdStartStartAdditionParameters])
		return [ amplitudes ] + [datasetRatios] + list(np.real(paraLlhd[self.paraLlhdStartStartAdditionParameters:]))

	def paraNegLlhd2Llhd(self, paraNegLlhd):
		return np.hstack([np.hstack(paraNegLlhd[0])] + [paraNegLlhd[1]] + paraNegLlhd[2:])


	def gradLlhd2Fitter(self, gradLlhd, gradFitter):
		'''
		@param gradLlhd: complex-valued array of derivatives in the likelihood parameter space from autgrad.grad
		@param gradLlhd: real-valued array of derivatives in the fitter parameter space
		'''
		# complex conj: autograd requires this for complex valued parameters (complex valued parameters, but real optimization ;) )
		gradFitter[self.indicesRealFitterParaNonNone] =  np.real(gradLlhd[self.indicesLlhdParaNonNoneReal])
		gradFitter[self.indicesImagFitterParaNonNone] = -np.imag(gradLlhd[self.indicesLlhdParaNonNoneImag])


	def hessianLlhd2Fitter(self, hessianLlhd, hessianFitter):
		'''
		@param hessianLlhd: nmbLlhdParameters*2 x nmbLlhdParameters*2 real-valued matrix
		@param hessianFitter: nmbParameters x nmbParameters real-valued matrix
		'''
		indices = []
		for iFitterPara in xrange(self.nmbParameters):
			if iFitterPara in self.indicesRealFitterParaNonNone:
				indices.append(self.indicesRealFitterPara.index(iFitterPara)*2)
			elif iFitterPara in self.indicesImagFitterParaNonNone:
				indices.append(self.indicesImagFitterPara.index(iFitterPara)*2+1)
			else:
				raise Exception("Cannot find fitter parameter index in index mapping :(")
		hessianFitter[:,:] = hessianLlhd[:,indices][indices]


	def paraFitter2AmpsForRpwaFitresult(self, paraFitter):
		return self.paraFitter2Llhd(paraFitter)[:self.paraLlhdStartSectors[-1]]

	def paraFitter2DatasetRatiosForRpwaFitresult(self, paraFitter):
		datasetRatioParameters = self.paraLlhd2negLlhd(self.paraFitter2Llhd(paraFitter))[1]
		datasetRatiosNorm = np.sum(np.exp(datasetRatioParameters))
		datasetRatios = np.exp(datasetRatioParameters)/datasetRatiosNorm
		return datasetRatios

	def paraFitterCovMatrixIndicesForRpwaFitresult(self):
		indices = []
		for i in xrange(self.paraLlhdStartSectors[-1]):
			indices.append((self.indicesRealFitterPara[i] if self.indicesRealFitterPara[i] is not None else -1,
			                self.indicesImagFitterPara[i] if self.indicesImagFitterPara[i] is not None else -1))
		return indices


class ParameterMappingConnected(ParameterMapping):


	def __init__(self, model):
		ParameterMapping.__init__(self, model)

		self.parameterMappings = [likelihood.parameterMapping for likelihood in self.model.likelihood.likelihoods]
		self.nmbModels = len(self.parameterMappings)
		self.nmbParameters = np.sum([pm.nmbParameters for pm in self.parameterMappings])
		self.nmbLlhdParameters = np.sum([pm.nmbLlhdParameters for pm in self.parameterMappings])

		self.offsetsLlhdForBins = [0]
		self.offsetsFitterForBins = [0]
		for parameterMapping in self.parameterMappings:
			self.offsetsLlhdForBins.append(self.offsetsLlhdForBins[-1] + parameterMapping.nmbLlhdParameters)
			self.offsetsFitterForBins.append(self.offsetsFitterForBins[-1] + parameterMapping.nmbParameters)
		self.offsetsLlhdForBins = np.array(self.offsetsLlhdForBins, dtype=np.int64)
		self.offsetsFitterForBins = np.array(self.offsetsFitterForBins, dtype=np.int64)

	def paraFitterOfBin(self, paraFitter, iBin):
		'''
		@param paraFitter: Array of all fitter parameters
		@return fitter parameters of it i'th bin
		'''
		return paraFitter[self.offsetsFitterForBins[iBin]:self.offsetsFitterForBins[iBin+1]]


	def paraLlhdOfBin(self, paraLlhd, iBin):
		'''
		@param paraLlhd: Array of all Llhd parameters
		@return Llhd parameters of it i'th bin
		'''
		return paraLlhd[self.offsetsLlhdForBins[iBin]:self.offsetsLlhdForBins[iBin+1]]

	def paraFitter2Llhd(self, paraFitter):

		parasLlhd = np.empty(self.nmbLlhdParameters, dtype=np.complex128)
		for i,parameterMapping in enumerate(self.parameterMappings):
			parasLlhd[self.offsetsLlhdForBins[i]:self.offsetsLlhdForBins[i+1]] = \
			  parameterMapping.paraFitter2Llhd(paraFitter[self.offsetsFitterForBins[i]:self.offsetsFitterForBins[i+1]])

		return parasLlhd


	def paraLlhd2Fitter(self, paraLlhd):

		parasFitter = np.empty(self.nmbParameters)

		for i,parameterMapping in enumerate(self.parameterMappings):
			parasFitter[self.offsetsFitterForBins[i]:self.offsetsFitterForBins[i+1]]\
			  = parameterMapping.paraLlhd2Fitter(paraLlhd[self.offsetsLlhdForBins[i]:self.offsetsLlhdForBins[i+1]])
		return parasFitter


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
		for i, parameterMapping in enumerate(self.parameterMappings):
			parameterMapping.gradLlhd2Fitter(gradLlhd[self.offsetsLlhdForBins[i]:self.offsetsLlhdForBins[i+1]], gradFitter[self.offsetsFitterForBins[i]:self.offsetsFitterForBins[i+1]])


	def hessianLlhd2Fitter(self, hessianLlhd, hessianFitter):
		indices = []
		for iParaMapping, parameterMapping in enumerate(self.parameterMappings):
			for iFitterPara in xrange(parameterMapping.nmbParameters):
				if iFitterPara in parameterMapping.indicesRealFitterParaNonNone:
					indices.append(self.offsetsFitterForBins[iParaMapping] + parameterMapping.indicesRealFitterPara.index(iFitterPara) * 2)
				elif iFitterPara in parameterMapping.indicesImagFitterParaNonNone:
					indices.append(self.offsetsFitterForBins[iParaMapping] + parameterMapping.indicesImagFitterPara.index(iFitterPara) * 2 + 1)
				else:
					raise Exception("Cannot find fitter parameter index in index mapping :(")
		hessianFitter[:, :] = hessianLlhd[:, indices][indices]
