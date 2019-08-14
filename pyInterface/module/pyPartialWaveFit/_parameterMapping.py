'''
@author: F. Kaspar, S. Wallner
'''
# pylint: disable=E1101,W0212

from __future__ import absolute_import
from __future__ import division

import copy
import inspect
import autograd.numpy as np
import yaml

try: # python 2
# pylint: disable=W0622
	range = xrange
# pylint: enable=W0622
except NameError: # python 3
	pass



class ParameterMapping(object):
	'''
	The parameters transforms between three parameter spaces
		- The fitter parameters: array of real-valued parameters
		- The likelihood parameters: array of complex-valued parameters
			- The space in which autograd works
		- The parameters of the likelihood function (paraNegLlhd): List of parameter arrays, specific to the likelihood function
	'''

	def __init__(self):

		self.nmbParameters = None
		self.indicesReferenceWaveFitterPara = []
		self.indicesReferenceWaveInSectors = []

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

	def datasetRatioParameters2DatasetRatios(self, datasetRatioParameters):
		raise NotImplementedError("Needs to be implemented in specialized class!")

	def datasetDatasetRatios2RatioParameters(self, datasetRatios):
		raise NotImplementedError("Needs to be implemented in specialized class!")


	def _prepareToYaml(self):
		for key, attribute in self.__dict__.items():
			if isinstance(attribute, np.ndarray):
				self.__dict__[key] = map(int, list(attribute))

	def _recoverAfterYaml(self):
		for key, attribute in self.__dict__.items():
			if isinstance(attribute, list) and all([isinstance(i, int) for i in attribute]):
				self.__dict__[key] = np.array(attribute, dtype=np.int64)

	def toYaml(self):
		'''
		@return: String representing the parameter mapping in YAML format
		'''
		obj = copy.deepcopy(self)
		obj._prepareToYaml()
		return yaml.dump(obj)

	@classmethod
	def fromYaml(cls, yamlString):
		'''
		@return: ParameterMapping from the given YAML string
		@rtype: ParameterMapping
		'''
		parameterMapping = yaml.load(yamlString, Loader = yaml.FullLoader)
		parameterMapping._recoverAfterYaml()
		return parameterMapping


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

	def __init__(self, model, realWaves = None, zeroWaves = None, zeroDatasetRatios = None):
		'''
		@param realWaves: List lists of waves whose amplitude should be fixed to real, i.e. the imaginary part is fixed to zero. One list for each incoherent sector.
		@param zeroWaves: List of waves whose amplitude should be fixed to zero. Fixing of the reference wave is ignored.
		@param zeroDatasetRatios: List of data set indices whose ratio should be 0
		'''
		ParameterMapping.__init__(self)

		self.wavesInSectors = model.wavesInSectors
		self.nmbWavesInSectors = np.array([ len(wavesInSector) for wavesInSector in self.wavesInSectors ], dtype=np.int64)
		self.nmbSectors = len(self.nmbWavesInSectors)
		self.nmbDatasets = model.nmbDatasets

		if realWaves is None:
			realWaves = [ [] ] * self.nmbSectors

		if len(realWaves) != self.nmbSectors:
			raise ValueError("List of real-valued wave does not have n sectors.")

		if zeroWaves is None:
			zeroWaves = []

		if zeroDatasetRatios is None:
			zeroDatasetRatios = []
		if 0 in zeroDatasetRatios:
			raise ValueError("Cannot set data set ratio of 0th data set to 0!")
		if zeroDatasetRatios and max(zeroDatasetRatios) >= model.nmbDatasets:
			raise ValueError("Cannot set data set ratio of {0}nd data set to 0! Index out of range!".format(max(zeroDatasetRatios)))
		self.zeroDatasetRatios = np.array(zeroDatasetRatios, dtype=np.int64)

		negLlhdFunction = model.clsLikelihood.negLlhd
		argsNegLlhdFunction = inspect.getargspec(negLlhdFunction).args
		self.nmbAdditionalParameters = len(argsNegLlhdFunction) -3 # the first one is self, the second one are the transition amplitudes, the third one are the data-set ratios
		self.nmbAmplitudeParameters = int(np.sum([ 2*n-len(realWaves[i]) for i,n in enumerate(self.nmbWavesInSectors)]) - 2*len(zeroWaves))
		self.nmbDatasetRatioParameters = self.nmbDatasets-1-len(self.zeroDatasetRatios)
		self.nmbParameters = self.nmbAmplitudeParameters + self.nmbDatasetRatioParameters + self.nmbAdditionalParameters


		# build indices of paraLlhd where a new sector starts
		self.paraLlhdStartSectors = [0]
		for i, nWavesInSector in enumerate(map(int, self.nmbWavesInSectors)):
			self.paraLlhdStartSectors.append(nWavesInSector + self.paraLlhdStartSectors[i])
		self.paraLlhdStartDatasetRatios = self.paraLlhdStartSectors[-1]
		self.paraLlhdStartStartAdditionParameters = self.paraLlhdStartDatasetRatios + self.nmbDatasets
		self.nmbLlhdParameters = self.paraLlhdStartSectors[-1] + self.nmbDatasets + self.nmbAdditionalParameters

		self.paraLlhdStartSectors = np.array(self.paraLlhdStartSectors, dtype=np.int64)


		# build index mapping between paraLlhd and paraFitter

		self._buildIndexMapping(model, zeroWaves, realWaves, zeroDatasetRatios)
		self._buildParameterNames(zeroWaves, realWaves, argsNegLlhdFunction)


	def _buildIndexMapping(self, model, zeroWaves, realWaves, zeroDatasetRatios):
		iFitterPara = 0
		self.indicesRealFitterPara = []  # mapping paraLlhd index -> real value fitter parameter index
		self.indicesImagFitterPara = []  # mapping paraLlhd index -> imag value fitter parameter index
		self.indicesReferenceWaveFitterPara = []  # list of fitter parameter indices of the reference waves
		self.indicesReferenceWaveInSectors = []   # list of paraNegLlhd indices in sectors of the reference waves
		for iSector, nWavesInSector in enumerate(self.nmbWavesInSectors):
			for iWave in range(nWavesInSector):
				if self.wavesInSectors[iSector][iWave] not in zeroWaves:
					self.indicesRealFitterPara.append(iFitterPara) # add real part to fitter parameters
					if self.wavesInSectors[iSector][iWave] == model.referenceWaves[iSector]:
						self.indicesReferenceWaveFitterPara.append(iFitterPara)
						self.indicesReferenceWaveInSectors.append(iWave)
					iFitterPara += 1

					if self.wavesInSectors[iSector][iWave] not in realWaves[iSector]:
						self.indicesImagFitterPara.append(iFitterPara) # add imaginary part to fitter parameters
						iFitterPara += 1
					else:
						self.indicesImagFitterPara.append(None)
				else:
					self.indicesRealFitterPara.append(None)
					self.indicesImagFitterPara.append(None)

		for iDataset in range(self.nmbDatasets):
			if iDataset == 0 or iDataset in zeroDatasetRatios:
				self.indicesRealFitterPara.append(None)
			else:
				self.indicesRealFitterPara.append(iFitterPara)
				iFitterPara += 1
			self.indicesImagFitterPara.append(None)

		for _ in range(self.nmbAdditionalParameters):
			self.indicesRealFitterPara.append(iFitterPara)
			iFitterPara += 1
			self.indicesImagFitterPara.append(None)

		assert iFitterPara == self.nmbParameters

		# build list of indices which are not None, i.e. which are real parameters
		self.indicesRealFitterParaNonNone = np.array([i for i in self.indicesRealFitterPara if i is not None], dtype=np.int64)
		self.indicesLlhdParaNonNoneReal = np.array([j for j, i in enumerate(self.indicesRealFitterPara) if i is not None], dtype=np.int64)
		self.indicesImagFitterParaNonNone = np.array([i for i in self.indicesImagFitterPara if i is not None], dtype=np.int64)
		self.indicesLlhdParaNonNoneImag = np.array([j for j, i in enumerate(self.indicesImagFitterPara) if i is not None], dtype=np.int64)

		self.indicesReferenceWaveFitterPara = np.array(self.indicesReferenceWaveFitterPara, dtype=np.int64)


	def _buildParameterNames(self, zeroWaves, realWaves, argsNegLlhdFunction):
		# build parameter names
		self.paraNamesFitter = []
		for iSector, wavesInSector in enumerate(self.wavesInSectors):
			for wave in wavesInSector:
				if wave not in zeroWaves:
					self.paraNamesFitter.append("V{s}_{w}_re".format(s=iSector, w=wave))
					if wave not in realWaves[iSector]:
						self.paraNamesFitter.append("V{s}_{w}_im".format(s=iSector, w=wave))

		for iDatasetRatioParameter in range(self.nmbDatasetRatioParameters):
			self.paraNamesFitter.append("r_{0}".format(iDatasetRatioParameter+1))
		for iAdditionalParameter in range(self.nmbAdditionalParameters):
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
		for iDataset in self.zeroDatasetRatios:
			paraLlhdReal[self.paraLlhdStartDatasetRatios+iDataset] = -np.inf # only for exponential mapping
		paraLlhd = paraLlhdReal + 1j*paraLlhdImag
		return paraLlhd

	def paraLlhd2Fitter(self, paraLlhd):
		paraFitter = np.empty(self.nmbParameters)
		paraFitter[self.indicesRealFitterParaNonNone] = np.real(paraLlhd[self.indicesLlhdParaNonNoneReal])
		paraFitter[self.indicesImagFitterParaNonNone] = np.imag(paraLlhd[self.indicesLlhdParaNonNoneImag])
		return paraFitter


	def paraLlhd2negLlhd(self, paraLlhd):
		'''
		Split paraLlhd (one long parameter list) into [ [<list of T-arrays in sectors>], <arrayof likelihood ratios>, <additional parameters>, ...]) for the `negLlhd` function.
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
		for iFitterPara in range(self.nmbParameters):
			if iFitterPara in self.indicesRealFitterParaNonNone:
				indices.append(self.indicesRealFitterPara.index(iFitterPara)*2)
			elif iFitterPara in self.indicesImagFitterParaNonNone:
				indices.append(self.indicesImagFitterPara.index(iFitterPara)*2+1)
			else:
				raise Exception("Cannot find fitter parameter index in index mapping :(")
		hessianFitter[:,:] = hessianLlhd[:,indices][indices]


	def paraFitter2AmpsForRpwaFitresult(self, paraFitter):
		paraNegLlhd = self.paraLlhd2negLlhd(self.paraFitter2Llhd(paraFitter))
		# flip reference waves to positive
		for iSector in range(self.nmbSectors):
			iRef = self.indicesReferenceWaveInSectors[iSector]
			if paraNegLlhd[0][iSector][iRef] < 0.0:
				ampsInSector = paraNegLlhd[0][iSector]
				ampsInSector *= -1
		return self.paraNegLlhd2Llhd(paraNegLlhd)[:self.paraLlhdStartSectors[-1]]

	def paraFitter2DatasetRatiosForRpwaFitresult(self, paraFitter):
		datasetRatioParameters = self.paraLlhd2negLlhd(self.paraFitter2Llhd(paraFitter))[1]
		datasetRatios, _ = self.datasetRatioParameters2DatasetRatios(datasetRatioParameters)
		return datasetRatios


	def paraFitterCovMatrixIndicesOfDatasetRatiosForRpwaFitresult(self):
		indices = []
		for i in range(self.paraLlhdStartDatasetRatios, self.paraLlhdStartStartAdditionParameters):
			indices.append(self.indicesRealFitterPara[i] if self.indicesRealFitterPara[i] is not None else -1)
		return indices


	def paraFitterCovMatrixIndicesForRpwaFitresult(self):
		indices = []
		for i in range(self.paraLlhdStartSectors[-1]):
			indices.append((self.indicesRealFitterPara[i] if self.indicesRealFitterPara[i] is not None else -1,
			                self.indicesImagFitterPara[i] if self.indicesImagFitterPara[i] is not None else -1))
		return indices

	def paraFitterJacobianMatrixForRpwa(self, paraFitter):
		'''
		@return: Jacobian matrix J[i,j] from Rpwa parameters (i) to fitter parameter (j)
		'''
		jacobian = np.zeros((self.nmbParameters, self.nmbParameters))
		np.fill_diagonal(jacobian, 1.0)

		# build transformation for data set ratios
		indices = []
		mappingIndex2DatasetRatioIndex = {}
		for iLlhd in range(self.paraLlhdStartDatasetRatios, self.paraLlhdStartStartAdditionParameters):
			if self.indicesRealFitterPara[iLlhd] is not None:
				indices.append(self.indicesRealFitterPara[iLlhd])
				mappingIndex2DatasetRatioIndex[indices[-1]] = iLlhd-self.paraLlhdStartDatasetRatios
		datasetRatios = self.paraFitter2DatasetRatiosForRpwaFitresult(paraFitter)
		for i in indices:
			ratioI = datasetRatios[mappingIndex2DatasetRatioIndex[i]]
			for j in indices:
				ratioJ = datasetRatios[mappingIndex2DatasetRatioIndex[j]]
				jacobian[i,j] = -ratioJ*ratioI
			jacobian[i,i] += ratioI
		return jacobian


	def datasetRatioParameters2DatasetRatios(self, datasetRatioParameters):
		datasetRatiosNorm = np.sum(np.exp(datasetRatioParameters))
		datasetRatiosLog = datasetRatioParameters-np.log(datasetRatiosNorm)
		datasetRatios = np.exp(datasetRatiosLog)
		return datasetRatios, datasetRatiosLog

	def datasetDatasetRatios2RatioParameters(self, datasetRatios):
		ratioParameters = np.empty_like(datasetRatios)
		for i, value in enumerate(datasetRatios):
			if value != 0.0:
				ratioParameters[i] = np.log(value)-np.log(datasetRatios[0])
			else:
				ratioParameters[i] = - np.inf
		return ratioParameters


class ParameterMappingConnected(ParameterMapping):

	def __init__(self, parameterMappings):
		ParameterMapping.__init__(self)

		self.parameterMappings = parameterMappings
		self.nmbModels = len(self.parameterMappings)
		self.nmbParameters = int(np.sum([pm.nmbParameters for pm in self.parameterMappings]))
		self.nmbLlhdParameters = int(np.sum([pm.nmbLlhdParameters for pm in self.parameterMappings]))

		self.offsetsLlhdForBins = [0]
		self.offsetsFitterForBins = [0]
		self.indicesImagFitterPara = []
		self.indicesRealFitterPara = []
		self.indicesReferenceWaveFitterPara = []

		for parameterMapping in self.parameterMappings:
			self.indicesImagFitterPara += [idx + self.offsetsFitterForBins[-1] if idx is not None else None for idx in parameterMapping.indicesImagFitterPara]
			self.indicesRealFitterPara += [idx + self.offsetsFitterForBins[-1] if idx is not None else None for idx in parameterMapping.indicesRealFitterPara]

			self.indicesReferenceWaveFitterPara += [idx + self.offsetsFitterForBins[-1] for idx in parameterMapping.indicesReferenceWaveFitterPara]

			self.offsetsLlhdForBins.append(self.offsetsLlhdForBins[-1] + parameterMapping.nmbLlhdParameters)
			self.offsetsFitterForBins.append(self.offsetsFitterForBins[-1] + parameterMapping.nmbParameters)
		self.indicesReferenceWaveFitterPara = np.array(self.indicesReferenceWaveFitterPara, dtype=np.int64)

		self.offsetsLlhdForBins = np.array(self.offsetsLlhdForBins, dtype=np.int64)
		self.offsetsFitterForBins = np.array(self.offsetsFitterForBins, dtype=np.int64)

	def _prepareToYaml(self):
		ParameterMapping._prepareToYaml(self)
		for parameterMapping in self.parameterMappings:
			parameterMapping._prepareToYaml()

	def _recoverAfterYaml(self):
		ParameterMapping._recoverAfterYaml(self)
		for parameterMapping in self.parameterMappings:
			parameterMapping._recoverAfterYaml()

	def paraFitterOfBin(self, paraFitter, iBin):
		'''
		@param paraFitter: Array of all fitter parameters
		@return fitter parameters of it i'th bin
		'''
		return paraFitter[self.offsetsFitterForBins[iBin]:self.offsetsFitterForBins[iBin+1]]

	def individualParaFitter2paraFitter(self, parasFitter):
		'''
		@return: complete list of fitter parameter from bin-individual fitter parameters
		'''
		if len(parasFitter) != self.nmbModels:
			raise Exception("Wrong length of `parasFitter`!")

		paraFitter = np.empty(self.nmbParameters)
		for i in range(self.nmbModels):
			paraFitter[self.offsetsFitterForBins[i]:self.offsetsFitterForBins[i+1]] = parasFitter[i]
		return paraFitter


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
			for iFitterPara in range(parameterMapping.nmbParameters):
				if iFitterPara in parameterMapping.indicesRealFitterParaNonNone:
					indices.append(self.offsetsFitterForBins[iParaMapping] + parameterMapping.indicesRealFitterPara.index(iFitterPara) * 2)
				elif iFitterPara in parameterMapping.indicesImagFitterParaNonNone:
					indices.append(self.offsetsFitterForBins[iParaMapping] + parameterMapping.indicesImagFitterPara.index(iFitterPara) * 2 + 1)
				else:
					raise Exception("Cannot find fitter parameter index in index mapping :(")
		hessianFitter[:, :] = hessianLlhd[:, indices][indices]
