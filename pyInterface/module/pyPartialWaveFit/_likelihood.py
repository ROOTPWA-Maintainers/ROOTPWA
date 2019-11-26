'''
@author: F. Kaspar, S. Wallner
'''
# pylint: disable=E1120,E1101,W0221,W0122,W0212

from __future__ import absolute_import
from __future__ import division

import sys
import os
import inspect
import collections
import glob

import autograd
import autograd.numpy as np

from .. import utils
from ._parameterMapping import ParameterMapping

try: # python 2
# pylint: disable=W0622
	range = xrange
# pylint: enable=W0622
except NameError: # python 3
	pass


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

		self.datasetIndecesToProcess = [i for i in range(self.nmbDatasets) if self.nmbEvents[i] > 0]

		self.parameterMapping = parameterMapping

		self._valueAndGrad = autograd.value_and_grad(self._negLlhd)
		self._grad = autograd.grad(self._negLlhd)
		self._hessian = autograd.jacobian(self._gradientForHessian)
		self.countFdF = 0 # number of calls to f with gradient
		self.countF = 0 # number of calls to f without gradient

	def setParameters(self):
		pass # no parameters for this likelihood class


	def _getDumpData(self):
		'''
		@return: Dictionary with data dumped to file for `dumpToFile`
		'''
		dumpData = {}
		for varName in ['decayAmplitudes', 'accMatrices', 'normMatrices', 'normIntegrals']:
			for i, data in enumerate(getattr(self, varName)):
				for j, amps in enumerate(data):
					key = "{varName}_{i}_{j}".format(varName=varName, i=i, j=j)
					dumpData[key] = amps
		dumpData['parameterMapping'] = self.parameterMapping.toYaml()
		return dumpData


	def dumpToFile(self, filepath):
		'''
		Dump the current likelihood object to the given file
		@param filepath: Path to the output file, which must not exist
		'''
		if os.path.exists(filepath):
			utils.printErr("'{0}' exists already!".format(filepath))
			raise IOError()
		dumpData = self._getDumpData()
		dumpData['__class__'] = self.__class__.__name__
		np.savez(filepath, **dumpData)


	@classmethod
	def _fromDumpData(cls, dumpData):
		'''
		@return: Likelihood object from dump data
		'''
		kwargs = collections.defaultdict(list)
		for varName in ['decayAmplitudes', 'accMatrices', 'normMatrices', 'normIntegrals']:
			for key in sorted([f for f in dumpData.files if varName in f]):
				i, j = map(int, key.split('_')[1:])
				if i >= len(kwargs[varName]):
					kwargs[varName].append([])
				kwargs[varName][i].append(dumpData[key])
				assert len(kwargs[varName]) == i+1
				assert len(kwargs[varName][i]) == j+1
		parameterMappingStr = "{0}".format(dumpData['parameterMapping'].astype(str))
		kwargs['parameterMapping'] = ParameterMapping.fromYaml(parameterMappingStr)
		likelihood = cls(**kwargs)
		return likelihood


	@classmethod
	def fromFile(cls, filepath):
		'''
		@param filepath: Path to the input file
		@return: Likelihood object from the information dumped to the file
		@rtype: Likelihood
		'''
		if not os.path.exists(filepath):
			utils.printErr("'{0}' does not exist!".format(filepath))
			raise IOError()
		dumpData = np.load(filepath)
		clsLikelihood = None

		clsLikelihoodName = "{0}".format(dumpData['__class__'].astype(str))
		varOut = {}
		exec("clsLikelihood = {0}".format(clsLikelihoodName), globals(), varOut)
		clsLikelihood = varOut['clsLikelihood']
		return clsLikelihood._fromDumpData(dumpData)


	def negLlhd(self, transitionAmps, datasetRatioParameters):
		'''
		@return Negative log-likelihood for the given set of transition amplitudes
		'''
		datasetRatios, datasetRatiosLog = self.parameterMapping.datasetRatioParameters2DatasetRatios(datasetRatioParameters)
		negLlhd = 0.0
		for iDataset in self.datasetIndecesToProcess:
			nBar = 0.0
			dataTerm = 0.0
			for iSector in range(self.nmbSectors):
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

	def valueAndGradient(self, paraLlhd):
		negLL, gradLlhd = self._valueAndGrad(paraLlhd)
		gradLlhd = np.conj(gradLlhd)
		return negLL, gradLlhd


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

	def useNEvents(self, nmbEventsTotal, generator):
		'''
		Use only nEvents of the given set of data events.
		The nEvents are distributed over the data sets according to the initial data-set ratio.
		The nEvents are drawn for each data set from the data set using the given random-number generator.
		The order of the events is shuffled.
		@param nmbEventsTotal: Total number of events (all data sets) that should be used
		@param generator: Random number generator that is used
		'''
		nmbEventsOrig = self.nmbEvents
		if np.sum(nmbEventsOrig) < nmbEventsTotal:
			utils.printErr("Cannot set {0} events for likelihood when there are only {1} events".format(nmbEventsTotal, np.sum(nmbEventsOrig)))
			raise Exception()
		ratios = nmbEventsOrig/np.sum(nmbEventsOrig)

		nmbEvents = np.around(ratios * nmbEventsTotal).astype(int)
		nmbEvents[-1] = nmbEventsTotal - np.sum(nmbEvents[:-1]) # to have exactly nmbEventsTotal events
		nmbEvents = np.maximum(nmbEvents, 0)

		self.nmbEvents = nmbEvents
		self.datasetIndecesToProcess = [i for i in range(self.nmbDatasets) if self.nmbEvents[i] > 0]
		for iDataset, datasetAmplitudes in enumerate(self.decayAmplitudes):
			eventIndices = generator.randint(0, nmbEventsOrig[iDataset], nmbEvents[iDataset])
			for iSector, sectorAmplitudesOrig in enumerate(datasetAmplitudes):
				if sectorAmplitudesOrig.shape[1] == 1 and eventIndices.size > 0: # the decay amplitudes of all events have the same vale (e.g. flat wave)
					continue
				sectorAmplitudes = np.empty((sectorAmplitudesOrig.shape[0], eventIndices.size), dtype=sectorAmplitudesOrig.dtype)
				sectorAmplitudes[:,:] = sectorAmplitudesOrig[:,eventIndices]
				datasetAmplitudes[iSector] = sectorAmplitudes


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
		self.sectors = list(range(self.nmbSectors))[:-1] # apply regularization to all but the last one, which corresponds to the flat wave usually
		self.width = [0.5]*self.nmbSectors


	def _getDumpData(self):
		dumpData = Likelihood._getDumpData(self)
		for iSector, widths in enumerate(self.width):
			key = "width_{0}".format(iSector)
			dumpData[key] = np.array(widths)
		dumpData['sectors'] = np.array(self.sectors, dtype=np.int64)
		return dumpData


	@classmethod
	def _fromDumpData(cls, dumpData):
		likelihood = super(LikelihoodCauchy, cls)._fromDumpData(dumpData)
		for iSector in range(likelihood.nmbSectors):
			key = "width_{0}".format(iSector)
			likelihood.width[iSector] = dumpData[key]
		likelihood.sectors = dumpData['sectors']
		return likelihood


	def setParameters(self, width, correctAcceptance=False, sectors= None):
		if isinstance(width, float):
			if not correctAcceptance:
				self.width = [width] * len(self.nmbSectors)
			else:
				self.width = [width / np.sqrt(np.abs(accMatrix.diagonal())) for accMatrix in self.accMatrices[0]]
		else:
			if len(self.width) != self.nmbSectors:
				utils.printErr("Wrong number of sectors in Cauchy width!")
				raise Exception()
			if correctAcceptance:
				utils.printErr("'widthTimesAcceptance' can only be used with scalar width!")
				raise Exception()
			self.width = width

		if sectors is not None:
			if True in [ i < 0 or i >= self.nmbSectors for i in sectors]:
				utils.printErr("One of the sectors for the Cauchy regularization is out of range")
				raise Exception()
			self.sectors = sectors


	def negLlhd(self, transitionAmps, datasetRatios):
		negLlhd = Likelihood.negLlhd(self, transitionAmps, datasetRatios)
		for i in self.sectors:
			absT = np.abs(transitionAmps[i])
			negLlhd = negLlhd - np.sum(np.log(1.0/(1.0+absT**2/self.width[i]**2)))
		return negLlhd


class LikelihoodConnected(object):

	def __init__(self, likelihoods, binWidths, parameterMapping):

		self.likelihoods = likelihoods
		self.binWidths = binWidths
		self.binWidthsNormalization = np.empty((len(likelihoods), likelihoods[0].parameterMapping.nmbLlhdParameters-likelihoods[0].parameterMapping.nmbDatasets))
		for i in range(self.binWidthsNormalization.shape[0]):
			self.binWidthsNormalization[i,:] = np.sqrt(1.0/binWidths[i]*0.02)

		self.countFdF = 0  # number of calls to f with gradient
		self.countF = 0  # number of calls to f without gradient

		# this is a bad hack :(
		self.nmbEvents = np.sum([likelihood.nmbEvents for likelihood in self.likelihoods])
		self.nmbParameters = np.sum([likelihood.parameterMapping.nmbParameters for likelihood in self.likelihoods])

		self.parameterMapping = parameterMapping

		self._connectionValueAndGrad = autograd.value_and_grad(self._connection)
		self._connectionGrad = autograd.grad(self._connection)
		self._connectionHessian = autograd.jacobian(self._gradientForConnectionHessian)


	def setParameters(self, **kwargs):
		'''
		All keyword arguments are passed to the individual likelihoods
		'''
		if kwargs:
			for likelihood in self.likelihoods:
				likelihood.setParameters(**kwargs)


	def _getDumpData(self):
		'''
		@return: Dictionary with data dumped to file for `dumpToFolder`
		'''
		dumpData = {}
		dumpData['binWidths'] = np.array(self.binWidths)
		dumpData['parameterMapping'] = self.parameterMapping.toYaml()
		return dumpData


	def dumpToFolder(self, folderpath):
		'''
		Dump the current likelihood object to the given folder
		@param folderpath: Path to the output folder, which must not exist
		'''
		if os.path.exists(folderpath):
			utils.printErr("'{0}' exists already!".format(folderpath))
			raise IOError()
		os.makedirs(folderpath)
		dumpData = self._getDumpData()
		dumpData['__class__'] = self.__class__.__name__
		np.savez(os.path.join(folderpath, 'connected_data.npz'), **dumpData)
		for iBin, likelihood in enumerate(self.likelihoods):
			filepath = os.path.join(folderpath, "bin{0:04d}_data.npz".format(iBin))
			likelihood.dumpToFile(filepath)


	@classmethod
	def _fromDumpData(cls, likelihoods, dumpData):
		'''
		@return: Likelihood object from dump data
		'''
		kwargs = {}
		kwargs['likelihoods'] = likelihoods
		kwargs['binWidths'] = dumpData['binWidths']
		parameterMappingStr = "{0}".format(dumpData['parameterMapping'].astype(str))
		kwargs['parameterMapping'] = ParameterMapping.fromYaml(parameterMappingStr)
		likelihood = cls(**kwargs)
		return likelihood


	@classmethod
	def fromFolder(cls, folderpath):
		'''
		@param folderpath: Path to the input folder
		@return: Likelihood object from the information dumped to the folder
		@rtype: Likelihood
		'''
		if not os.path.exists(folderpath):
			utils.printErr("'{0}' does not exist!".format(folderpath))
			raise IOError()
		likelihoods = []
		for iBin, likelihoodDumpFile in enumerate(sorted(glob.glob(os.path.join(folderpath, "bin????_data.npz")))):
			iBinFromFile = int(os.path.basename(likelihoodDumpFile).split('_')[0].lstrip('bin'))
			if iBin != iBinFromFile:
				utils.printErr("Bin index and filename for likelihood dump does not agree!")
				raise Exception()
			likelihoods.append(Likelihood.fromFile(likelihoodDumpFile))
		clsLikelihood = None
		dumpData = np.load(os.path.join(folderpath, 'connected_data.npz'))
		clsLikelihoodName = "{0}".format(dumpData['__class__'].astype(str))
		varOut = {}
		exec("clsLikelihood = {0}".format(clsLikelihoodName), globals(), varOut)
		clsLikelihood = varOut['clsLikelihood']
		return clsLikelihood._fromDumpData(likelihoods, dumpData)


# pylint: disable=W0613,R0201
	def _connection(self, paraLlhd):
		return 0.
# pylint: enable=W0613,R0201


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


	def valueAndGradient(self, paraLlhd):
		negLL = 0
		grads = []
		for iLikelihood, likelihood in enumerate(self.likelihoods):
			paraLlhdBin = self.parameterMapping.paraLlhdOfBin(paraLlhd, iLikelihood)
			negLLBin, gradBin = likelihood.valueAndGradient(paraLlhdBin)
			negLL += negLLBin
			grads.append(gradBin)
		return negLL, np.hstack(grads)


	def value(self, paraLlhd):
		negLL = 0
		for iLikelihood, likelihood in enumerate(self.likelihoods):
			paraLlhdBin = self.parameterMapping.paraLlhdOfBin(paraLlhd, iLikelihood)
			negLLBin = likelihood._negLlhd(paraLlhdBin)
			negLL += negLLBin
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


class LikelihoodConnectedGauss(LikelihoodConnected):

	def __init__(self, likelihoods, binWidths, parameterMapping, strength = 0.08, scaleStrengthByEvents=False):
		LikelihoodConnected.__init__(self, likelihoods, binWidths, parameterMapping)

		# set connection strength on init
		self.strength = None
		self.setParameters(strength, scaleStrengthByEvents)


	def _connection(self, paraLlhd):
		paras = np.reshape(paraLlhd, (len(self.likelihoods), -1))[:,:self.likelihoods[0].parameterMapping.paraLlhdStartSectors[-1]]
		# normalize to 20 MeV bins
		parasNormed = paras*self.binWidthsNormalization
		return np.sum(self.strength * (np.abs(parasNormed[1:] - parasNormed[:-1]) ** 2))


	def setParameters(self, strength=0.08, scaleStrengthByEvents=False, **kwargs):
		'''
		@param strength: Strength parameter of the connection term
		@param scaleStrengthByEvents: Scale the strength parameter by the number of events in each bin
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

		super(LikelihoodConnectedGauss, self).setParameters(**kwargs)

	def _getDumpData(self):
		dumpData = LikelihoodConnected._getDumpData(self)
		dumpData['strength'] = self.strength
		return dumpData

	@classmethod
	def _fromDumpData(cls, likelihoods, dumpData):
		likelihood = super(LikelihoodConnectedGauss, cls)._fromDumpData(likelihoods, dumpData)
		likelihood.strength = dumpData['strength']
		return likelihood


class LikelihoodConnectedCauchy(LikelihoodConnected):

	def __init__(self, likelihoods, binWidths, parameterMapping, strength = 0.08):
		LikelihoodConnected.__init__(self, likelihoods, binWidths, parameterMapping)

		# set connection strength on init
		self.strength = None
		self.setParameters(strength)

	def _connection(self, paraLlhd):
		paras = np.reshape(paraLlhd, (len(self.likelihoods), -1))[:,:self.likelihoods[0].parameterMapping.paraLlhdStartSectors[-1]]
		# normalize to 20 MeV bins
		parasNormed = paras*self.binWidthsNormalization
		return np.sum(np.log(1.0+self.strength * (np.abs(parasNormed[1:] - parasNormed[:-1]) ** 2)))

	def setParameters(self, strength=0.08, **kwargs):
		'''
		@param strength: Strength parameter of the connection term
		@param scaleStrengthByEvents: Scale the strength parameter by the number of events in each bin
		All other keyword arguments are passed to the individual likelihoods
		'''
		self.strength = strength
		super(LikelihoodConnectedCauchy, self).setParameters(**kwargs)

	def _getDumpData(self):
		dumpData = LikelihoodConnected._getDumpData(self)
		dumpData['strength'] = self.strength
		return dumpData

	@classmethod
	def _fromDumpData(cls, likelihoods, dumpData):
		likelihood = super(LikelihoodConnectedCauchy, cls)._fromDumpData(likelihoods, dumpData)
		likelihood.strength = dumpData['strength']
		return likelihood


class LikelihoodConnectedFFT(LikelihoodConnected):

	def __init__(self, likelihoods, binWidths, parameterMapping, strength=0.08, ran=5):

		if np.sum(np.abs(np.array(binWidths) - binWidths[0]) > 1e-7) != 0:
			raise Exception("binWidths must all be equal for FFT likelihood!")

		LikelihoodConnected.__init__(self, likelihoods, binWidths, parameterMapping)

		# factor of 2 due to artifical periodictiy (compare connection term)
		self.freq = np.fft.fftfreq(2*len(self.likelihoods), d=self.binWidths[0])
		# set connection strength on init
		'''
		@todo: please explain this parameters
		'''
		self.strength = None
		self.ran = None
		self.setParameters(strength, ran)

	def setParameters(self, strength=0.08, ran=5, **kwargs):
		'''
		@param strength: Strength parameter of the connection term
		All other keyword arguments are passed to the individual likelihoods
		'''
		self.strength = strength
		self.ran = ran

		super(LikelihoodConnectedFFT, self).setParameters(**kwargs)


	# suppress frequencies above a certain frequency threshold -> enforce smoothness of the amplitudes
	# this is trying to connect the amplitudes between bins via a prior/penalty on frequency
	# the idea is to provide a "backward" version the smooth fits that IFT can produce
	def _connection(self, paraLlhd):
		paras = np.reshape(paraLlhd, (len(self.likelihoods), -1))[:,:self.likelihoods[0].parameterMapping.paraLlhdStartSectors[-1]]

		negConnection = 0
		for i in range(paras.shape[1]):
			# idea by S. Wallner: add mirrored values to prevent problems with DFT periodictiy
			traf = np.fft.fft(np.concatenate( [paras[:,i],paras[:,i][::-1]],axis=0)  )
			negConnection = negConnection + np.sum( self.strength*np.abs(traf[np.abs(self.freq) > self.ran])**2 )

		return negConnection

	def _getDumpData(self):
		raise NotImplementedError("Needs to be implemented")


	@classmethod
	def _fromDumpData(cls, likelihoods, dumpData):
		raise NotImplementedError("Needs to be implemented")
