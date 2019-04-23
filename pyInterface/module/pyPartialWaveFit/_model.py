'''
@author: F. Kaspar, S. Wallner
'''

import numpy as np
import pyRootPwa
import pyRootPwa.utils

from _likelihood import Likelihood, LikelihoodConnected
from _parameterMapping import ParameterMappingRpwa, ParameterMappingConnected



class Model(object):
	'''
	classdocs
	'''


	def __init__(self, clsLikelihood, clsParameterMapping):
		'''
		Constructor
		@param clsLikelihood: Likelihood class
		@param clsParameterMapping: Parameter mapping class
		'''

		self.clsLikelihood = clsLikelihood
		self.clsParameterMapping = clsParameterMapping
		self.likelihood = None
		self.parameterMapping = None
		self.wavesInSectors = []
		self.nmbWaves = None
		self.waveNames = None
		self.nmbDatasets = None
		self.datasetLabels = []


	def amplitudeNames(self):
		raise NotImplementedError("Must be implemented in specialized class")


class ModelRpwa(Model):


	def __init__(self, clsLikelihood, clsParameterMapping):
		Model.__init__(self, clsLikelihood, clsParameterMapping)
		self.rankPosRefl = None
		self.rankNegRefl = None
		self.waveNames = None
		self.waveNamesPosRefl = None
		self.waveNamesNegRefl = None
		self.waveNamesFixed = None
		self.integralsPosRefl = None  # list of tuple of (normIntegralMatrix, accIntegralMatrix, normIntegrals) one for each data set
		self.integralsNegRefl = None  # list of tuple of (normIntegralMatrix, accIntegralMatrix, normIntegrals) one for each data set
		self.multibin = None


	def initModelInBin(self, fileManagerOrConfigfile, multiBin, waveListFileName, rankPosRefl, rankNegRefl, datasets = None):
		'''
		@param fileManagerOrConfigfile: Can be a fileManager object or a path to the config file
		@param multiBin: Can be a multibin or the integer bin-id
		@param waveListFileName: Path to the wavelist file
		@param rankPosRefl: Rank of positive reflectivity waves
		@param rankNegRefl: Rank of negative reflectivity waves
		@param datasets: List of data-sets to be fitted. If None, all data-sets are used.
		'''

		if isinstance(fileManagerOrConfigfile, str):
			config = pyRootPwa.rootPwaConfig()
			if not config.initialize(fileManagerOrConfigfile):
				pyRootPwa.utils.printErr("loading config file '" + fileManagerOrConfigfile + "' failed. Aborting...")
				raise Exception()
			pyRootPwa.core.particleDataTable.readFile(config.pdgFileName)
			fileManager = pyRootPwa.loadFileManager(config.fileManagerPath)
			if not fileManager:
				pyRootPwa.utils.printErr("loading the file manager failed. Aborting...")
				raise Exception()
		else:
			fileManager = fileManagerOrConfigfile

		if isinstance(multiBin, int):
			if multiBin < 0:
				pyRootPwa.utils.printErr("bin < 0 (" + str(multiBin) + "). Aborting...")
				raise Exception()
			elif multiBin >= len(fileManager.binList):
				pyRootPwa.utils.printErr("bin out of range (" + str(multiBin) + ">=" + str(len(fileManager.binList)) + "). Aborting...")
				raise Exception()
			multiBin = fileManager.binList[multiBin]
		if not multiBin in fileManager.binList:
			pyRootPwa.utils.printErr("Cannot find multibin (" + str(multiBin) + ") in binlist! Aborting...")
			raise Exception()

		if datasets is None:
			datasets = fileManager.datasetLabels
		self.nmbDatasets = len(datasets)
		self.datasetLabels = datasets

		self.initModelWithFiles(eventAndAmpFileDicts=[fileManager.getEventAndAmplitudeFilePathsInBin(multiBin, pyRootPwa.core.eventMetadata.REAL, d ) for d in datasets],
		                        normIntegralFileNames=[fileManager.getIntegralFilePath(multiBin, pyRootPwa.core.eventMetadata.GENERATED, d) for d in datasets],
		                        accIntegralFileNames=[fileManager.getIntegralFilePath(multiBin, pyRootPwa.core.eventMetadata.ACCEPTED, d) for d in datasets],
		                        multiBin=multiBin,
		                        waveListFileName=waveListFileName,
		                        waveDescriptions=fileManager.getWaveDescriptions(),
		                        rankPosRefl=rankPosRefl,
		                        rankNegRefl=rankNegRefl,
		                        datasets = datasets)


	def initModelWithFiles(self,
	                       eventAndAmpFileDicts,
	                       normIntegralFileNames, accIntegralFileNames,
	                       multiBin,
	                       waveListFileName, waveDescriptions,
	                       rankPosRefl,
	                       rankNegRefl,
	                       datasets):

		self.rankPosRefl = rankPosRefl
		self.rankNegRefl = rankNegRefl
		self.multibin = multiBin

		self.referenceWaves = self._initWaveNames(multiBin, waveListFileName, waveDescriptions)
		emptyDatasets = []
		normIntegralMatrices = []
		accIntegralMatrices = []
		normIntegrals = []
		decayAmps = []
		self.integralsPosRefl = []
		self.integralsNegRefl = []
		for iDataset, (eventAndAmpFileDict, normIntegralFileName, accIntegralFileName, dataset) in \
		   enumerate(zip(eventAndAmpFileDicts, normIntegralFileNames, accIntegralFileNames, datasets)):
			normIntegralMatrices.append(list())
			accIntegralMatrices.append(list())
			normIntegrals.append(list())
			decayAmps.append(list())
			normIntegralMatrixFull, accIntegralMatrixFull, normIntegralsFull, totAcc = loadMatrices(normIntegralFileName, accIntegralFileName, self.waveNames)
			for tag, waveNames, rank in zip(['pos', 'neg'], [self.waveNamesPosRefl, self.waveNamesNegRefl], [self.rankPosRefl, self.rankNegRefl]):
				if not waveNames:
					continue
				integrals = (buildSubMatrix(normIntegralMatrixFull, self.waveNames, waveNames),
							buildSubMatrix(accIntegralMatrixFull, self.waveNames, waveNames),
							buildSubList(normIntegralsFull, self.waveNames, waveNames))
				if tag == 'pos':
					self.integralsPosRefl.append(integrals)
				if tag == 'neg':
					self.integralsNegRefl.append(integrals)
				decayAmpsFull = loadAmplitudes(eventAndAmpFileDict, waveNames, multiBin, integrals[2])
				if decayAmpsFull.shape[1] == 0 and iDataset not in emptyDatasets:
					pyRootPwa.utils.printWarn("No events to fit in multibin: {0} for dataset {1}".format(multiBin, dataset))
					emptyDatasets.append(iDataset)

				for iRank in xrange(rank):
					decayAmps[-1].append(decayAmpsFull[iRank:])
					normIntegralMatrices[-1].append(integrals[0][iRank:, iRank:])
					accIntegralMatrices[-1].append(integrals[1][iRank:, iRank:])
					normIntegrals[-1].append(integrals[2][iRank:])
			# add flat wave
			decayAmps[-1].append(np.ones((1, 1), dtype=np.complex128))
			normIntegralMatrices[-1].append(np.ones((1, 1), dtype=np.complex128))
			accIntegralMatrices[-1].append(np.full((1, 1), totAcc, dtype=np.complex128))
			normIntegrals[-1].append(np.ones((1), dtype=np.float64))

		# needs to be added afterwards as the integral matrices do not know a flat wave
		self.wavesInSectors.append(["flat"])
		self.referenceWaves.append("flat")
		self.waveNames.append("flat")
		self.nmbWaves = len(self.waveNames)


		# buildParameterMapping
		self.parameterMapping = self.clsParameterMapping(self, [[w] for w in self.referenceWaves], self.waveNamesFixed, emptyDatasets)

		self.likelihood = self.clsLikelihood(decayAmps, accIntegralMatrices, normIntegralMatrices, normIntegrals, self.parameterMapping)

		return True


	def _initWaveNames(self, multiBin, waveListFileName, waveDescriptions):
		waveDescriptionActive, referenceWavesDefinitions = pyRootPwa.utils.getWaveDescriptionActiveFromWavelist(waveListFileName, waveDescriptions, multiBin)

		self.waveNames = [waveName for waveName, _, _ in waveDescriptionActive]

		# group waves by reflectivity
		self.waveNamesPosRefl = []
		self.waveNamesNegRefl = []
		for waveName, waveDescription, _ in waveDescriptionActive:
			status, topology = waveDescription.constructDecayTopology()
			if status:
				refl = topology.XParticle().reflectivity
				if refl > 0:
					self.waveNamesPosRefl.append(waveName)
				elif refl < 0:
					self.waveNamesNegRefl.append(waveName)
				else:
					pyRootPwa.utils.printErr("Reflectivity cannot be '{0}' for wave '{1}'!".format(refl, waveName))
					raise Exception()
			else:
				pyRootPwa.utils.printErr("Cannot build decay topology for wave '{0}'".format(waveName))
				raise Exception()

		self.waveNamesFixed = [waveName for waveName, _, isActive in waveDescriptionActive if not isActive ]

		# here we have to distinguish negatives from positives
		referenceWaves = []
		for tag, waveNames, rank in zip(['pos', 'neg'], [self.waveNamesPosRefl, self.waveNamesNegRefl], [self.rankPosRefl, self.rankNegRefl]):
			if not waveNames:
				continue
			for iRank in xrange(rank):
				self.wavesInSectors.append(waveNames[iRank:])
				if referenceWavesDefinitions is not None:
					referenceWaves.append(findReferenceWave(referenceWavesDefinitions, tag, iRank))
					if referenceWaves[-1] not in self.wavesInSectors[-1]:
						pyRootPwa.utils.printErr("Reference wave '{0}' not in list of waves of reflectivity '{1}' and rank '{2}'.".format(referenceWaves[-1], tag, iRank))
						raise Exception()
				else:
					referenceWaves.append(self.wavesInSectors[-1][0])
		return referenceWaves

	def amplitudeNames(self):
		amplitudeNames = []
		for iSector, wavesInSector in enumerate(self.wavesInSectors):
			for wave in wavesInSector:
				if wave != 'flat':
					name = "V{s}_{w}".format(s=iSector, w=wave)
				else:
					name = "V_flat"
				amplitudeNames.append(name)
		return amplitudeNames


	def getNormSubmatrices(self):
		matrices = []
		if self.integralsPosRefl:
			matrices += [np.sum(np.dstack([i[0] for i in self.integralsPosRefl]),axis=2)/len(self.integralsPosRefl)]
		if self.integralsNegRefl:
			matrices += [np.sum(np.dstack([i[0] for i in self.integralsNegRefl]),axis=2)/len(self.integralsNegRefl)]
		return matrices + [np.sum(np.dstack([i[-1] for i in self.likelihood.normMatrices]), axis=2)/len(self.likelihood.normMatrices)]


	def getAccSubmatrices(self):
		matrices = []
		if self.integralsPosRefl:
			matrices += [np.sum(np.dstack([i[1] for i in self.integralsPosRefl]),axis=2)/len(self.integralsPosRefl)]
		if self.integralsNegRefl:
			matrices += [np.sum(np.dstack([i[1] for i in self.integralsNegRefl]),axis=2)/len(self.integralsNegRefl)]
		matrices = matrices + [np.sum(np.dstack([i[-1] for i in self.likelihood.accMatrices]), axis=2)/len(self.likelihood.accMatrices)]
		if self.nmbDatasets > 1:
			pyRootPwa.utils.printWarn("Return empty accept matrices because of multiple data sets!")
			matrices = [0*m for m in matrices]
		return matrices


	def getNormIntegrals(self):
		matrices = []
		if self.integralsPosRefl:
			matrices += [np.sum(np.dstack([i[2] for i in self.integralsPosRefl]),axis=2)]
		if self.integralsNegRefl:
			matrices += [np.sum(np.dstack([i[2] for i in self.integralsNegRefl]),axis=2)]
		matrices = matrices + [np.sum(np.dstack([i[-1] for i in self.likelihood.normIntegrals]), axis=2)]
		matrices = [m.reshape(m.size) for m in matrices]
		return matrices


def findReferenceWave(referenceWavesDefinitions, refl, rank):
	'''
	Find the reference wave from the given list of reference waves for the given rank and reflectivity.
	@return: Wave name of reference wave
	'''
	foundReferences = []
	for ref in referenceWavesDefinitions:
		status, topology = ref['description'].constructDecayTopology()
		if status:
			refRefl = topology.XParticle().reflectivity
			if refRefl > 0:
				refRefl = 'pos'
			elif refRefl < 0:
				refRefl = 'neg'
			else:
				pyRootPwa.utils.printErr("Reflectivity cannot be '{0}' for wave '{1}'!".format(refRefl, ref['name']))
				raise Exception()
		if refRefl == refl:
			if 'rank' not in refl or refl['rank'] == rank:
				foundReferences.append(ref['name'])
	if not foundReferences:
		pyRootPwa.utils.printErr("Cannot find reference wave for {0} reflectivity and rank {1}!".format(refl, rank))
		raise Exception()
	elif len(foundReferences) > 1:
		pyRootPwa.utils.printErr("Found {2} reference waves for {0} reflectivity and rank {1}!".format(refl, rank, len(foundReferences)))
		raise Exception()

	return foundReferences[0]


def loadMatrices(normIntegralFileName, accIntegralFileName, waveNames):
	'''
	@return: Normalized norm and acceptance-corrected integral matrices for the given waves
	'''
	pyRootPwa.utils.printInfo("Load integral matrices")
	# load matrices
	normIntegralFile = pyRootPwa.ROOT.TFile.Open(normIntegralFileName)
	accIntegralFile = pyRootPwa.ROOT.TFile.Open(accIntegralFileName)

	normIntegralMeta = pyRootPwa.core.ampIntegralMatrixMetadata.readIntegralFile(normIntegralFile)
	accIntegralMeta = pyRootPwa.core.ampIntegralMatrixMetadata.readIntegralFile(accIntegralFile)

	normIntegralMatrixRpwa = normIntegralMeta.getAmpIntegralMatrix()
	accIntegralMatrixRpwa = accIntegralMeta.getAmpIntegralMatrix()

	totAcc = accIntegralMatrixRpwa.nmbEvents() * 1. / normIntegralMatrixRpwa.nmbEvents()

	normIntegralMatrix = np.empty((len(waveNames), len(waveNames)), dtype=np.complex128)
	accIntegralMatrix = np.empty((len(waveNames), len(waveNames)), dtype=np.complex128)

	for i, wavei in enumerate(waveNames):
		iMatrix = normIntegralMatrixRpwa.waveIndex(wavei)
		for j, wavej in enumerate(waveNames):
			jMatrix = normIntegralMatrixRpwa.waveIndex(wavej)
			normIntegralMatrix[i, j] = normIntegralMatrixRpwa.element(iMatrix, jMatrix)
	for i, wavei in enumerate(waveNames):
		iMatrix = accIntegralMatrixRpwa.waveIndex(wavei)
		for j, wavej in enumerate(waveNames):
			jMatrix = accIntegralMatrixRpwa.waveIndex(wavej)
			accIntegralMatrix[i, j] = accIntegralMatrixRpwa.element(iMatrix, jMatrix)

	normIntegralFile.Close()
	accIntegralFile.Close()

	# get PS matrix diag elements
	normIntegrals = np.real(np.copy(normIntegralMatrix.diagonal()))

	# normalize matrices
	normIntegralMatrix = normIntegralMatrix / np.outer(np.sqrt(normIntegrals), np.sqrt(normIntegrals))
	accIntegralMatrix = totAcc * accIntegralMatrix / np.outer(np.sqrt(normIntegrals), np.sqrt(normIntegrals))

	return normIntegralMatrix, accIntegralMatrix, normIntegrals, totAcc


def buildSubMatrix(fullMatrix, waveNamesAll, waveNames):
	indices = [waveNamesAll.index(w) for w in waveNames]
	return np.copy(fullMatrix[:, indices][indices, :])


def buildSubList(fullList, waveNamesAll, waveNames):
	indices = [waveNamesAll.index(w) for w in waveNames]
	return np.copy(fullList[indices])


def loadAmplitudes(eventAndAmpFileDict, waveNames, multibin, normIntegrals=None):
	'''
	@param normIntegrals: If given, the decay amplitudes will be normalized to normIntegrals
	'''
	pyRootPwa.utils.printInfo("Load decay amplitudes using on-the-fly bin:")
	pyRootPwa.utils.printInfo("\t" + str(multibin))
	amps = []
	for eventFileName, amplitudeFilenames in eventAndAmpFileDict.iteritems():
		ampsInEventfile = pyRootPwa.core.loadAmplitudes([amplitudeFilenames[w] for w in waveNames], waveNames, eventFileName, multibin.boundaries)
		amps.append(ampsInEventfile)
	amps = np.hstack(amps)

	if normIntegrals is not None:
		for i in xrange(amps.shape[0]):
			amps[i, :] /= np.sqrt(normIntegrals[i])

	return amps


class ModelConnected(Model):


	def __init__(self,clsLikelihood = LikelihoodConnected, clsParameterMapping = ParameterMappingConnected, clsLikelihoodInBin = Likelihood, clsParameterMappingInBin = ParameterMappingRpwa):
		'''
		@param clsLikelihoodInBin: Likelihood class for the likelihood of the individual bins
		@param clsParameterMappingInBin: Parameter-mapping class for the parameter mappint in the individual bins
		'''
		Model.__init__(self, clsLikelihood, clsParameterMapping)
		self.clsLikelihoodInBin = clsLikelihoodInBin
		self.clsParameterMappingInBin = clsParameterMappingInBin
		self.models = None # list of cells in the individual bins


	def initModelInBins(self, fileManagerOrConfigfile, binIDs, waveListFileName, rankPosRefl, rankNegRefl, datasets = None):

		if isinstance(fileManagerOrConfigfile, str):
			config = pyRootPwa.rootPwaConfig()
			if not config.initialize(fileManagerOrConfigfile):
				pyRootPwa.utils.printErr("loading config file '" + fileManagerOrConfigfile + "' failed. Aborting...")
				raise Exception()
			pyRootPwa.core.particleDataTable.readFile(config.pdgFileName)
			fileManager = pyRootPwa.loadFileManager(config.fileManagerPath)
			if not fileManager:
				pyRootPwa.utils.printErr("loading the file manager failed. Aborting...")
				raise Exception()
		else:
			fileManager = fileManagerOrConfigfile

		self.models = []
		for binID in binIDs:
			print "Reading in bin: " + str(binID)
			model = ModelRpwa(self.clsLikelihoodInBin, self.clsParameterMappingInBin)
			model.initModelInBin(fileManager, binID, waveListFileName, rankPosRefl, rankNegRefl, datasets=datasets)
			self.models.append(model)

		binWidths = np.array([model.multibin.getBinWidths()['mass'] for model in self.models])
		self.likelihood = self.clsLikelihood([model.likelihood for model in self.models], binWidths)
		self.parameterMapping = self.clsParameterMapping(self)
		self.likelihood.parameterMapping = self.parameterMapping


	def amplitudeNames(self):
		raise Exception("Must be implemented")
