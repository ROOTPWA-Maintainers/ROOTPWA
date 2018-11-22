'''
@author: F. Kaspar, S. Wallner
'''

import numpy as np
import pyRootPwa
ROOT = pyRootPwa.ROOT


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


	def amplitudeNames(self):
		raise NotImplementedError("Must be implemented in specialized class")



class ModelRpwa(Model):
	def __init__(self, clsLikelihood, clsParameterMapping):
		Model.__init__(self, clsLikelihood, clsParameterMapping)
		self.rank = None
		self.multibin = None

	def initModelInBin(self, fileManagerOrConfigfile, multiBin, waveListFileName, rank):
		'''
		@param fileManagerOrConfigfile: Can be a fileManager object or a path to the config file
		@param multiBin: Can be a multibin or the integer bin-id
		@param waveListFileName: Path to the wavelist file
		@param rank: Rank
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
			fileManager = fileManager

		if isinstance(multiBin, int):
			multiBin = fileManager.binList[multiBin]

		self.initModelWithFiles(eventAndAmpFileDict = fileManager.getEventAndAmplitudeFilePathsInBin(multiBin, pyRootPwa.core.eventMetadata.REAL),
		                        normIntegralFileName = fileManager.getIntegralFilePath(multiBin, pyRootPwa.core.eventMetadata.GENERATED),
		                        accIntegralFileName = fileManager.getIntegralFilePath(multiBin, pyRootPwa.core.eventMetadata.ACCEPTED),
		                        multiBin = multiBin,
		                        waveListFileName = waveListFileName,
		                        waveDescriptions = fileManager.getWaveDescriptions(),
		                        rank = rank)



	def initModelWithFiles(self,
	                       eventAndAmpFileDict,
	                       normIntegralFileName, accIntegralFileName,
	                       multiBin,
	                       waveListFileName, waveDescriptions,
	                       rank):

		self.rank = rank
		self.multibin = multiBin

		waveDescThres = pyRootPwa.utils.getWaveDescThresFromWaveList(waveListFileName, waveDescriptions)

		self.waveNames = [waveName for waveName, _, _ in waveDescThres]

		normIntegralMatrixFull, accIntegralMatrixFull, normIntegralsFull, totAcc = loadMatrices(normIntegralFileName, accIntegralFileName, self.waveNames)

		decayAmpsFull = loadAmplitudes(eventAndAmpFileDict, self.waveNames, multiBin, normIntegralsFull)

		# here we have to distinguish negatives from positives
		normIntegralMatrices = []
		accIntegralMatrices = []
		normIntegrals = []
		decayAmps = []
		for iRank in xrange(rank):
			self.wavesInSectors.append(self.waveNames[iRank:])
			decayAmps.append(decayAmpsFull[iRank:])
			normIntegralMatrices.append(normIntegralMatrixFull[iRank:,iRank:])
			accIntegralMatrices.append( accIntegralMatrixFull[iRank:,iRank:])
			normIntegrals.append(       normIntegralsFull[iRank:])
		# add flat wave
		self.wavesInSectors.append(["flat"])
		decayAmps.append(np.ones((1,decayAmpsFull.shape[1]), dtype=np.complex128))
		normIntegralMatrices.append(np.ones((1,1), dtype=np.complex128))
		accIntegralMatrices.append(np.full((1,1), totAcc, dtype=np.complex128))
		normIntegrals.append(      np.ones((1), dtype=np.float64))

		self.waveNames.append("flat")
		self.nmbWaves = len(self.waveNames)


		# buildParameterMapping
		self.parameterMapping = self.clsParameterMapping(self)

		self.likelihood = self.clsLikelihood(decayAmps, accIntegralMatrices, normIntegralMatrices,normIntegrals, self.parameterMapping)


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
		# TODO: Add negatives here
		return [self.likelihood.normMatrices[0], self.likelihood.normMatrices[-1]]

	def getAccSubmatrices(self):
		# TODO: Add negatives here
		return [self.likelihood.accMatrices[0], self.likelihood.accMatrices[-1]]

	def getNormIntegrals(self):
		# TODO: Add negatives here
		return [self.likelihood.normIntegrals[0], self.likelihood.normIntegrals[-1]]


def loadMatrices(normIntegralFileName, accIntegralFileName, waveNames):
	'''
	@return: Normalized norm and acceptance-corrected integral matrices for the given waves
	'''
	pyRootPwa.utils.printInfo("Load integral matrices")
	# load matrices
	normIntegralFile      = pyRootPwa.ROOT.TFile.Open(normIntegralFileName)
	accIntegralFile     = pyRootPwa.ROOT.TFile.Open(accIntegralFileName)

	normIntegralMeta  = pyRootPwa.core.ampIntegralMatrixMetadata.readIntegralFile(normIntegralFile)
	accIntegralMeta = pyRootPwa.core.ampIntegralMatrixMetadata.readIntegralFile(accIntegralFile)

	normIntegralMatrixRpwa    = normIntegralMeta.getAmpIntegralMatrix()
	accIntegralMatrixRpwa     = accIntegralMeta.getAmpIntegralMatrix()

	totAcc = accIntegralMatrixRpwa.nmbEvents() * 1. / normIntegralMatrixRpwa.nmbEvents()

	normIntegralMatrix    = np.empty((len(waveNames), len(waveNames)), dtype=np.complex128)
	accIntegralMatrix     = np.empty((len(waveNames), len(waveNames)), dtype=np.complex128)

	for i, wavei in enumerate(waveNames):
		iMatrix = normIntegralMatrixRpwa.waveIndex(wavei)
		for j, wavej in enumerate(waveNames):
			jMatrix = normIntegralMatrixRpwa.waveIndex(wavej)
			normIntegralMatrix[i,j] = normIntegralMatrixRpwa.element(iMatrix, jMatrix)
	for i, wavei in enumerate(waveNames):
		iMatrix = accIntegralMatrixRpwa.waveIndex(wavei)
		for j,wavej in enumerate(waveNames):
			jMatrix = accIntegralMatrixRpwa.waveIndex(wavej)
			accIntegralMatrix[i,j] = accIntegralMatrixRpwa.element(iMatrix, jMatrix)

	normIntegralFile.Close()
	accIntegralFile.Close()

	# get PS matrix diag elements
	normIntegrals = np.real(np.copy(normIntegralMatrix.diagonal()))

	# normalize matrices
	normIntegralMatrix = normIntegralMatrix/np.outer(np.sqrt(normIntegrals), np.sqrt(normIntegrals))
	accIntegralMatrix = totAcc*accIntegralMatrix/np.outer(np.sqrt(normIntegrals), np.sqrt(normIntegrals))

	return normIntegralMatrix, accIntegralMatrix, normIntegrals, totAcc


def loadAmplitudes(eventAndAmpFileDict, waveNames, multibin, normIntegrals = None):
	'''
	@param normIntegrals: If given, the decay amplitudes will be normalized to normIntegrals
	'''
	pyRootPwa.utils.printInfo("Load decay amplitudes using on-the-fly bin:")
	pyRootPwa.utils.printInfo("\t" + str(multibin))
	amps = []
	for eventFileName, amplitudeFilenames in eventAndAmpFileDict.iteritems():
		eventFile, eventMeta = pyRootPwa.utils.openEventFile(eventFileName)
		if not eventFile or not eventMeta:
			pyRootPwa.utils.printErr("could not open event file '" + eventFileName + "'. Aborting...")
			raise Exception()

		ampMetas = []
		ampFiles = []
		for waveName in waveNames:
			ampFileName = amplitudeFilenames[waveName]
			ampFile = ROOT.TFile.Open(ampFileName, "READ")
			if not ampFile:
				pyRootPwa.utils.printErr("could not open amplitude file '" + ampFileName + "'.")
				raise Exception()
			meta = pyRootPwa.core.amplitudeMetadata.readAmplitudeFile(ampFile, waveName)
			if not meta:
				pyRootPwa.utils.printErr("could not get metadata for waveName '" + waveName + "'.")
				raise Exception()
			ampMetas.append(meta)
			ampFiles.append(ampFile)

		ampsInEventfile = pyRootPwa.core.loadAmplitudes(ampMetas, eventMeta, multibin.boundaries)
		amps.append(ampsInEventfile)
		for ampFile in ampFiles:
			ampFile.Close()
	amps = np.hstack(amps)

	if normIntegrals is not None:
		for i in xrange(amps.shape[0]):
			amps[i,:] /= np.sqrt(normIntegrals[i])

	return amps
