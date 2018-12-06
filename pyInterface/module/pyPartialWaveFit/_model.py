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
		self.rankPosRefl = None
		self.rankNegRefl = None
		self.waveNames = None
		self.waveNamesPosRefl = None
		self.waveNamesNegRefl = None
		self.waveNamesFixed = None
		self.integralsPosRefl = None # tuple of (normIntegralMatrix, accIntegralMatrix, normIntegrals)
		self.integralsNegRefl = None # tuple of (normIntegralMatrix, accIntegralMatrix, normIntegrals)
		self.multibin = None

	def initModelInBin(self, fileManagerOrConfigfile, multiBin, waveListFileName, rankPosRefl, rankNegRefl):
		'''
		@param fileManagerOrConfigfile: Can be a fileManager object or a path to the config file
		@param multiBin: Can be a multibin or the integer bin-id
		@param waveListFileName: Path to the wavelist file
		@param rankPosRefl: Rank of positive reflectivity waves
		@param rankNegRefl: Rank of negative reflectivity waves
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

		self.initModelWithFiles(eventAndAmpFileDict = fileManager.getEventAndAmplitudeFilePathsInBin(multiBin, pyRootPwa.core.eventMetadata.REAL),
		                        normIntegralFileName = fileManager.getIntegralFilePath(multiBin, pyRootPwa.core.eventMetadata.GENERATED),
		                        accIntegralFileName = fileManager.getIntegralFilePath(multiBin, pyRootPwa.core.eventMetadata.ACCEPTED),
		                        multiBin = multiBin,
		                        waveListFileName = waveListFileName,
		                        waveDescriptions = fileManager.getWaveDescriptions(),
		                        rankPosRefl = rankPosRefl,
		                        rankNegRefl = rankNegRefl)



	def initModelWithFiles(self,
	                       eventAndAmpFileDict,
	                       normIntegralFileName, accIntegralFileName,
	                       multiBin,
	                       waveListFileName, waveDescriptions,
	                       rankPosRefl,
	                       rankNegRefl):

		self.rankPosRefl = rankPosRefl
		self.rankNegRefl = rankNegRefl
		self.multibin = multiBin

		waveDescriptionActive = pyRootPwa.utils.getWaveDescriptionActiveFromWavelist(waveListFileName, waveDescriptions, multiBin)

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
		normIntegralMatrices = []
		accIntegralMatrices = []
		normIntegrals = []
		decayAmps = []

		normIntegralMatrixFull, accIntegralMatrixFull, normIntegralsFull, totAcc = loadMatrices(normIntegralFileName, accIntegralFileName, self.waveNames)

		for tag, waveNames, rank in zip(['pos', 'neg'], [self.waveNamesPosRefl, self.waveNamesNegRefl], [self.rankPosRefl, self.rankNegRefl]):
			if not waveNames:
				continue
			integrals = (buildSubMatrix(normIntegralMatrixFull, self.waveNames, waveNames),
			            buildSubMatrix(accIntegralMatrixFull,   self.waveNames, waveNames),
			            buildSubList(  normIntegralsFull,       self.waveNames, waveNames))
			decayAmpsFull = loadAmplitudes(eventAndAmpFileDict, waveNames, multiBin, integrals[2])
			if tag == 'pos':
				self.integralsPosRefl = integrals
			if tag == 'neg':
				self.integralsNegRefl = integrals

			for iRank in xrange(rank):
				self.wavesInSectors.append(waveNames[iRank:])
				decayAmps.append(decayAmpsFull[iRank:])
				normIntegralMatrices.append(integrals[0][iRank:,iRank:])
				accIntegralMatrices.append( integrals[1][iRank:,iRank:])
				normIntegrals.append(       integrals[2][iRank:])

		# add flat wave
		self.wavesInSectors.append(["flat"])
		decayAmps.append(np.ones((1,1), dtype=np.complex128))
		normIntegralMatrices.append(np.ones((1,1), dtype=np.complex128))
		accIntegralMatrices.append(np.full((1,1), totAcc, dtype=np.complex128))
		normIntegrals.append(      np.ones((1), dtype=np.float64))

		self.waveNames.append("flat")
		self.nmbWaves = len(self.waveNames)


		# buildParameterMapping
		self.parameterMapping = self.clsParameterMapping(self, self.waveNamesFixed)

		self.likelihood = self.clsLikelihood(decayAmps, accIntegralMatrices, normIntegralMatrices,normIntegrals, self.parameterMapping)

		return True


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
			matrices += [self.integralsPosRefl[0]]
		if self.integralsNegRefl:
			matrices += [self.integralsNegRefl[0]]
		return matrices + [self.likelihood.normMatrices[-1]]

	def getAccSubmatrices(self):
		matrices = []
		if self.integralsPosRefl:
			matrices += [self.integralsPosRefl[1]]
		if self.integralsNegRefl:
			matrices += [self.integralsNegRefl[1]]
		return matrices + [self.likelihood.accMatrices[-1]]

	def getNormIntegrals(self):
		matrices = []
		if self.integralsPosRefl:
			matrices += [self.integralsPosRefl[2]]
		if self.integralsNegRefl:
			matrices += [self.integralsNegRefl[2]]
		return matrices + [self.likelihood.normIntegrals[-1]]


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

def buildSubMatrix(fullMatrix, waveNamesAll, waveNames):
	indices = [waveNamesAll.index(w) for w in waveNames]
	return np.copy(fullMatrix[:,indices][indices,:])

def buildSubList(fullList, waveNamesAll, waveNames):
	indices = [waveNamesAll.index(w) for w in waveNames]
	return np.copy(fullList[indices])


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
