'''
@author: F. Kaspar, S. Wallner
'''

import numpy as np
from numpy.random import RandomState

class StartParameterGenerator(object):
	'''
	Generic class to generate start parameters
	'''


	def __init__(self, model, seed = None):
		'''
		Constructor
		'''
		self.generator = RandomState(seed)
		self.model = model
		self.nmbEvents = np.sum(self.model.likelihood.nmbEvents)
		self.approxAcceptance = self.model.likelihood.accMatrices[0][-1][0,0].real


	def __call__(self):
		raise Exception("Call operator needs to be implemented")



class StartParameterGeneratorRpwaUniform(StartParameterGenerator):
	'''
	Draw amplitudes uniformly in the range [-sqrt(nevents), sqrt(nevents)]
	'''
	def __init__(self, model, seed=None):
		StartParameterGenerator.__init__(self, model, seed=seed)
		self.scaleToNbar = True


	def __call__(self):
		randomRatios=False
		parameters = []
		amplVectors = []
		amplMax = self.nmbEvents**0.5 / self.approxAcceptance**0.5
		for nWavesInSector in self.model.parameterMapping.nmbWavesInSectors:
			ampl = np.random.uniform(0.01, amplMax, nWavesInSector)*(2*np.random.randint(0,2, size=nWavesInSector)-1) + \
			       1j*np.random.uniform(0.01, amplMax, nWavesInSector)*(2*np.random.randint(0,2, size=nWavesInSector)-1)
			amplVectors.append(ampl)
		# this sets the amplitues of real-valued waves (e.g. reference wave) to real
		paraNegLlhd = [amplVectors, np.array([1]*self.model.nmbDatasets)] + [1]*self.model.parameterMapping.nmbAdditionalParameters
		paraFitter = self.model.parameterMapping.paraLlhd2Fitter(self.model.parameterMapping.paraNegLlhd2Llhd(paraNegLlhd))
		paraNegLlh = self.model.parameterMapping.paraLlhd2negLlhd(self.model.parameterMapping.paraFitter2Llhd(paraFitter))
		amplVectors = paraNegLlh[0]

		if not self.scaleToNbar:
			# add ratio parameters
			ratios = np.empty(self.model.nmbDatasets)
			if randomRatios:
				tot = 1.0
				for i in xrange(ratios.size-1):
					ratios[i] = self.generator.rand()*tot # [0, tot)
					tot -= ratios[i]
				ratios[-1] = tot
			else: # from data
				nmbEvents = self.model.likelihood.nmbEvents
				ratios[:] = nmbEvents.astype(ratios.dtype)/np.sum(nmbEvents)
		else:
			nBars = [] # not the actual n-bar, but estimates using the generated amplitudes (without correct scaling)
			for iDataset in xrange(self.model.likelihood.nmbDatasets):
				nBar = 0.0
				for iSector in xrange(self.model.likelihood.nmbSectors):
					nBar = nBar + np.real(np.dot( np.dot( self.model.likelihood.accMatrices[iDataset][iSector],np.conj(amplVectors[iSector]) ), amplVectors[iSector]))
				nBars.append(nBar)
			nBars = np.array(nBars)
			scalings = self.model.likelihood.nmbEvents/nBars # scaling factors sucht that nBar[i] == nmbEvents[i]
			ratios = scalings/np.sum(scalings) # using nmbEvents[i] = r[i] * amplitudeScaling^2 * nBar[i]
			amplitudeScaling = np.sqrt(np.sum(scalings))
			amplVectors = [ amplitudeScaling*v for v in amplVectors]


		parameters += [amplVectors]
		datasetRatioParameters = self.model.parameterMapping.datasetDatasetRatios2RatioParameters(ratios)
		parameters.append(datasetRatioParameters)

		# add additional parameters
		for _ in xrange(self.model.parameterMapping.nmbAdditionalParameters):
			parameters.append(self.generator.rand())

		return self.model.parameterMapping.paraLlhd2Fitter(self.model.parameterMapping.paraNegLlhd2Llhd(parameters))


class StartParameterGeneratorRpwaEllipsoid(StartParameterGenerator):
	'''
	Draw amplitudes uniformly on the surface of a d-dimensional hyperellipsoid in the amplitude space.
	Start values for additional parameters are drawn from a uniform distribution in [0,1)

	The transofmation algorithm was taken from here: https://stackoverflow.com/questions/5408276/sampling-uniformly-distributed-random-points-inside-a-spherical-volume
	A detailed explanation can be found here: http://www-alg.ist.hokudai.ac.jp/~jan/randsphere.pdf
	'''
	def __init__(self, model, seed=None):
		StartParameterGenerator.__init__(self, model, seed=seed)

		self.accMatrices = self.model.likelihood.accMatrices

		# build real-valued integral matrix over all sectors
		realMatrices = []
		self.totNmbRealAmplitudes = 0
		for accMatrix in self.accMatrices:
			realAccMatrixTemp = np.zeros(shape=(2*accMatrix.shape[0],2*accMatrix.shape[0]),dtype=np.float64)
			realAccMatrixTemp[:-1:2,:] = accMatrix.view(np.float64)
			realAccMatrixTemp[1::2,:] = (1.j*accMatrix).view(np.float64)

			realMatrices.append(realAccMatrixTemp)
			self.totNmbRealAmplitudes += realAccMatrixTemp.shape[0]
		self.realAccMatrix = np.zeros(shape=(self.totNmbRealAmplitudes,self.totNmbRealAmplitudes),dtype=np.float64)
		offset = 0
		for realMatrix in realMatrices:
			self.realAccMatrix[offset:offset+realMatrix.shape[0],offset:offset+realMatrix.shape[0]] = realMatrix
			offset += realMatrix.shape[0]

		self.inverseRealAccMatrix = np.linalg.inv(self.realAccMatrix)
		# make matrix symmetric by construction
		self.inverseRealAccMatrix = 0.5*(self.inverseRealAccMatrix+self.inverseRealAccMatrix.T)


	def __call__(self):
		realAmplVector = np.random.multivariate_normal(np.zeros(self.totNmbRealAmplitudes), self.inverseRealAccMatrix)
		realAmplVectorLength = np.sqrt(np.dot(realAmplVector,np.dot(self.realAccMatrix,realAmplVector)))

		# normalize amplitude to the total intensity
		realAmplVectorNormalized = np.sqrt(self.nmbEvents)*realAmplVector/realAmplVectorLength

		# transform to complex-valued arrays
		amplVectors = []
		offset = 0
		for nWavesInSector in [accMatrix.shape[0] for accMatrix in self.accMatrices]:
			amplVector =  realAmplVectorNormalized[offset:offset+nWavesInSector*2:2] + 1j*realAmplVectorNormalized[offset+1:offset+nWavesInSector*2+1:2]
			offset += nWavesInSector*2
			# rotate phase to zero
			amplVector = amplVector * np.conj(amplVector[0])/np.abs(amplVector[0])
			amplVectors.append(amplVector)

		parameters = [amplVectors]
		# add ratio parameters
		ratios = np.empty(self.model.nmbDatasets)
		tot = 1.0
		for i in xrange(ratios.size):
			ratios[i] = tot
			tot -= self.generator.rand()*tot # [0, tot)
		ratios = np.array([0.5, 0.25, 0.25])
		parameters.append(ratios)

		# add additional parameters
		for _ in xrange(self.model.parameterMapping.nmbAdditionalParameters):
			parameters.append(self.generator.rand())

		return self.model.parameterMapping.paraLlhd2Fitter(self.model.parameterMapping.paraNegLlhd2Llhd(parameters))


# this is gonna be tricky to realize!!!
# maybe sampling in fourier-space works better
class StartParameterGeneratorUniform(object):


	def __init__(self, model, seed=None):
		self.generator = RandomState(seed)
		self.model = model
		self.startParameterGenerators = [StartParameterGeneratorRpwaUniform(subbinModel, seed=1) for subbinModel in self.model.models]
		for spg in self.startParameterGenerators:
			spg.generator = self.generator


	def __call__(self):
		pars = []

		# there is probably a better way to do this ...
		for spg in self.startParameterGenerators:
			pars += spg().tolist()

		return np.array(pars)
