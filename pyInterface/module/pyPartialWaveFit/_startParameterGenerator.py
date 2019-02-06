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
		self.nmbEvetns = self.model.likelihood.nmbEvents


	def __call__(self):
		raise Exception("Call operator needs to be implemented")



class StartParameterGeneratorRpwaUniform(StartParameterGenerator):
	'''
	Draw amplitudes uniformly in the range [-sqrt(nevents), sqrt(nevents)]
	'''
	def __init__(self, model, seed=None):
		StartParameterGenerator.__init__(self, model, seed=seed)


	def __call__(self):
		amplVectors = []
		for nWavesInSector in self.model.parameterMapping.nmbWavesInSectors:
			ampl = np.random.uniform(0.01, self.nmbEvetns**0.5, nWavesInSector)*(2*np.random.randint(0,2, size=nWavesInSector)-1) + \
			       1j*np.random.uniform(0.01, self.nmbEvetns**0.5, nWavesInSector)*(2*np.random.randint(0,2, size=nWavesInSector)-1)
			ampl[0] = np.real(ampl[0])
			amplVectors.append(ampl)

		# add additional parameters
		parameters = [amplVectors]
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
		realAmplVectorNormalized = np.sqrt(self.nmbEvetns)*realAmplVector/realAmplVectorLength

		# transform to complex-valued arrays
		amplVectors = []
		offset = 0
		for nWavesInSector in [accMatrix.shape[0] for accMatrix in self.accMatrices]:
			amplVector =  realAmplVectorNormalized[offset:offset+nWavesInSector*2:2] + 1j*realAmplVectorNormalized[offset+1:offset+nWavesInSector*2+1:2]
			offset += nWavesInSector*2
			# rotate phase to zero
			amplVector = amplVector * np.conj(amplVector[0])/np.abs(amplVector[0])
			amplVectors.append(amplVector)

		# add additional parameters
		parameters = [amplVectors]
		for _ in xrange(self.model.parameterMapping.nmbAdditionalParameters):
			parameters.append(self.generator.rand())

		return self.model.parameterMapping.paraLlhd2Fitter(self.model.parameterMapping.paraNegLlhd2Llhd(parameters))


# this is gonna be tricky to realize!!!
# maybe sampling in fourier-space works better
class StartParameterGeneratorUniform(StartParameterGenerator):
    def __init__(self, model, seed=None):
       StartParameterGenerator.__init__(self, model, seed=seed)
       self.startParameterGenerators = [StartParameterGeneratorRpwaUniform(m,seed=int(self.generator.rand()*999999)) for m in self.model.likelihood.models]   
        
    def __call__(self):
        pars = []

        # there is probably a better way to do this ...
        for spg in self.startParameterGenerators:
            pars += spg().tolist()

        return np.array(pars)

