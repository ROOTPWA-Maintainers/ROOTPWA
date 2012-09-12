
import pyRootPwa
import pyRootPwa.exception
import pyRootPwa.utils

class inputFile():

	tree = None
	prodKinParticles = None
	decayKinParticles = None

	_inFileName = ""

	def __init__(self, inFileName):
		self._inFileName = inFileName
		if inFileName.endswith('.root'):
			self._inFile = pyRootPwa.ROOT.TFile.Open(inFileName)
			self.prodKinParticles = self._inFile.Get(pyRootPwa.config.prodKinPartNamesObjName)
			self.decayKinParticles = self._inFile.Get(pyRootPwa.config.decayKinPartNamesObjName)
			self.tree = self._inFile.Get(pyRootPwa.config.inTreeName)
		elif inFileName.endswith('.evt'):
			(self.prodKinParticles, self.decayKinParticles, self.tree) = pyRootPwa.utils.getTreeFromEvtFile(inFileName, inFileName)
		else:
			raise pyRootPwa.exception.pyRootPwaException("Unknown file extension")

