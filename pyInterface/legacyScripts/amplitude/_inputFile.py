
import pyRootPwa
import pyRootPwa.utils

class inputFile():

	inFileName = ""

	tree = None
	prodKinParticles = None
	decayKinParticles = None

	_readRootFile = None
	_entriesInTree = 0
	_inFile = None

	_inWith = False

	def __init__(self, inFileName):
		self.inFileName = inFileName
		if pyRootPwa.config is None:
			raise pyRootPwa.rootPwaException("pyRootPwa configuration not initialized")
		if inFileName.endswith('.root'):
			self._readRootFile = True
		elif inFileName.endswith('.evt'):
			self._readRootFile = False
			(self.prodKinParticles, self.decayKinParticles, self.tree) = pyRootPwa.utils.getTreeFromEvtFile(inFileName, inFileName)
		else:
			raise pyRootPwa.rootPwaException("Unknown file extension")

	def __len__(self):
		if not self._inWith:
			raise pyRootPwa.rootPwaException("Not in with statement")
		return self._entriesInTree

	def __enter__(self):
		if self._readRootFile:
			self._inFile = pyRootPwa.ROOT.TFile.Open(self.inFileName)
			self.prodKinParticles = self._inFile.Get(pyRootPwa.config.prodKinPartNamesObjName)
			self.decayKinParticles = self._inFile.Get(pyRootPwa.config.decayKinPartNamesObjName)
			self.tree = self._inFile.Get(pyRootPwa.config.inTreeName)
			self._entriesInTree = int(self.tree.GetEntries())
		else:
			self._entriesInTree = int(self.tree.GetEntries())
		self._inWith = True
		return self._inFile

	def __exit__(self, type, value, traceback):
		if self._readRootFile:
			self._inFile.Close()
		self._inWith = False

	def __str__(self):
		retval = ""
		retval += "pyRootPwa.inputFile\n"
		retval += "-------------------\n\n"
		retval += "inFileName = " + self.inFileName + "\n"
		retval += str(self.tree) + "\n"
		retval += str(self.prodKinParticles) + "\n"
		retval += str(self.decayKinParticles) + "\n\n"
		return retval
