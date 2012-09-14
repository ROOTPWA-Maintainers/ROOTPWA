
import multiprocessing

import pyRootPwa
import pyRootPwa.exception
import pyRootPwa.utils

class inputFile():

	lock = multiprocessing.Lock()

	tree = None
	prodKinParticles = None
	decayKinParticles = None

	treeEntries = []

	_inFileName = ""

	def __init__(self, inFileName):
		self._inFileName = inFileName
		if pyRootPwa.config is None:
			raise pyRootPwa.exception.pyRootPwaException("pyRootPwa configuration not initialized")
		if inFileName.endswith('.root'):
			self._inFile = pyRootPwa.ROOT.TFile.Open(inFileName)
			with pyRootPwa.config.lock:
				self.prodKinParticles = self._inFile.Get(pyRootPwa.config.prodKinPartNamesObjName)
				self.decayKinParticles = self._inFile.Get(pyRootPwa.config.decayKinPartNamesObjName)
			self._inFile.Close()
		elif inFileName.endswith('.evt'):
			(self.prodKinParticles, self.decayKinParticles, self.tree) = pyRootPwa.utils.getTreeFromEvtFile(inFileName, inFileName)
		else:
			raise pyRootPwa.exception.pyRootPwaException("Unknown file extension")

	def readTree(self):
		if not self.treeEntries:
			if self._inFileName.endswith('.root'):
				self._inFile = pyRootPwa.ROOT.TFile.Open(self._inFileName)
				with pyRootPwa.config.lock:
					self.tree = self._inFile.Get(pyRootPwa.config.inTreeName)
			prodKinMomenta = pyRootPwa.ROOT.TClonesArray("TVector3")
			decayKinMomenta = pyRootPwa.ROOT.TClonesArray("TVector3")
			with pyRootPwa.config.lock:
				self.tree.SetBranchAddress(pyRootPwa.config.prodKinMomentaLeafName, prodKinMomenta)
				self.tree.SetBranchAddress(pyRootPwa.config.decayKinMomentaLeafName, decayKinMomenta)
			for treeIndex in range(self.tree.GetEntries()):
				self.tree.GetEntry(treeIndex)
				self.treeEntries.append((prodKinMomenta, decayKinMomenta))
		return self.treeEntries

	def __len__(self):
		if not self.treeEntries:
			return -1
		else:
			return len(self.treeEntries)

	def __iter__(self):
		if not self.treeEntries:
			raise pyRootPwa.exception.pyRootPwaException('Trying to iterate over uninitialized inputFile')
		return iter(self.treeEntries)

	def __getitem__(self, index):
		if not self.treeEntries:
			return None
		return self.treeEntries[index]

	def __str__(self):
		retval = ""
		retval += "pyRootPwa.inputFile\n"
		retval += "-------------------\n\n"
		retval += "_inFileName = " + self._inFileName + "\n"
		retval += str(self.tree) + "\n"
		retval += str(self.prodKinParticles) + "\n"
		retval += str(self.decayKinParticles) + "\n\n"
		return retval
