
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

	inFileName = ""
	_prodKinMomentaLeafName = ""
	_decayKinMomentaLeafName = ""
	_prodKinPartNamesObjName = ""
	_decayKinPartNamesObjName = ""
	_inTreeName = ""
	_readRootFile = None
	_entriesInTree = 0

	def __init__(self, inFileName):
		self.inFileName = inFileName
		if pyRootPwa.config is None:
			raise pyRootPwa.exception.pyRootPwaException("pyRootPwa configuration not initialized")
		with pyRootPwa.config.lock:
			self._prodKinMomentaLeafName =  pyRootPwa.config.prodKinMomentaLeafName
			self._decayKinMomentaLeafName = pyRootPwa.config.decayKinMomentaLeafName
			self._prodKinPartNamesObjName = pyRootPwa.config.prodKinPartNamesObjName
			self._decayKinPartNamesObjName = pyRootPwa.config.decayKinPartNamesObjName
			self._inTreeName = pyRootPwa.config.inTreeName
		if inFileName.endswith('.root'):
			self._readRootFile = True
			self._inFile = pyRootPwa.ROOT.TFile.Open(inFileName)
			self.prodKinParticles = self._inFile.Get(self._prodKinPartNamesObjName)
			self.decayKinParticles = self._inFile.Get(self._decayKinPartNamesObjName)
			self._entriesInTree = int(self._inFile.Get(self._inTreeName).GetEntries())
			self._inFile.Close()
		elif inFileName.endswith('.evt'):
			(self.prodKinParticles, self.decayKinParticles, self.tree) = pyRootPwa.utils.getTreeFromEvtFile(inFileName, inFileName)
			prodKinMomenta = pyRootPwa.ROOT.TClonesArray("TVector3")
			decayKinMomenta = pyRootPwa.ROOT.TClonesArray("TVector3")
			self.tree.SetBranchAddress(self._prodKinMomentaLeafName, prodKinMomenta)
			self.tree.SetBranchAddress(self._decayKinMomentaLeafName, decayKinMomenta)
			for treeIndex in range(self.tree.GetEntries()):
				self.tree.GetEntry(treeIndex)
				self.treeEntries.append((prodKinMomenta, decayKinMomenta))
			self._entriesInTree = int(self.tree.GetEntries())
		else:
			raise pyRootPwa.exception.pyRootPwaException("Unknown file extension")

	def __len__(self):
		return self._entriesInTree

	def __getitem__(self, index):
		if self._readRootFile:
			if not self._inFile.IsOpen():
				self._inFile = pyRootPwa.ROOT.TFile.Open(self.inFileName)
			self.tree = self._inFile.Get(self._inTreeName)
			prodKinMomenta = pyRootPwa.ROOT.TClonesArray("TVector3")
			decayKinMomenta = pyRootPwa.ROOT.TClonesArray("TVector3")
			self.tree.SetBranchAddress(self._prodKinMomentaLeafName, prodKinMomenta)
			self.tree.SetBranchAddress(self._decayKinMomentaLeafName, decayKinMomenta)
			if isinstance(index, slice):
				retval = []
				for treeIndex in range(index.start, index.stop):
					self.tree.GetEntry(treeIndex)
					retval.append((prodKinMomenta, decayKinMomenta))
				return retval
			else:
				self.tree.GetEntry(index)
				return ((prodKinMomenta, decayKinMomenta))
		else:
			return self.treeEntries[index]

	def __str__(self):
		retval = ""
		retval += "pyRootPwa.inputFile\n"
		retval += "-------------------\n\n"
		retval += "inFileName = " + self.inFileName + "\n"
		retval += str(self.tree) + "\n"
		retval += str(self.prodKinParticles) + "\n"
		retval += str(self.decayKinParticles) + "\n\n"
		return retval
