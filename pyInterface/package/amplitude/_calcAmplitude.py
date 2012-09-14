
import array
import multiprocessing
import sys

import pyRootPwa
import pyRootPwa.exception
import pyRootPwa.utils

class AmplitudeCalculator(multiprocessing.Process):

	queue = None

	def __init__(self, queue):
		self.queue = queue
		multiprocessing.Process.__init__(self)

	def run(self):
		while True:
			inTuple = self.queue.get()
			outFileName = inTuple[2]
			pyRootPwa.utils.printInfo('Calculating amplitutes with inuput file "' + inTuple[0].inFileName + '" and key file "' + str(inTuple[1]) + '".')
			if outFileName.endswith('.amp'):
				with open(outFileName, 'w') as outFile:
					outTuple = (inTuple[0], inTuple[1], outFile)
					self.calcAmplitudes(outTuple)
			else:
				outTuple = (inTuple[0], inTuple[1], outFile)
				self.calcAmplitudes(outTuple)
			self.queue.task_done()

	def calcAmplitudes(self, inTuple):

		inFile = inTuple[0]
		keyfile = inTuple[1]
		outFile = inTuple[2]

		if pyRootPwa.config is None:
			raise pyRootPwa.exception.pyRootPwaException("pyRootPwa configuration not initialized")

		with pyRootPwa.config.lock:
			prodKinMomentaLeafName = pyRootPwa.config.prodKinMomentaLeafName
			decayKinMomentaLeafName = pyRootPwa.config.decayKinMomentaLeafName
			outputFileFormat = pyRootPwa.config.outputFileFormat
			amplitudeLeafName = pyRootPwa.config.amplitudeLeafName
			nTreeEntriesToCache = pyRootPwa.config.nTreeEntriesToCache

		pythonAdmin = pyRootPwa.pythonAdministrator()

		writeRootFile = False
		if outputFileFormat == "root":
			writeRootFile = True

		if writeRootFile:
			outFile.cd()
			ampTreeName = keyfile.rsplit('/',1)[-1].replace('.key', '.amp')
			outTree = pyRootPwa.ROOT.TTree(ampTreeName, ampTreeName)
			amplitudeTreeLeaf = pyRootPwa.amplitudeTreeLeaf()
			pythonAdmin.branch(outTree, amplitudeTreeLeaf, amplitudeLeafName)
		with keyfile.lock:
			if not pythonAdmin.constructAmplitude(str(keyfile)):
				pyRootPwa.utils.printWarn('Could not construct amplitude for keyfile "' + keyfile + '".')
				return False
		sys.stdout.write(str(pythonAdmin))

		prodKinMomenta = pyRootPwa.ROOT.TClonesArray("TVector3")
		decayKinMomenta = pyRootPwa.ROOT.TClonesArray("TVector3")

		with inFile.lock:
			nEntries = len(inFile)
			if not pythonAdmin.initKinematicsData(inFile.prodKinParticles, inFile.decayKinParticles):
				pyRootPwa.utils.printErr('Could not initialize kinematics Data "' + keyfile + '".')
				return False

		progressbar = pyRootPwa.utils.progressBar(0, nEntries)
		progressbar.start()
		try:
			treeIndex = 0
			while treeIndex < nEntries:
				upperBound = treeIndex + nTreeEntriesToCache
				if upperBound > nEntries:
					upperBound = nEntries
				with inFile.lock:
					data = inFile[treeIndex:upperBound]
				treeIndex = upperBound
				for datum in data:
					if not pythonAdmin.readKinematicsData(datum[0], datum[1]):
						progressbar.cancel()
						pyRootPwa.utils.printErr('Could not read kinematics data.')
						return False
					amp = pythonAdmin()
					if outputFileFormat == "ascii":
						outFile.write("(" + str(amp.real) + "," + str(amp.imag) + ")\n")
					elif outputFileFormat == "binary":
						arrayAmp = array.array('d', [amp.real, amp.imag])
						arrayAmp.tofile(outFile)
					elif writeRootFile:
						amplitudeTreeLeaf.setAmp(amp)
						outTree.Fill()
					else:
						raise Exception('Something is wrong, this should have been checked in the initialization of the configuration!')
				progressbar.update(treeIndex)
		except:
			progressbar.cancel()
			raise

		if writeRootFile:
			outTree.Write()
			outFile.Close()

		pyRootPwa.utils.printSucc('Created amplitude file for ' + str(nEntries) + ' events.')

		return True

