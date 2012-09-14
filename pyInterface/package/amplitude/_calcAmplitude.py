
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
			inFile.readTree()
			nEntries = len(inFile)
			if not pythonAdmin.initKinematicsData(inFile.prodKinParticles, inFile.decayKinParticles):
				pyRootPwa.utils.printErr('Could not initialize kinematics Data "' + keyfile + '".')
				return False

		progressbar = pyRootPwa.utils.progressBar(0, nEntries)
		progressbar.start()
		try:
			treeIndex = 0
			while treeIndex < nEntries:
				with inFile.lock:
					data = inFile[treeIndex]
				treeIndex += 1
				if not pythonAdmin.readKinematicsData(data[0], data[1]):
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

#				for cacheIndex in range(nTreeEntriesToCache):
#					with inFile.lock:
##						print(inFile.tree.GetEntry(treeIndex))
##						prodKinMomentaCache = pyRootPwa.ROOT.TClonesArray(prodKinMomenta)
##						decayKinMomentaCache = pyRootPwa.ROOT.TClonesArray(decayKinMomenta)
#						data = inFile[treeIndex]
#					treeIndex += 1
#					data[0].Print()
#					dataCache.append(data)
#
#				while len(dataCache) > 0:
#					data = dataCache.pop(0)
#					if not pythonAdmin.readKinematicsData(data[0], data[1]):
#						progressbar.cancel()
#						pyRootPwa.utils.printErr('Could not read kinematics data.')
#						return False
#					amp = pythonAdmin()
#					if outputFileFormat == "ascii":
#						outFile.write("(" + str(amp.real) + "," + str(amp.imag) + ")\n")
#					elif outputFileFormat == "binary":
#						arrayAmp = array.array('d', [amp.real, amp.imag])
#						arrayAmp.tofile(outFile)
#					elif writeRootFile:
#						amplitudeTreeLeaf.setAmp(amp)
#						outTree.Fill()
#					else:
#						raise Exception('Something is wrong, this should have been checked in the initialization of the configuration!')
#				progressbar.update(treeIndex)
		except:
			progressbar.cancel()
			raise

		if writeRootFile:
			outTree.Write()
			outFile.Close()

		return True

