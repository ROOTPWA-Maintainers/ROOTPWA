
import array
import multiprocessing
import os
import sys
import traceback

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
			try:
				inTuple = self.queue.get()
				with pyRootPwa.utils.Silencer() as silencer:
					processedEvents = 0
					outFileName = inTuple[2]
					pyRootPwa.utils.printInfo('Calculating amplitutes with inuput file "' + inTuple[0].inFileName + '" and key file "' + str(inTuple[1]) + '".')
					try:
						if outFileName.endswith('.amp'):
							with open(outFileName, 'w') as outFile:
								outTuple = (inTuple[0], inTuple[1], outFile)
								processedEvents = self.calcAmplitudes(outTuple)
						else:
							outFile = pyRootPwa.ROOT.TFile.Open(outFileName, 'RECREATE')
							outTuple = (inTuple[0], inTuple[1], outFile)
							processedEvents = self.calcAmplitudes(outTuple)
					except:
						processedEvents = -1
						if os.path.exists(outFileName):
							os.remove(outFileName)
						raise
				with pyRootPwa.utils.stdoutLock:
					sys.stdout.write(silencer.output)
					if processedEvents > 0:
						pyRootPwa.utils.printSucc('Created amplitude file "' + outFileName + '" with ' + str(processedEvents) + ' events.\n')
					elif processedEvents == 0:
						pyRootPwa.utils.printWarn('Created amplitude file "' + outFileName + '", but wrote 0 events to it.\n')
					else:
						if os.path.exists(outFileName):
							os.remove(outFileName)
						pyRootPwa.utils.printErr('Amplitude calculation failed for input file "' + inputFile + '" and keyfile "' + keyfile + '".\n')
			except KeyboardInterrupt:
				pyRootPwa.utils.printInfo('Process ' + str(self.pid) + ' caught keyboard interrupt. Terminating...')
			except:
				pyRootPwa.utils.printErr('Process ' + str(self.pid) + ' caught exception. Terminating...')
				traceback.print_exc()
				break
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
				return -1
		sys.stdout.write(str(pythonAdmin))

		prodKinMomenta = pyRootPwa.ROOT.TClonesArray("TVector3")
		decayKinMomenta = pyRootPwa.ROOT.TClonesArray("TVector3")

		with inFile.lock:
			nEntries = len(inFile)
			if not pythonAdmin.initKinematicsData(inFile.prodKinParticles, inFile.decayKinParticles):
				pyRootPwa.utils.printErr('Could not initialize kinematics Data "' + keyfile + '".')
				return -1

		progressbar = pyRootPwa.utils.progressBar(0, nEntries, sys.stdout)
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
						return -1
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

		return nEntries

