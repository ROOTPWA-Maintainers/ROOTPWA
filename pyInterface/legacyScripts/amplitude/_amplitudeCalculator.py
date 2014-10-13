
import array
import multiprocessing
import os
import sys
import traceback
import Queue

import pyRootPwa.core
import pyRootPwa.utils

class amplitudeCalculator(multiprocessing.Process):

	queue = None
	silence = True
	progressBar = True
	maxNmbEvents = -1

	def __init__(self, queue, silence, progressBar, maxNmbEvents):
		self.queue = queue
		self.silence = silence
		self.progressBar = progressBar
		self.maxNmbEvents = maxNmbEvents
		multiprocessing.Process.__init__(self)

	def run(self, runDirectly=False):
		while True:
			try:
				if not runDirectly:
					inTuple = self.queue.get()
				else:
					# This is needed for the profiler
					try:
						inTuple = self.queue.get(True, 5)
					except Queue.Empty:
						pyRootPwa.utils.printInfo('Queue seems to be empty and I was called directly. Terminating now...')
						break
				with pyRootPwa.utils.silencer(self.silence) as silencer:
					processedEvents = 0
					outFileName = inTuple[2]
					pyRootPwa.utils.printInfo('Calculating amplitudes with input file "' + inTuple[0].inFileName + '" and key file "' + str(inTuple[1]) + '".')
					try:
						with inTuple[0]:
							if outFileName.endswith('.amp'):
								with open(outFileName, 'w') as outFile:
									outTuple = (inTuple[0], inTuple[1], outFile)
									processedEvents = self.calcAmplitudes(outTuple)
							else:
								outFile = pyRootPwa.ROOT.TFile.Open(outFileName, 'RECREATE')
								outTuple = (inTuple[0], inTuple[1], outFile)
								processedEvents = self.calcAmplitudes(outTuple)
								outFile.Close()
					except:
						processedEvents = -1
						if os.path.exists(outFileName):
							os.remove(outFileName)
						raise
				sys.stdout.write(silencer.output)
				if processedEvents > 0:
					pyRootPwa.utils.printSucc('Created amplitude file "' + outFileName + '" with ' + str(processedEvents) + ' events.\n')
				elif processedEvents == 0:
					pyRootPwa.utils.printWarn('Created amplitude file "' + outFileName + '", but wrote 0 events to it.\n')
				else:
					if os.path.exists(outFileName):
						os.remove(outFileName)
					pyRootPwa.utils.printErr('Amplitude calculation failed for input file "' + inTuple[0].inFileName + '" and keyfile "' + inTuple[1] + '".\n')
			except:
				sys.stdout.flush()
				if sys.exc_type == KeyboardInterrupt:
					pyRootPwa.utils.printInfo('Process ' + str(self.pid) + ' caught keyboard interrupt. Terminating...')
				else:
					pyRootPwa.utils.printErr('Process ' + str(self.pid) + ' caught exception. Terminating...')
					traceback.print_exc()
				break
			self.queue.task_done()

	def calcAmplitudes(self, inTuple):

		inFile = inTuple[0]
		keyfile = inTuple[1]
		outFile = inTuple[2]

		if pyRootPwa.config is None:
			raise pyRootPwa.rootPwaException("pyRootPwa configuration not initialized")

		outputFileFormat = pyRootPwa.config.outputFileFormat
		amplitudeLeafName = pyRootPwa.config.amplitudeLeafName
		outputCacheSize = pyRootPwa.config.outputCacheSize

		waveDescription = pyRootPwa.core.waveDescription()
		waveDescription.parseKeyFile(keyfile)
		(result, amplitude) = waveDescription.constructAmplitude()
		if not result:
			pyRootPwa.utils.printWarn('Could not construct amplitude for keyfile "' + keyfile + '".')
			return -1
		amplitude.init()
		print(amplitude)

		writeRootFile = False
		if outputFileFormat == "root":
			writeRootFile = True

		if writeRootFile:
			outFile.cd()
			ampTreeName = keyfile.rsplit('/',1)[-1].replace('.key', '.amp')
			outTree = pyRootPwa.ROOT.TTree(ampTreeName, ampTreeName)
			amplitudeTreeLeaf = pyRootPwa.core.amplitudeTreeLeaf()
			amplitudeTreeLeaf.branch(outTree, amplitudeLeafName)

		prodKinMomenta = pyRootPwa.ROOT.TClonesArray("TVector3")
		decayKinMomenta = pyRootPwa.ROOT.TClonesArray("TVector3")

		nEntries = len(inFile)
		if self.maxNmbEvents > 0:
			nEntries = min(nEntries, self.maxNmbEvents)
		inTree = inFile.tree
		prodKinMomenta = pyRootPwa.ROOT.TClonesArray("TVector3")
		decayKinMomenta = pyRootPwa.ROOT.TClonesArray("TVector3")
		inTree.SetBranchAddress(pyRootPwa.config.prodKinMomentaLeafName, prodKinMomenta)
		inTree.SetBranchAddress(pyRootPwa.config.decayKinMomentaLeafName, decayKinMomenta)

		decayTopology = amplitude.decayTopology()
		if not decayTopology.initKinematicsData(inFile.prodKinParticles, inFile.decayKinParticles):
			pyRootPwa.utils.printErr('Could not initialize kinematics Data "' + keyfile + '".')
			return -1

		if self.progressBar:
			progressbar = pyRootPwa.utils.progressBar(0, nEntries-1, sys.stdout)
			progressbar.start()
		vals = []
		if outputFileFormat == "binary":
			arrayAmp = array.array('d')

		try:
			for treeIndex in range(nEntries):
				inTree.GetEntry(treeIndex)
				if not decayTopology.readKinematicsData(prodKinMomenta, decayKinMomenta):
					if self.progressBar:
						progressbar.cancel()
					pyRootPwa.utils.printErr('Could not read kinematics data.')
					return -1
				amp = amplitude()
				vals.append(amp)
				if len(vals) > outputCacheSize:
					if outputFileFormat == "binary":
						arrayAmp.extend([number for complexNumber in vals for number in [complexNumber.real, complexNumber.imag]])
						arrayAmp.tofile(outFile)
						arrayAmp = array.array('d')
					elif outputFileFormat == "ascii":
						for val in vals:
							outFile.write("(" + str(val.real) + "," + str(val.imag) + ")\n")
					elif writeRootFile:
						for val in vals:
							amplitudeTreeLeaf.setAmp(val)
							outTree.Fill()
					else:
						raise Exception('Something is wrong, this should have been checked in the initialization of the configuration!')
					vals = []
				if self.progressBar:
					progressbar.update(treeIndex)
		except:
			if self.progressBar:
				progressbar.cancel()
			raise

		if outputFileFormat == "binary":
			arrayAmp.extend([number for complexNumber in vals for number in [complexNumber.real, complexNumber.imag]])
			arrayAmp.tofile(outFile)
		elif outputFileFormat == "ascii":
			for val in vals:
				outFile.write("(" + str(val.real) + "," + str(val.imag) + ")\n")
		elif writeRootFile:
			for val in vals:
				amplitudeTreeLeaf.setAmp(val)
				outTree.Fill()
			outTree.Write()
			outFile.Close()
		else:
			raise Exception('Something is wrong, this should have been checked in the initialization of the configuration!')

		return nEntries
