
import array
import sys

import pyRootPwa
import pyRootPwa.exception
import pyRootPwa.utils

def calcAmplitudes(inFile, keyfile, outFile):

	if pyRootPwa.config is None:
		raise pyRootPwa.exception.pyRootPwaException("pyRootPwa configuration not initialized")

	prodKinMomenta = pyRootPwa.ROOT.TClonesArray("TVector3")
	decayKinMomenta = pyRootPwa.ROOT.TClonesArray("TVector3")

	inFile.tree.SetBranchAddress(pyRootPwa.config.prodKinMomentaLeafName, prodKinMomenta)
	inFile.tree.SetBranchAddress(pyRootPwa.config.decayKinMomentaLeafName, decayKinMomenta)

	pythonAdmin = pyRootPwa.pythonAdministrator()

	writeRootFile = False
	if pyRootPwa.config.outputFileFormat == "root":
		writeRootFile = True

	if writeRootFile:
		outFile.cd()
		ampTreeName = keyfile.rsplit('/',1)[-1].replace('.key', '.amp')
		outTree = pyRootPwa.ROOT.TTree(ampTreeName, ampTreeName)
		amplitudeTreeLeaf = pyRootPwa.amplitudeTreeLeaf()
		pythonAdmin.branch(outTree, amplitudeTreeLeaf, pyRootPwa.config.amplitudeLeafName)

	if not pythonAdmin.constructAmplitude(keyfile):
		pyRootPwa.utils.printWarn('Could not construct amplitude for keyfile "' + keyfile + '".')
		return False
	sys.stdout.write(str(pythonAdmin))
	if not pythonAdmin.initKinematicsData(inFile.prodKinParticles, inFile.decayKinParticles):
		pyRootPwa.utils.printErr('Could not initialize kinematics Data "' + keyfile + '".')
		return False

	progressbar = pyRootPwa.utils.progressBar(0, inFile.tree.GetEntries())
	progressbar.start()
	for treeIndex in range(inFile.tree.GetEntries()):
		inFile.tree.GetEntry(treeIndex)
		if not pythonAdmin.readKinematicsData(prodKinMomenta, decayKinMomenta):
			progressbar.cancel()
			pyRootPwa.utils.printErr('Could not read kinematics data.')
			return False
		amp = pythonAdmin()
		if pyRootPwa.config.outputFileFormat == "ascii":
			outFile.write("(" + str(amp.real) + "," + str(amp.imag) + ")\n")
		elif pyRootPwa.config.outputFileFormat == "binary":
			arrayAmp = array.array('d', [amp.real, amp.imag])
			arrayAmp.tofile(outFile)
		elif writeRootFile:
			amplitudeTreeLeaf.setAmp(amp)
			outTree.Fill()
		else:
			raise Exception('Something is wrong, this should have been checked in the initialization of the configuration!')
		progressbar.update(treeIndex)

	if writeRootFile:
		outTree.Write()
		outFile.Close()

	return True

