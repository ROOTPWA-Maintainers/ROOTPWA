
import array
import sys

import pyRootPwa
import pyRootPwa.utils

def calcAmplitudes(inFileName, keyfile, outFile):

	if inFileName.endswith('.root'):
		inFile = pyRootPwa.ROOT.TFile.Open(inFileName)
		prodKinParticles = inFile.Get(pyRootPwa.config.prodKinPartNamesObjName)
		decayKinParticles = inFile.Get(pyRootPwa.config.decayKinPartNamesObjName)
		inTree = inFile.Get(pyRootPwa.config.inTreeName)
	else:
		(prodKinParticles, decayKinParticles, inTree) = pyRootPwa.utils.getTreeFromEvtFile(inFileName, inFileName)

	prodKinMomenta = pyRootPwa.ROOT.TClonesArray("TVector3")
	decayKinMomenta = pyRootPwa.ROOT.TClonesArray("TVector3")

	inTree.SetBranchAddress(pyRootPwa.config.prodKinMomentaLeafName, prodKinMomenta)
	inTree.SetBranchAddress(pyRootPwa.config.decayKinMomentaLeafName, decayKinMomenta)

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
	if not pythonAdmin.initKinematicsData(prodKinParticles, decayKinParticles):
		pyRootPwa.utils.printErr('Could not initialize kinematics Data "' + keyfile + '".')
		return False

	progressbar = pyRootPwa.utils.progressBar(0, inTree.GetEntries())
	progressbar.start()
	for treeIndex in range(inTree.GetEntries()):
		inTree.GetEntry(treeIndex)
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

