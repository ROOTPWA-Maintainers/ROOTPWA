
import array

import pyRootPwa
import pyRootPwa.utils

def calcAmplitudes(inFile, keyfile, outfile):

	prodKinParticles = inFile.Get(pyRootPwa.config.prodKinPartNamesObjName)
	decayKinParticles = inFile.Get(pyRootPwa.config.decayKinPartNamesObjName)

	inTree = inFile.Get(pyRootPwa.config.inTreeName)

	prodKinMomenta = pyRootPwa.ROOT.TClonesArray("TVector3")
	decayKinMomenta = pyRootPwa.ROOT.TClonesArray("TVector3")

	inTree.SetBranchAddress(pyRootPwa.config.prodKinMomentaLeafName, prodKinMomenta)
	inTree.SetBranchAddress(pyRootPwa.config.decayKinMomentaLeafName, decayKinMomenta)

	pythonAdmin = pyRootPwa.pythonAdministrator()
	if not pythonAdmin.constructAmplitude(keyfile):
		pyRootPwa.utils.printWarn('Could not construct amplitude for keyfile "' + keyfile + '".')
		return False
	if not pythonAdmin.initKinematicsData(prodKinParticles, decayKinParticles):
		pyRootPwa.utils.printErr('Could not initialize kinematics Data "' + keyfile + '".')
		return False

	for treeIndex in range(inTree.GetEntries()):
		inTree.GetEntry(treeIndex)
		if not pythonAdmin.readKinematicsData(prodKinMomenta, decayKinMomenta):
			pyRootPwa.utils.printErr('Could not read kinematics data.')
			return False
		amp = pythonAdmin()
		if pyRootPwa.config.outputFileFormat == "ascii":
			outfile.write("(" + str(amp.real) + "," + str(amp.imag) + ")\n")
		elif pyRootPwa.config.outputFileFormat == "binary":
			arrayAmp = array.array('d', [amp.real, amp.imag])
			arrayAmp.tofile(outfile)
		else:
			raise Exception('Something is wrong, this should have been checked in the initialization of the configuration!')
	return True

