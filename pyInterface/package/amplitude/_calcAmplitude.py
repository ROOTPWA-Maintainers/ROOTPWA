
import pyRootPwa
import pyRootPwa.utils

def calcAmplitudes(inFile, keyfile, outfile):

	prodKinParticles = inFile.Get(pyRootPwa.config.get('amplitudes', 'prodKinPartNamesObjName'))
	decayKinParticles = inFile.Get(pyRootPwa.config.get('amplitudes', 'decayKinPartNamesObjName'))

	inTree = inFile.Get(pyRootPwa.config.get('amplitudes', 'inTreeName'))

	prodKinMomenta = pyRootPwa.ROOT.TClonesArray("TVector3")
	decayKinMomenta = pyRootPwa.ROOT.TClonesArray("TVector3")

	inTree.SetBranchAddress(pyRootPwa.config.get('amplitudes', 'prodKinMomentaLeafName'), prodKinMomenta)
	inTree.SetBranchAddress(pyRootPwa.config.get('amplitudes', 'decayKinMomentaLeafName'), decayKinMomenta)

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
		outfile.write("(" + str(amp.real) + ", " + str(amp.imag) + ")\n")
	return True

