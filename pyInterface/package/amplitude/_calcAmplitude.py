
import pyRootPwa

def calcAmplitudes(inFile, keyfile, outfile):

	prodKinParticles = inFile.Get(pyRootPwa.config.get('amplitudes', 'prodKinPartNamesObjName'))
	decayKinParticles = inFile.Get(pyRootPwa.config.get('amplitudes', 'decayKinPartNamesObjName'))

	inTree = inFile.Get(pyRootPwa.config.get('amplitudes', 'inTreeName'))

	prodKinMomenta = pyRootPwa.ROOT.TVector3()
	decayKinMomenta = pyRootPwa.ROOT.TVector3()

	inTree.SetBranchAddress(pyRootPwa.config.get('amplitudes', 'prodKinMomentaLeafName'), prodKinMomenta)
	inTree.SetBranchAddress(pyRootPwa.config.get('amplitudes', 'decayKinMomentaLeafName'), decayKinMomenta)

	waveDesc = pyRootPwa.waveDescription()
	waveDesc.parseKeyFile(keyfile)
	(waveDescConstructionSuccess, amplitude) = waveDesc.constructAmplitude()

	amplitude.decayTopology().initKinematicsData(prodKinParticles, decayKinParticles)

	for treeIndex in range(inTree.GetEntries()):
		inTree.getEntry(treeIndex)
		amplitude.decayTopology().readKinematicsData(prodKinMomenta, decayKinMomenta)
		print(amplitude())

