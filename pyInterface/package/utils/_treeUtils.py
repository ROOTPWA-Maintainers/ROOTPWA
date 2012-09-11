
import sys

import pyRootPwa

_geantParticleCodes = {}
_geantParticleCodes[0] = "unknown"
_geantParticleCodes[1] = "gamma"
_geantParticleCodes[2] = "e"
_geantParticleCodes[3] = "e"
_geantParticleCodes[7] = "pi0"
_geantParticleCodes[8] = "pi"
_geantParticleCodes[9] = "pi"
_geantParticleCodes[11] = "K"
_geantParticleCodes[12] = "K"
_geantParticleCodes[13] = "n"
_geantParticleCodes[14] = "p"
_geantParticleCodes[15] = "pbar"
_geantParticleCodes[16] = "K0"
_geantParticleCodes[17] = "eta"
_geantParticleCodes[18] = "lambda"
_geantParticleCodes[57] = "rho(770)"
_geantParticleCodes[58] = "rho(770)"
_geantParticleCodes[59] = "rho(770)"
_geantParticleCodes[60] = "omega(782)"
_geantParticleCodes[61] = "eta'(958)"
_geantParticleCodes[62] = "phi(1020)"
_geantParticleCodes[45] = "d"


class _EventFile:

	evtfile = None

	def __init__(self, infile):
		self.evtfile = infile


	def get_event(self):
		n_lines_to_read = self.evtfile.readline()[:-1]
		if n_lines_to_read == "":
			return None
		lines = [n_lines_to_read]
		for i in range(0, int(n_lines_to_read)):
			lines.append(self.evtfile.readline()[:-1])
			if lines[-1] == "":
				pyRootPwa.utils.printErr("Unexpected end of event file. Aborting...")
				sys.exit(1)
		return _Event(lines)


	def write_event(self, event):
		for line in event.lines:
			self.evtfile.write(line + "\n")



class _Event:

	lines = []

	particleNames = []

	def __init__(self, lines):
		self.lines = lines


	def sort(self):
		new_lines = self.lines[2:]
		new_lines = sorted(new_lines, key=lambda entry: int(entry.split()[0]))
		self.lines = self.lines[0:2] + new_lines

	def __convertLineToPartProps(self, line):
		part = line.split(' ')
		(part[0], part[1]) = (int(part[0]), int(part[1]))
		for j in range(2, 6):
			part[j] = float(part[j])
		partname = _geantParticleCodes[part[0]]
		if part[1] > 0:
			partname += "+"
		elif part[1] < 0:
			partname += "-"
		else:
			partname += "0"
		part[0] = partname
		part.pop(1)
		part.pop(len(part)-1)
		return part

	def get_physics_event(self):
		retval = []
		nmbParticles = int(self.lines[0])
		part = self.__convertLineToPartProps(self.lines[1])
		retval.append(pyRootPwa.ROOT.TVector3(part[1], part[2], part[3]))
		fillPN = False
		if self.particleNames == []:
			fillPN = True
			self.particleNames.append(part[0])
		for i in range(2, nmbParticles + 1):
			part = self.__convertLineToPartProps(self.lines[i])
			if fillPN:
				self.particleNames.append(part[0])
			retval.append(pyRootPwa.ROOT.TVector3(part[1], part[2], part[3]))
		return retval

	def get_particle_names(self):
		if self.particleNames == []:
			self.get_physics_event()
		return self.particleNames

	def __str__(self):
		retval = ""
		for line in self.lines:
			retval += line + '\n'
		return retval[:-1]


def getTreeFromEvtFile(filename, treename = ""):

	if pyRootPwa.config is None:
		pyRootPwa.utils.printErr("Config singleton not initialized. Aborting...")
		sys.exit(1)

	if treename == "":
		treename = str(hash(filename))

	outTree = pyRootPwa.ROOT.TTree(treename, treename)

	prodKinMomenta = pyRootPwa.ROOT.TClonesArray("TVector3")
	decayKinMomenta = pyRootPwa.ROOT.TClonesArray("TVector3")

	prodKinPartName = pyRootPwa.ROOT.TClonesArray("TObjString")
	decayKinPartName = pyRootPwa.ROOT.TClonesArray("TObjString")

	prodKinMomentaLeafName = pyRootPwa.config.get('amplitudes', 'prodKinMomentaLeafName')
	decayKinMomentaLeafName= pyRootPwa.config.get('amplitudes', 'decayKinMomentaLeafName')

	outTree.Branch(prodKinMomentaLeafName, "TClonesArray", prodKinMomenta)
	outTree.Branch(decayKinMomentaLeafName, "TClonesArray", decayKinMomenta)

	first = True
	with open(filename, 'r') as infile:

		inEventFile = _EventFile(infile)
		event = inEventFile.get_event()
		while event is not None:

			event.sort()
			physicsVectors = event.get_physics_event()

			# check for the correct ordering of the names
			particleNames = event.get_particle_names()
			assert(len(physicsVectors) == len(particleNames))
			if first:
				prodKinPartName[0] = pyRootPwa.ROOT.TObjString(particleNames[0])
				for i in range(1, len(particleNames)):
					decayKinPartName[i-1] = pyRootPwa.ROOT.TObjString(particleNames[i])
			else:
				assert(len(particleNames) == prodKinPartName.GetEntriesFast() + decayKinPartName.GetEntriesFast())
				assert(prodKinPartName[0].GetString() == particleNames[0])
				for i in range(1, len(particleNames)):
					assert(decayKinPartName[i-1].GetString() == particleNames[i])

			# set the physics vectors in the tree
			prodKinMomenta[0] = physicsVectors[0]
			for i in range(len(physicsVectors[1:])):
				decayKinMomenta[i] = physicsVectors[i+1]
			outTree.Fill()
			event = inEventFile.get_event()

	return (prodKinPartName, decayKinPartName, outTree)


