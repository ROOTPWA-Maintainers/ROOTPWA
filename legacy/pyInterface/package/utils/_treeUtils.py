
import sys

import pyRootPwa
import pyRootPwa.utils

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
	lineCounter = 0
	nLines = 0

	def __init__(self, infile):
		self.evtfile = infile
		beginning = self.evtfile.tell()
		nLinesPerParticle = int(self.evtfile.readline()[:-1]) + 1
		i = 1
		while self.evtfile.readline() != "":
			i += 1
		self.nLines = int(i / nLinesPerParticle)
		self.evtfile.seek(beginning)

	def __iter__(self):
		return self

	def __len__(self):
		return self.nLines

	def next(self):
		n_lines_to_read = self.evtfile.readline()[:-1]
		if n_lines_to_read == "":
			raise StopIteration()
		lines = [n_lines_to_read]
		for i in range(0, int(n_lines_to_read)):
			lines.append(self.evtfile.readline()[:-1])
			if lines[-1] == "":
				pyRootPwa.utils.printErr("Unexpected end of event file. Aborting...")
				sys.exit(1)
		return _Event(lines)

	def writeEvent(self, event):
		for line in event.lines:
			self.evtfile.write(line + "\n")


class _Event:

	lines = []

	particleNames = []
	physicsEvent = []

	def __init__(self, lines):
		self.lines = lines

	def sort(self):
		new_lines = self.lines[2:]
		new_lines = sorted(new_lines, key=lambda entry: int(entry.split()[0]))
		self.lines = self.lines[0:2] + new_lines
		self.physicsEvent = []
		self.particleNames = []

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

	def getPhysicsEvent(self):
		if self.physicsEvent:
			return self.physicsEvent
		nmbParticles = int(self.lines[0])
		part = self.__convertLineToPartProps(self.lines[1])
		self.physicsEvent.append(pyRootPwa.ROOT.TVector3(part[1], part[2], part[3]))
		fillPN = False
		if self.particleNames == []:
			fillPN = True
			self.particleNames.append(part[0])
		for i in range(2, nmbParticles + 1):
			part = self.__convertLineToPartProps(self.lines[i])
			if fillPN:
				self.particleNames.append(part[0])
			self.physicsEvent.append(pyRootPwa.ROOT.TVector3(part[1], part[2], part[3]))
		return self.physicsEvent

	def getParticleNames(self):
		if not self.particleNames:
			self.getPhysicsEvent()
		return self.particleNames

	def __str__(self):
		retval = ""
		for line in self.lines:
			retval += line + '\n'
		return retval[:-1]


def getTreeFromEvtFile(filename, treename = ""):

	if pyRootPwa.config is None:
		raise pyRootPwa.rootPwaException("pyRootPwa configuration not initialized")

	if treename == "":
		treename = str(hash(filename))

	outTree = pyRootPwa.ROOT.TTree(treename, treename)

	prodKinMomenta = pyRootPwa.ROOT.TClonesArray("TVector3")
	decayKinMomenta = pyRootPwa.ROOT.TClonesArray("TVector3")

	prodKinPartName = pyRootPwa.ROOT.TClonesArray("TObjString")
	decayKinPartName = pyRootPwa.ROOT.TClonesArray("TObjString")

	prodKinMomentaLeafName = pyRootPwa.config.prodKinMomentaLeafName
	decayKinMomentaLeafName= pyRootPwa.config.decayKinMomentaLeafName

	outTree.Branch(prodKinMomentaLeafName, "TClonesArray", prodKinMomenta)
	outTree.Branch(decayKinMomentaLeafName, "TClonesArray", decayKinMomenta)

	pyRootPwa.utils.printInfo('Converting "' + filename + '" to memory residing TTree...')

	with open(filename, 'r') as infile:

		inEventFile = _EventFile(infile)
		index = 0
		events = len(inEventFile)
		progressbar = pyRootPwa.utils.progressBar(0, events)
		progressbar.start()
		try:
			first = True
			for event in inEventFile:

				event.sort()
				physicsVectors = event.getPhysicsEvent()

				# check for the correct ordering of the names
				particleNames = event.getParticleNames()
				if first:
					prodKinPartName[0] = pyRootPwa.ROOT.TObjString(particleNames[0])
					for i in range(1, len(particleNames)):
						decayKinPartName[i-1] = pyRootPwa.ROOT.TObjString(particleNames[i])
				else:
					if len(particleNames) != prodKinPartName.GetEntriesFast() + decayKinPartName.GetEntriesFast():
						progressbar.cancel()
						raise pyRootPwa.rootPwaException("Mismatch between number of particle names in TClonesArray and number of particles in event")
					if prodKinPartName[0].GetString() != particleNames[0]:
						progressbar.cancel()
						raise pyRootPwa.rootPwaException("Inconsistent production particle types")
					for i in range(1, len(particleNames)):
						if decayKinPartName[i-1].GetString() != particleNames[i]:
							progressbar.cancel()
							raise pyRootPwa.rootPwaException("Inconsistent decay particle types")

				# set the physics vectors in the tree
				prodKinMomenta[0] = physicsVectors[0]
				for i in range(len(physicsVectors[1:])):
					decayKinMomenta[i] = physicsVectors[i+1]
				outTree.Fill()
				index += 1
				progressbar.update(index)
		except:
			progressbar.cancel()
			raise

	pyRootPwa.utils.printSucc('Successfully created TTree with ' + str(events) + ' events.')
	return (prodKinPartName, decayKinPartName, outTree)
