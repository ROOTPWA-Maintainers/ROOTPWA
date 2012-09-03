
import sys

success = True
skip = False

# Some functions
# ---------------------------------------------------------

def print_yellow(string):
	print('\033[93m' + string + '\033[0m')

def print_green(string):
	print('\033[92m' + string + '\033[0m')

def print_red(string):
	print('\033[91m' + string + '\033[0m')

def do_test(function, name, skip_test = False):
	sys.stdout.write(name + "...")
	if skip_test:
		global skip
		skip = True
		print_yellow("skipped")
		return None
	try:
		retval = (function)()
	except:
		print_red("error")
		global success
		success = False
		return None
	print_green("success")
	return retval

def do_test_raw(function, name):
	print(name)
	retval = (function)()

# ---------------------------------------------------------

# ---------------------------------------------------------
#
#	General Stuff
#
# ---------------------------------------------------------

def impLib(): import pyRootPwa
do_test(impLib, "Importing pyRootPwa")

import pyRootPwa

print
print("########################################################################")
print

# ---------------------------------------------------------
#
#	particleDataTable
#
# ---------------------------------------------------------

def getpDTinstance(): return pyRootPwa.particleDataTable.instance
particleTable = do_test( getpDTinstance, "Testing particleDataTable.instance")

def tPDTreadFile():
	print("\n")
	pyRootPwa.particleDataTable.readFile("../../amplitude/particleDataTable.txt")
	print
do_test(tPDTreadFile, "Testing particleDataTable.readFile()")

print
print("########################################################################")
print

# ---------------------------------------------------------
#
#	particlePropertiesTable
#
# ---------------------------------------------------------

def partPropTestConst(): return pyRootPwa.particleProperties()
partProp = do_test(partPropTestConst, "Testing particleProperties constructor")

def partPropTestCopyConst(): return pyRootPwa.particleProperties(partProp)
partProp2 = do_test(partPropTestCopyConst, "Testing particleProperties copy constructor")

def partPropTestOps():
	assert(partProp == partProp2)
	old_name = partProp2.name
	partProp2.name = "bliblablup"
	assert(partProp != partProp2)
	partProp2.name = old_name
do_test(partPropTestOps, "Testing particleProperties \"==\"/\"!=\" operators...")

print
print("########################################################################")
print

# ---------------------------------------------------------
#
#	particle	
#
# ---------------------------------------------------------

def particleConst(): return pyRootPwa.particle()
part = do_test(particleConst, "Testing particle default constructor")

def particleRead(): part.read("Delta(1910)+     Deltabar(1910)-  1.91        0.25        +1           3    0    0    0     0     1    +1     0")
do_test(particleRead, "Testing particle.read()")

def printPart():
	print("\n\n" + str(part) + "\n")
do_test(printPart, "Testing print(particle)")

def partCopyConst():
	p2 = pyRootPwa.particle(part)
	assert(part == p2)
do_test(partCopyConst, "Testing particle copy constructor")

def partTestConsts():
	pP = pyRootPwa.particleProperties()
	pP.read("Delta(1910)+     Deltabar(1910)-  1.91        0.25        +1           3    0    0    0     0     1    +1     0")
	t = pyRootPwa.ROOT.TVector3(1., 1., 1.)
	
	p3 = pyRootPwa.particle(pP)
	p3 = pyRootPwa.particle(pP, -1)
	p3 = pyRootPwa.particle(pP, -1, 0)
	p3 = pyRootPwa.particle(pP, -1, 0, 0)
	p3 = pyRootPwa.particle(pP, -1, 0, 0, t)
	
	p3 = pyRootPwa.particle("pi+")
	p3 = pyRootPwa.particle("pi+", True)
	p3 = pyRootPwa.particle("pi+", True, -1)
	p3 = pyRootPwa.particle("pi+", True, -1, 0)
	p3 = pyRootPwa.particle("pi+", True, -1, 0, 0)
	p3 = pyRootPwa.particle("pi+", True, -1, 0, 0, t)
	
	p3 = pyRootPwa.particle("pi-", 0, 1, 2, 3, 4, 5)
	p3.qnSummary()
	p3 = pyRootPwa.particle("pi-", 0, 1, 2, 3, 4, 5, 6)
	p3 = pyRootPwa.particle("pi-", 0, 1, 2, 3, 4, 5, 6, 7)
do_test(partTestConsts, "Testing particle other constructors")

def partTestClone(): p2 = part.clone()
do_test(partTestClone, "Testing particle.clone()")

def partTestqnSummary(): part.qnSummary()
do_test(partTestqnSummary, "Testing particle.qnSummary()")

def partTestLabel(): assert(part.label() == "Delta(1910)+[1.5(0.5+)]")
do_test(partTestLabel, "Testing particle.label()")

def partTestMomentum():
	part.momentum = pyRootPwa.ROOT.TVector3(10, 10, 10)
	assert(part.momentum.X() == 10.)
do_test(partTestMomentum, "Testing particle.momentum")

def partTestlzVec():
	lz = pyRootPwa.ROOT.TLorentzVector(1., 1., 1., 1.)
	part.lzVec = lz
	assert(part.lzVec.E() == 1.)
do_test(partTestlzVec, "Testing particle.lzVec")

def partTestIndex():
	part.index = 12
	assert(part.index == 12)
do_test(partTestIndex, "Testing particle.index")

def partTestReflectivity():
	part.reflectivity = 2
	assert(part.reflectivity == 1)
do_test(partTestReflectivity, "Testing part.reflectivity")

def partTestSetProperties():
	pP = pyRootPwa.particleProperties()
	pP.read("Kstar2(1430)+    Kstar2bar(1430)- 1.4256      0.0985      0            1   +1    0    0     0     4    +1     0")
	part.setProperties(pP)
	assert(part.name == "Kstar2(1430)+")
do_test(partTestSetProperties, "Testing particle.setProperties()")

def partTestTransformTV3():
	tv3 = pyRootPwa.ROOT.TVector3(1., 1., 1.)
	dump = part.transform(tv3)
	assert(dump is not None)
do_test(partTestTransformTV3, "Testing particle.transform(TVector3)")

def partTestTransformTLR():
	t = pyRootPwa.ROOT.TLorentzRotation()
	dump = part.transform(t)
	assert(dump is not None)
do_test(partTestTransformTLR, "Testing particle.transform(TLorentzRotation)")

def partTestTransformFAIL():
	print("\n")
	t = pyRootPwa.ROOT.TH1D()
	assert(part.transform(t) is None)
	t = [0, 1, 2]
	assert(part.transform(t) is None)
	print
do_test(partTestTransformFAIL, "Testing particle.transform(UNSUPPORTED TYPE)")

def partTestDebugFlag():
	old_pP_debug = pyRootPwa.particleProperties.debug
	old_debug = pyRootPwa.particle.debug
	assert(old_pP_debug == old_debug)
	pyRootPwa.particle.debug = (not old_debug)
	assert(pyRootPwa.particle.debug == (not old_debug))
	assert(pyRootPwa.particleProperties.debug == old_pP_debug)
	pyRootPwa.particle.debug = old_debug
	assert(part.debug == old_debug)
	part.debug = (not old_debug)
	assert(part.debug == (not old_debug))
	assert(pyRootPwa.particleProperties.debug == old_pP_debug)
	part.debug = old_debug
do_test(partTestDebugFlag, "Testing particle debug flag", True)
# This test fails for unknown reason

print
print("########################################################################")
print

# ---------------------------------------------------------
#
#	interactionVertex
#
# ---------------------------------------------------------

def defConst(): return pyRootPwa.interactionVertex()
iV = do_test(defConst, "Testing interactionVertex default constructor")

def copyConst(): return pyRootPwa.interactionVertex(iV)
do_test(copyConst, "Testing interactionVertex copy constructor")

def tClone():
	iV2 = iV.clone()
	iV2 = iV.clone(True)
	iV2 = iV.clone(True, True)
do_test(tClone, "Testing interactionVertex::clone")

def tPrint(): print("\n\n" + str(iV) + "\n")
do_test(tPrint, "Testing \"print(interactionVertex)\"")

def tClear():
	class tiV(pyRootPwa.interactionVertex):
		def clear(self):
			return "testString"
	iV.clear()
	iV2 = tiV(iV)
	assert(iV2.clear() == "testString")
do_test(tClear, "Testing interactionVertex::clear()")

def tAiP():
	p = pyRootPwa.particle()
	assert(iV.addInParticle(p))
do_test(tAiP, "Testing interactionVertex::addInParticle()")

def tAoP():
	p = pyRootPwa.particle()
	assert(iV.addOutParticle(p))
do_test(tAoP, "Testing interactionVertex::addOutParticle()")

def ttOP():
	rot = pyRootPwa.ROOT.TLorentzRotation(1, 1, 1)
	iV.transformOutParticles(rot)
do_test(ttOP, "Testing interactionVertex::transformOutParticles()")

def tNmbIP(): assert(iV.nmbInParticles == 1)
do_test(tNmbIP, "Testing interactionVertex::nmbInParticles")

def tNmbOP(): assert(iV.nmbOutParticles == 1)
do_test(tNmbOP, "Testing interactionVertex::nmbOutParticles")

def tInParts(): assert(len(iV.inParticles()) == 1)
do_test(tInParts, "Testing interactionVertex::inParticles()")

def tOutParts(): assert(len(iV.outParticles()) == 1)
do_test(tOutParts, "Testing interactionVertex::outParticles()")

def tName(): assert(iV.name() == "interactionVertex")
do_test(tName, "Testing interactionVertex::name()")

def tDebug():
	old_debug = iV.debug
	iV.debug = (not old_debug)
	assert(iV.debug == (not old_debug))
	iV.debug = old_debug
do_test(tDebug, "Testing \"debug\" property")

print
print("########################################################################")
print

# ---------------------------------------------------------
#
#	fsVertex
#
# ---------------------------------------------------------

def fsVertexTestConsts():
	fsVert = pyRootPwa.fsVertex(part)
	fsVert2 = pyRootPwa.fsVertex(fsVert)
	return fsVert
fsVert = do_test(fsVertexTestConsts, "Testing fsVertex constructors")

def fsVertexTestClone():
	fsVert2 = fsVert.clone()
	fsVert2 = fsVert.clone(True)
	fsVert2 = fsVert.clone(True, True)
do_test(fsVertexTestClone, "Testing fsVertex.clone()")

def fsVertexTestPrint(): print("\n\n" + str(fsVert) + "\n")
do_test(fsVertexTestPrint, "Testing print(fsVertex)")

def fsVertexTestAddInOutPart():
	assert(not fsVert.addInParticle(part))
	assert(not fsVert.addOutParticle(part))
do_test(fsVertexTestAddInOutPart, "Testing fsVertex.add{In/Out}Particle()")

def fsVertexTestfsParticle():
	part2 = fsVert.fsParticle()
	assert(part == part2)
do_test(fsVertexTestfsParticle, "Testing fsVertex.fsParticle")

def fsVertexTestName(): assert(fsVert.name() == "fsVertex")
do_test(fsVertexTestName, "Testing fsVertex.name()")

def fsVertexTestDebugFlag():
	old_iV_debug = pyRootPwa.interactionVertex.debug
	old_fsV_debug = pyRootPwa.fsVertex.debug
	pyRootPwa.fsVertex.debug = (not old_fsV_debug)
	assert(pyRootPwa.fsVertex.debug == (not old_fsV_debug))
	assert(pyRootPwa.interactionVertex.debug == old_iV_debug)
	pyRootPwa.fsVertex.debug = old_fsV_debug
do_test(fsVertexTestDebugFlag, "Testing particle debug flag", True)
# This test fails for unknown reason

print
print("########################################################################")
print

# ---------------------------------------------------------
#
#	massDependence
#
# ---------------------------------------------------------

def flatMassDepTestConst(): return pyRootPwa.flatMassDependence()
flatMassDep = do_test(flatMassDepTestConst, "Testing flatMassDependence default constructor")

def flatMassDepTestAmp():
	pass
do_test(flatMassDepTestAmp, "Testing flatMassDependence.amp()", True)
# isobarDecayVertex is missing at the time of this being written.

def flatMassDepTestName(): assert(flatMassDep.name() == "flatMassDependence")
do_test(flatMassDepTestName, "Testing flatMassDependence.name()")

def relBreitWigTestConst(): return pyRootPwa.relativisticBreitWigner()
relBreitWig = do_test(relBreitWigTestConst, "Testing relativisticBreitWigner default constructor")

def relBreitWigTestAmp():
	pass
do_test(relBreitWigTestAmp, "Testing relativisticBreitWigner.amp()", True)
# isobarDecayVertex is missing at the time of this being written.

def relBreitWigTestName(): assert(relBreitWig.name() == "relativisticBreitWigner")
do_test(relBreitWigTestName, "Testing relativisticBreitWigner.name()")

def SAuMoPenMTestConst(): return pyRootPwa.piPiSWaveAuMorganPenningtonM()
SAuMoPenM = do_test(SAuMoPenMTestConst, "Testing piPiSWaveAuMorganPenningtonM default constructor")

def SAuMoPenMTestAmp():
	pass
do_test(SAuMoPenMTestAmp, "Testing piPiSWaveAuMorganPenningtonM.amp()", True)
# isobarDecayVertex is missing at the time of this being written.

def SAuMoPenMTestName(): assert(SAuMoPenM.name() == "piPiSWaveAuMorganPenningtonM")
do_test(SAuMoPenMTestName, "Testing piPiSWaveAuMorganPenningtonM.name()")

def SAuMoPenVesTestConst(): return pyRootPwa.piPiSWaveAuMorganPenningtonVes()
SAuMoPenVes = do_test(SAuMoPenVesTestConst, "Testing piPiSWaveAuMorganPenningtonVes default constructor")

def SAuMoPenVesTestAmp():
	pass
do_test(SAuMoPenVesTestAmp, "Testing piPiSWaveAuMorganPenningtonVes.amp()", True)
# isobarDecayVertex is missing at the time of this being written.

def SAuMoPenVesTestName(): assert(SAuMoPenVes.name() == "piPiSWaveAuMorganPenningtonVes")
do_test(SAuMoPenVesTestName, "Testing piPiSWaveAuMorganPenningtonVes.name()")

def SAuMoPenKachaevTestConst(): return pyRootPwa.piPiSWaveAuMorganPenningtonKachaev()
SAuMoPenKachaev = do_test(SAuMoPenKachaevTestConst, "Testing piPiSWaveAuMorganPennigtonKachaev default constructor")

def SAuMoPenKachaevTestAmp():
	pass
do_test(SAuMoPenKachaevTestAmp, "Testing piPiSWaveAuMorganPennigtonKachaev.amp()", True)
# isobarDecayVertex is missing at the time of this being written.

def SAuMoPenKachaevTestName(): assert(SAuMoPenKachaev.name() == "piPiSWaveAuMorganPenningtonKachaev")
do_test(SAuMoPenKachaevTestName, "Testing piPiSWaveAuMorganPennigtonKachaev.name()")

print
print("########################################################################")
print

# ---------------------------------------------------------
#
#	isobarDecayVertex
#
# ---------------------------------------------------------

def isobDecVtxTestConstructor():
	lz = pyRootPwa.ROOT.TLorentzVector(1., 1., 1., 1.)
	part.lzVec = lz
	return pyRootPwa.isobarDecayVertex(part, part, part)
isobDecVtx = do_test(isobDecVtxTestConstructor, "Testing isobarDecayVertex constructor")

def isobarDecVtxTestPrint(): print("\n\n" + str(isobDecVtx) + "\n")
do_test(isobarDecVtxTestPrint, "Testing print(isobarDecayVertex)")

def isobarDecVtxTestClone():
	v2 = isobDecVtx.clone()
	v2 = isobDecVtx.clone(True)
	v2 = isobDecVtx.clone(True, True)
do_test(isobarDecVtxTestClone, "Testing isobarDecayVertex.clone()")

def isobarDecVtxTestAddInOutPart():
	assert(not isobDecVtx.addInParticle(part))
	assert(not isobDecVtx.addOutParticle(part))
do_test(isobarDecVtxTestAddInOutPart, "Testing isobarDecayVertex.add{In/Out}Particle()")

def isobarDecVtxTestParent(): assert(part == isobDecVtx.parent())
do_test(isobarDecVtxTestParent, "Testing isobarDecayVertex.parent()")

def isobarDecVtxTestDaughters():
	assert(part == isobDecVtx.daughter1())
	assert(part == isobDecVtx.daughter2())
do_test(isobarDecVtxTestDaughters, "Testing isobarDecayVertex.daughter{1/2}()")

def isobarDecVtxTestCalcLZ():
	lz = isobDecVtx.calcParentLzVec()
	lz2 = pyRootPwa.ROOT.TLorentzVector(2., 2., 2., 2.)
	assert(lz == lz2)
do_test(isobarDecVtxTestCalcLZ, "Testing isobarDecayVertex.calcParentLzVec()")

def isobarDecVtxTestParC(): assert(isobDecVtx.calcParentCharge() == 2)
do_test(isobarDecVtxTestParC, "Testing isobarDecayVertex.calcParentCharge()")

def isobarDecVtxTestBNC(): assert(isobDecVtx.calcParentBaryonNmb() == 0)
do_test(isobarDecVtxTestBNC, "Testing isobarDecayVertex.calcParentBaryonNmb()")

def isobarDecVtxTestL():
	old_L = isobDecVtx.L
	isobDecVtx.L = 12
	assert(isobDecVtx.L == 12)
	isobDecVtx.L = old_L
do_test(isobarDecVtxTestL, "Testing isobarDecayVertex.L")

def isobarDecVtxTestS():
	old_S = isobDecVtx.S
	isobDecVtx.S = 12
	assert(isobDecVtx.S == 12)
	isobDecVtx.S = old_S
do_test(isobarDecVtxTestS, "Testing isobarDecayVertex.S")

def isobarDecVtxTestMDA(): assert(isobDecVtx.massDepAmplitude() == (1+0j))
do_test(isobarDecVtxTestMDA, "Testing isobarDecayVertex.massDepAmplitude()")

def isobarDecVtxTestMD():
	mDname = isobDecVtx.massDependence().name() 
do_test(isobarDecVtxTestMD, "Testing isobarDecayVertex.massDependence()")

print
print("########################################################################")
print

# ---------------------------------------------------------
#
#	Summary
#
# ---------------------------------------------------------

print
if success:
	print_green("All tests successful.")
else:
	print_red("There were errors.")
if skip:
	print_yellow("Some tests were skipped.")
