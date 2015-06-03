
import math
import os
import sys

DIGITS_TO_ROUND_TO = 14
DIGITS_TO_ROUND_GJTRAFO_TO = 11

errors = 0
skip = 0

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
		skip += 1
		print_yellow("skipped")
		return None
	try:
		retval = (function)()
	except:
		print_red("error")
		global errors
		errors += 1
		return None
	print_green("success")
	return retval

def do_test_raw(function, name):
	sys.stdout.write(name + "...")
	retval = (function)()
	print_green("success")
	return retval

# ---------------------------------------------------------

# ---------------------------------------------------------
#
#	General Stuff
#
# ---------------------------------------------------------

def impLib(name):
	sys.stdout.write(name + "...")
	try:
		import pyRootPwa
		import pyRootPwa.core
	except:
		print_red("error")
		return False
	print_green("success")
	return True
if not impLib("Importing pyRootPwa"):
	print("Could not import library, aborting tests...")
	sys.exit(1)

import pyRootPwa
import pyRootPwa.core

def testPrintingUtils():
	print("\n")
	pyRootPwa.core.printCompilerInfo()
	pyRootPwa.core.printLibraryInfo()
	pyRootPwa.core.printGitHash()
	print
do_test(testPrintingUtils, "Testing printing utilities")

print
print("########################################################################")
print

# ---------------------------------------------------------
#
#	particleDataTable
#
# ---------------------------------------------------------

def getpDTinstance(): return pyRootPwa.core.particleDataTable.instance
particleTable = do_test( getpDTinstance, "Testing particleDataTable.instance")

def tPDTreadFile():
	print("\n")
	pyRootPwa.core.particleDataTable.readFile(os.environ['ROOTPWA'] + "/particleData/particleDataTable.txt")
	print
do_test(tPDTreadFile, "Testing particleDataTable.readFile()")

def tPDTiterator():
	for part in particleTable:
		assert(part.first == part.second.name)
do_test(tPDTiterator, "Testing particleDataTable iterator")

def tPDTentriesMatching():
	pP = particleTable.entry("rho(770)0")
	parts = pyRootPwa.core.particleDataTable.entriesMatching(pP, "allQn", 0, 0, [], [], [], False)
	assert(len(parts) == 5)
	for part in parts:
		assert(part.name[:3] == "rho")
do_test(tPDTentriesMatching, "Testing particleDataTable.entriesMatching()")

def tPDDebug():
	old_debug = particleTable.debugParticleDataTable
	particleTable.debugParticleDataTable = (not old_debug)
	assert(particleTable.debugParticleDataTable == (not old_debug))
	particleTable.debugParticleDataTable = old_debug
do_test(tPDDebug, "Testing particleDataTable debug flag")

print
print("########################################################################")
print

# ---------------------------------------------------------
#
#	particleProperties
#
# ---------------------------------------------------------

def partPropTestConst(): return pyRootPwa.core.particleProperties()
partProp = do_test(partPropTestConst, "Testing particleProperties constructor")

def partPropTestCopyConst(): return pyRootPwa.core.particleProperties(partProp)
partProp2 = do_test(partPropTestCopyConst, "Testing particleProperties copy constructor")

def partPropChargeFromName():
	tup = pyRootPwa.core.particleProperties.chargeFromName("bla+")
	assert(tup[0] == 'bla' and tup[1] == 1)
do_test(partPropChargeFromName, "Testing particleProperties::chargeFromName")

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

def particleConst(): return pyRootPwa.core.particle()
part = do_test(particleConst, "Testing particle default constructor")

def particleRead(): part.read("Delta(1910)+     Deltabar(1910)-  1.91        0.25        +1           3    0    0    0     0     1    +1     0")
do_test(particleRead, "Testing particle.read()")

def printPart():
	print("\n\n" + str(part) + "\n")
do_test(printPart, "Testing print(particle)")

def partCopyConst():
	p2 = pyRootPwa.core.particle(part)
	assert(part == p2)
do_test(partCopyConst, "Testing particle copy constructor")

def partTestConsts():
	pP = pyRootPwa.core.particleProperties()
	pP.read("Delta(1910)+     Deltabar(1910)-  1.91        0.25        +1           3    0    0    0     0     1    +1     0")
	t = pyRootPwa.ROOT.TVector3(1., 1., 1.)

	p3 = pyRootPwa.core.particle(pP)
	p3 = pyRootPwa.core.particle(pP, -1)
	p3 = pyRootPwa.core.particle(pP, -1, 0)
	p3 = pyRootPwa.core.particle(pP, -1, 0, 0)
	p3 = pyRootPwa.core.particle(pP, -1, 0, 0, t)

	p3 = pyRootPwa.core.particle("pi+")
	p3 = pyRootPwa.core.particle("pi+", True)
	p3 = pyRootPwa.core.particle("pi+", True, -1)
	p3 = pyRootPwa.core.particle("pi+", True, -1, 0)
	p3 = pyRootPwa.core.particle("pi+", True, -1, 0, 0)
	p3 = pyRootPwa.core.particle("pi+", True, -1, 0, 0, t)

	p3 = pyRootPwa.core.particle("pi-", 0, 1, 2, 3, 4, 5)
	p3.qnSummary()
	p3 = pyRootPwa.core.particle("pi-", 0, 1, 2, 3, 4, 5, 6)
	p3 = pyRootPwa.core.particle("pi-", 0, 1, 2, 3, 4, 5, 6, 7)
do_test(partTestConsts, "Testing particle other constructors")

def partTestClone(): p2 = part.clone()
do_test(partTestClone, "Testing particle.clone()")

def partTestqnSummary(): part.qnSummary()
do_test(partTestqnSummary, "Testing particle.qnSummary()")

def partTestLabel(): assert(part.label() == "Delta(1910)+[1.5(0.5+)]")
do_test(partTestLabel, "Testing particle.label()")

def partTestSpinProj():
	assert(part.spinProj == 0)
	part.spinProj = 2
	assert(part.spinProj == 2)
do_test(partTestSpinProj, "Testing particle.spinProj")

def partTestMomentum():
	assert(part.momentum.X() == 0.)
	part.momentum = pyRootPwa.ROOT.TVector3(10, 10, 10)
	assert(part.momentum.X() == 10.)
do_test(partTestMomentum, "Testing particle.momentum")

def partTestlzVec():
	assert(part.lzVec.E() == 17.42550142750560838)
	lz = pyRootPwa.ROOT.TLorentzVector(1., 1., 1., 1.)
	part.lzVec = lz
	assert(part.lzVec.E() == 1.)
do_test(partTestlzVec, "Testing particle.lzVec")

def partTestIndex():
	assert(part.index == -1)
	part.index = 12
	assert(part.index == 12)
do_test(partTestIndex, "Testing particle.index")

def partTestReflectivity():
	assert(part.reflectivity == 0)
	part.reflectivity = 2
	assert(part.reflectivity == 1)
do_test(partTestReflectivity, "Testing particle.reflectivity")

def partTestSetProperties():
	pP = pyRootPwa.core.particleProperties()
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
	old_pP_debug = pyRootPwa.core.particleProperties.debugParticleProperties
	old_debug = pyRootPwa.core.particle.debugParticle
	assert(old_pP_debug == old_debug)
	pyRootPwa.core.particle.debugParticle = (not old_debug)
	assert(pyRootPwa.core.particle.debugParticle == (not old_debug))
	assert(pyRootPwa.core.particleProperties.debugParticleProperties == old_pP_debug)
	pyRootPwa.core.particle.debugParticle = old_debug
	assert(part.debugParticle == old_debug)
	part.debugParticle = (not old_debug)
	assert(part.debugParticle == (not old_debug))
	assert(pyRootPwa.core.particleProperties.debugParticleProperties == old_pP_debug)
	part.debugParticle = old_debug
do_test(partTestDebugFlag, "Testing particle debug flag")
# This test fails for unknown reason

print
print("########################################################################")
print

# ---------------------------------------------------------
#
#	interactionVertex
#
# ---------------------------------------------------------

def defConst(): return pyRootPwa.core.interactionVertex()
iV = do_test(defConst, "Testing interactionVertex default constructor")

def copyConst(): return pyRootPwa.core.interactionVertex(iV)
do_test(copyConst, "Testing interactionVertex copy constructor")

def tClone():
	iV2 = iV.clone()
	iV2 = iV.clone(True)
	iV2 = iV.clone(True, True)
do_test(tClone, "Testing interactionVertex.clone")

def tPrint(): print("\n\n" + str(iV) + "\n")
do_test(tPrint, "Testing \"print(interactionVertex)\"")

def tClear():
	class tiV(pyRootPwa.core.interactionVertex):
		def clear(self):
			return "testString"
	iV.clear()
	iV2 = tiV(iV)
	assert(iV2.clear() == "testString")
do_test(tClear, "Testing interactionVertex.clear()")

def tAiP():
	p = pyRootPwa.core.particle()
	assert(iV.addInParticle(p))
do_test(tAiP, "Testing interactionVertex.addInParticle()")

def tAoP():
	p = pyRootPwa.core.particle()
	assert(iV.addOutParticle(p))
do_test(tAoP, "Testing interactionVertex.addOutParticle()")

def ttOP():
	rot = pyRootPwa.ROOT.TLorentzRotation(1, 1, 1)
	iV.transformOutParticles(rot)
do_test(ttOP, "Testing interactionVertex.transformOutParticles()")

def tNmbIP(): assert(iV.nmbInParticles == 1)
do_test(tNmbIP, "Testing interactionVertex.nmbInParticles")

def tNmbOP(): assert(iV.nmbOutParticles == 1)
do_test(tNmbOP, "Testing interactionVertex.nmbOutParticles")

def tInParts(): assert(len(iV.inParticles()) == 1)
do_test(tInParts, "Testing interactionVertex.inParticles()")

def tOutParts(): assert(len(iV.outParticles()) == 1)
do_test(tOutParts, "Testing interactionVertex.outParticles()")

def tName(): assert(iV.name() == "interactionVertex")
do_test(tName, "Testing interactionVertex.name()")

def tDebug():
	old_debug = iV.debugInteractionVertex
	iV.debugInteractionVertex = (not old_debug)
	assert(iV.debugInteractionVertex == (not old_debug))
	iV.debugInteractionVertex = old_debug
do_test(tDebug, "Testing interactionVertex debug flag")

print
print("########################################################################")
print

# ---------------------------------------------------------
#
#	fsVertex
#
# ---------------------------------------------------------

def fsVertexTestConsts():
	fsVert = pyRootPwa.core.fsVertex(part)
	fsVert2 = pyRootPwa.core.fsVertex(fsVert)
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
	old_iV_debug = pyRootPwa.core.interactionVertex.debugInteractionVertex
	old_fsV_debug = pyRootPwa.core.fsVertex.debugFsVertex
	pyRootPwa.core.fsVertex.debugFsVertex = (not old_fsV_debug)
	assert(pyRootPwa.core.fsVertex.debugFsVertex == (not old_fsV_debug))
	assert(pyRootPwa.core.interactionVertex.debugInteractionVertex == old_iV_debug)
	pyRootPwa.core.fsVertex.debugFsVertex = old_fsV_debug
do_test(fsVertexTestDebugFlag, "Testing fsVertex debug flag")

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
	retval = pyRootPwa.core.isobarDecayVertex(part, part, part, 0, 0, pyRootPwa.core.flatMassDependence())
	retval = pyRootPwa.core.isobarDecayVertex(part, part, part)
	return retval
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

def isobarDecVtxTestMD(): mDname = isobDecVtx.massDependence().name()
do_test(isobarDecVtxTestMD, "Testing isobarDecayVertex.massDependence()")

def isobarDecVtxTestSetMD():
	isobDecVtx.setMassDependence(pyRootPwa.core.flatMassDependence())
	assert(isobDecVtx.massDependence().name() == "flatMassDependence")
do_test(isobarDecVtxTestSetMD, "Testing isobarDecayVertex.setMassDependence()")

def isobarDecVtxTestCC():
	print("\n")
	assert(not isobDecVtx.checkConsistency())
	print
do_test(isobarDecVtxTestCC, "Testing isobarDecayVertex.checkConsistency()")

def isobarDecayVertexTestDebugFlag():
	old_iV_debug = pyRootPwa.core.interactionVertex.debugInteractionVertex
	old_fsV_debug = pyRootPwa.core.isobarDecayVertex.debugIsobarDecayVertex
	pyRootPwa.core.isobarDecayVertex.debugIsobarDecayVertex = (not old_fsV_debug)
	assert(pyRootPwa.core.isobarDecayVertex.debugIsobarDecayVertex == (not old_fsV_debug))
	assert(pyRootPwa.core.interactionVertex.debugInteractionVertex == old_iV_debug)
	pyRootPwa.core.isobarDecayVertex.debugIsobarDecayVertex = old_fsV_debug
do_test(isobarDecayVertexTestDebugFlag, "Testing isobarDecayVertex debug flag")


print
print("########################################################################")
print

# ---------------------------------------------------------
#
#	massDependence
#
# ---------------------------------------------------------

def flatMassDepTestConst(): return pyRootPwa.core.flatMassDependence()
flatMassDep = do_test(flatMassDepTestConst, "Testing flatMassDependence default constructor")

def flatMassDepTestDebug():
	old_debug = flatMassDep.debugMassDependence
	flatMassDep.debugMassDependence = (not old_debug)
	assert(flatMassDep.debugMassDependence == (not old_debug))
	flatMassDep.debugMassDependence = old_debug
do_test(flatMassDepTestDebug, "Testing massDependence debug flag.")

def flatMassDepTestAmp(): assert(flatMassDep.amp(isobDecVtx) == (1+0j))
do_test(flatMassDepTestAmp, "Testing flatMassDependence.amp()")

def flatMassDepTestName(): assert(flatMassDep.name() == "flatMassDependence")
do_test(flatMassDepTestName, "Testing flatMassDependence.name()")

def relBreitWigTestConst(): return pyRootPwa.core.relativisticBreitWigner()
relBreitWig = do_test(relBreitWigTestConst, "Testing relativisticBreitWigner default constructor")

def relBreitWigTestAmp():
	amp = relBreitWig.amp(isobDecVtx)
	zero = amp - (-0.0235304107169-0j)
	assert(zero.real < 1e-17)
	assert(zero.imag < 1e-17)
#do_test(relBreitWigTestAmp, "Testing relativisticBreitWigner.amp()")

def relBreitWigTestName(): assert(relBreitWig.name() == "relativisticBreitWigner")
do_test(relBreitWigTestName, "Testing relativisticBreitWigner.name()")

def constBreitWigTestConst(): return pyRootPwa.core.constWidthBreitWigner()
constBreitWig = do_test(constBreitWigTestConst, "Testing constWidthBreitWigner default constructor")

def constBreitWigTestAmp():
	amp = constBreitWig.amp(isobDecVtx)
	zero = amp - (-0.02351738960326379+0.00055337383635446j)
	assert(zero.real < 1e-17)
	assert(zero.imag < 1e-17)
do_test(constBreitWigTestAmp, "Testing constWidthBreitWigner.amp()")

def constBreitWigTestName(): assert(constBreitWig.name() == "constWidthBreitWigner")
do_test(constBreitWigTestName, "Testing constWidthBreitWigner.name()")

def rhoBreitWigTestConst(): return pyRootPwa.core.rhoBreitWigner()
rhoBreitWig = do_test(rhoBreitWigTestConst, "Testing rhoBreitWigner default constructor")

def rhoBreitWigTestAmp():
	amp = rhoBreitWig.amp(isobDecVtx)
	assert(math.isnan(amp.real) and math.isnan(amp.imag))
do_test(rhoBreitWigTestAmp, "Testing rhoBreitWigner.amp()")

def rhoBreitWigTestName(): assert(rhoBreitWig.name() == "rhoBreitWigner")
do_test(rhoBreitWigTestName, "Testing rhoBreitWigner.name()")

def f0980BreitWigTestConst(): return pyRootPwa.core.f0980BreitWigner()
f0980BreitWig = do_test(f0980BreitWigTestConst, "Testing f0980BreitWigner default constructor")

def f0980BreitWigTestAmp():
	amp = f0980BreitWig.amp(isobDecVtx)
	assert(math.isnan(amp.real) and math.isnan(amp.imag))
do_test(f0980BreitWigTestAmp, "Testing f0980BreitWigner.amp()")

def f0980BreitWigTestName(): assert(f0980BreitWig.name() == "f0980BreitWigner")
do_test(f0980BreitWigTestName, "Testing f0980BreitWigner.name()")

def SAuMoPenMTestConst(): return pyRootPwa.core.piPiSWaveAuMorganPenningtonM()
SAuMoPenM = do_test(SAuMoPenMTestConst, "Testing piPiSWaveAuMorganPenningtonM default constructor")

def SAuMoPenMTestAmp():
	amp = SAuMoPenM.amp(isobDecVtx)
	zero = amp - (0.00349779419823-1.21769219865e-05j)
	assert(zero.real < 1e-5)
	assert(zero.imag < 1e-5)
do_test(SAuMoPenMTestAmp, "Testing piPiSWaveAuMorganPenningtonM.amp()")

def SAuMoPenMTestName(): assert(SAuMoPenM.name() == "piPiSWaveAuMorganPenningtonM")
do_test(SAuMoPenMTestName, "Testing piPiSWaveAuMorganPenningtonM.name()")

def SAuMoPenVesTestConst(): return pyRootPwa.core.piPiSWaveAuMorganPenningtonVes()
SAuMoPenVes = do_test(SAuMoPenVesTestConst, "Testing piPiSWaveAuMorganPenningtonVes default constructor")

def SAuMoPenVesTestAmp():
	amp = SAuMoPenVes.amp(isobDecVtx)
	zero = amp - (0.00349779421172-1.2174984393e-05j)
	assert(zero.real < 1e-5)
	assert(zero.imag < 1e-5)
do_test(SAuMoPenVesTestAmp, "Testing piPiSWaveAuMorganPenningtonVes.amp()")

def SAuMoPenVesTestName(): assert(SAuMoPenVes.name() == "piPiSWaveAuMorganPenningtonVes")
do_test(SAuMoPenVesTestName, "Testing piPiSWaveAuMorganPenningtonVes.name()")

def SAuMoPenKachaevTestConst(): return pyRootPwa.core.piPiSWaveAuMorganPenningtonKachaev()
SAuMoPenKachaev = do_test(SAuMoPenKachaevTestConst, "Testing piPiSWaveAuMorganPennigtonKachaev default constructor")

def SAuMoPenKachaevTestAmp():
	amp = SAuMoPenKachaev.amp(isobDecVtx)
	zero = amp - (-0.00448845210035-2.00514420305e-05j)
	assert(zero.real < 1e-5)
	assert(zero.imag < 1e-5)
do_test(SAuMoPenKachaevTestAmp, "Testing piPiSWaveAuMorganPennigtonKachaev.amp()")

def SAuMoPenKachaevTestName(): assert(SAuMoPenKachaev.name() == "piPiSWaveAuMorganPenningtonKachaev")
do_test(SAuMoPenKachaevTestName, "Testing piPiSWaveAuMorganPennigtonKachaev.name()")

def rhoPrimeTestConst(): return pyRootPwa.core.rhoPrimeMassDep()
rhoPrime = do_test(rhoPrimeTestConst, "Testing rhoPrimeMassDep default constructor")

def rhoPrimeTestAmp():
	amp = rhoPrime.amp(isobDecVtx)
	zero = amp - (-0.00224006160232+0j)
	assert(zero.real < 1e-5)
	assert(zero.imag < 1e-5)
do_test(rhoPrimeTestAmp, "Testing rhoPrimeMassDep.amp()")

def rhoPrimeTestName(): assert(rhoPrime.name() == "rhoPrimeMassDep")
do_test(rhoPrimeTestName, "Testing rhoPrimeMassDep.name()")

print
print("########################################################################")
print

# ---------------------------------------------------------
#
#	diffractiveDissVertex
#
# ---------------------------------------------------------

def diDiVtxTestConsts():
	part2 = part.clone()
	part3 = part.clone()
	part4 = part.clone()
	v = pyRootPwa.core.diffractiveDissVertex(part, part2, part3, part4)
	v2 = pyRootPwa.core.diffractiveDissVertex(v)
	return v
diDiVtx = do_test(diDiVtxTestConsts, "Testing diffractiveDissVertex constructors")

def diDiVtxTestPrint(): print("\n\n" + str(diDiVtx) + "\n")
do_test(diDiVtxTestPrint, "Testing print(diffractiveDissVertex)")

def diDiVtxTestclone(): return diDiVtx.clone()
do_test(diDiVtxTestclone, "Testing diffractiveDissVertex.clone()")

def diDiVtxTestAddInOutPart():
	assert(not diDiVtx.addInParticle(part))
	assert(not diDiVtx.addOutParticle(part))
do_test(diDiVtxTestAddInOutPart, "Testing diffractiveDissVertex.add{In/Out}Particle()")

def diDiVtxTestrefLZVec(): assert(diDiVtx.referenceLzVec() == part.lzVec)
do_test(diDiVtxTestrefLZVec, "Testing diffractiveDissVertex.referenceLzVec()")

def diDiVtxTestXP(): assert(diDiVtx.XParticle() == part)
do_test(diDiVtxTestXP, "Testing diffractiveDissVertex.XParticle()")

def diDiVtxTestsXFXN():
	diDiVtx.XParticle().strangeness = 0
	diDiVtx.setXFlavorQN()
	assert(diDiVtx.XParticle().strangeness == 1)
do_test(diDiVtxTestsXFXN, "Testing diffractiveDissVertex.setXFlavorQN()")

def diDiVtxTestbeam():
	assert(diDiVtx.beam().name == "Kstar2(1430)2+")
	assert(diDiVtx.target().name == "Kstar2(1430)2+")
	assert(diDiVtx.recoil().name == "Kstar2(1430)2+")
do_test(diDiVtxTestbeam, "Testing diffractiveDissVertex.{beam/target/recoil}()")

def diDiVtxTestiKD():
	tCA = pyRootPwa.ROOT.TClonesArray("TObjString", 1)
	tCA[0] = pyRootPwa.ROOT.TObjString("Kstar2(1430)2+")
	assert(diDiVtx.initKinematicsData(tCA))
do_test(diDiVtxTestiKD, "Testing diffractiveDissVertex.initKinematicsData()")

def diDiVtxTestiKD2():
	print("\n")
	assert(not diDiVtx.initKinematicsData("bla"))
	print
do_test(diDiVtxTestiKD2, "Testing diffractiveDissVertex.initKinematicsData(UNSUPPORTED TYPE)")

def diDiVtxTestrKD():
	tCA = pyRootPwa.ROOT.TClonesArray("TVector3", 1)
	tCA[0] = pyRootPwa.ROOT.TVector3(12, 12, 12)
	assert(diDiVtx.readKinematicsData(tCA))
do_test(diDiVtxTestrKD, "Testing diffractiveDissVertex.readKinematicsData()")

def diDiVtxTestrKD2():
	print("\n")
	assert(not diDiVtx.readKinematicsData([0, 1, 2]))
	print
do_test(diDiVtxTestrKD2, "Testing diffractiveDissVertex.readKinematicsData(UNSUPPORTED TYPE)")

def diDiVtxTestRM(): assert(diDiVtx.revertMomenta())
do_test(diDiVtxTestRM, "Testing diffractiveDissVertex.revertMomenta()")

def diDiVtxTestName(): assert(diDiVtx.name() == "diffractiveDissVertex")
do_test(diDiVtxTestName, "Testing diffractiveDissVertex.name()")

def diDiVtxTestDebug():
	old_debug = diDiVtx.debugDiffractiveDissVertex
	diDiVtx.debugDiffractiveDissVertex = (not old_debug)
	assert(diDiVtx.debugDiffractiveDissVertex == (not old_debug))
	diDiVtx.debugDiffractiveDissVertex = old_debug
do_test(diDiVtxTestDebug, "Testing diffractiveDissVertex debug flag")

print
print("########################################################################")
print

# ---------------------------------------------------------
#
#	waveDescription
#
# ---------------------------------------------------------

def waveDescTestConst():
	return pyRootPwa.core.waveDescription()
waveDesc = do_test(waveDescTestConst, "Testing waveDescription constuctor")

def waveDescTestReadKeyFile():
	path = os.path.dirname(os.path.abspath(__file__))
	assert(waveDesc.parseKeyFile(path + "/test.key"))
	assert(waveDesc.keyFileParsed())
do_test(waveDescTestReadKeyFile, "Testing waveDescription.parseKeyFile()")

def waveDescTestConstuctTopo():
	print('\n')
	(result, topo) = waveDesc.constructDecayTopology()
	print
	assert(result)
	return topo
consistentIsobarTopo = do_test(waveDescTestConstuctTopo, "Testing waveDescription.constructDecayTopology()")

def waveDescTestConstructAmplitude():
	print('\n')
	(result, amp) = waveDesc.constructAmplitude()
	print
	assert(result)
	return amp
consistentIsobarAmp = do_test(waveDescTestConstructAmplitude, "Testing waveDescription.constructAmplitude()")

def waveDescWaveName():
	assert(pyRootPwa.core.waveDescription.waveNameFromTopology(consistentIsobarTopo) == '1-1+00+rho31690=a21320-=rho770_21_pi-_22_pi+_23_pi-')
do_test(waveDescWaveName, "Testing waveDescription::waveNameFromTopology()")

def waveDescWaveNameLatex():
	assert(pyRootPwa.core.waveDescription.waveLaTeXFromTopology(consistentIsobarTopo) == '1^{-}1^{+}0^{+}\\quad & \\rho_3^0(1690)\\rightarrow\\left\\{ a_2^-(1320)\\rightarrow\\left\\{ \\rho^0(770)\\ells{2}{1}\\pi^-\\right\\} \\ells{2}{2}\\pi^+\\right\\} \\ells{2}{3}\\pi^-')
do_test(waveDescWaveNameLatex, "Testing waveDescription::waveLaTeXFromTopology()")

def waveDescTestDebug():
	old_debug = waveDesc.debugWaveDescription
	waveDesc.debugWaveDescription = (not old_debug)
	assert(waveDesc.debugWaveDescription == (not old_debug))
	waveDesc.debugWaveDescription = old_debug
do_test(waveDescTestDebug, "Testing waveDescription debug flag")

print
print("########################################################################")
print

# ---------------------------------------------------------
#
#	decayTopology
#
# ---------------------------------------------------------

def dTTestConsts():
	t = pyRootPwa.core.decayTopology()
#	t = pyRootPwa.core.decayTopology(diDiVtx, [isobDecVtx], [part])
# There needs to be some work to produce a consistent topology
	t2 = pyRootPwa.core.decayTopology(t)
	return t
decTo = do_test(dTTestConsts, "Testing decayTopology constructors")

def dTTestPrint(): print("\n\n" + str(decTo) + "\n")
do_test(dTTestPrint, "Testing print(decayTopology)")

def dTTestClone(): return decTo.clone()
do_test(dTTestClone, "Testing decayTopology.clone()")

def dTTestnDV(): assert(decTo.nmbDecayVertices() == 0)
do_test(dTTestnDV, "Testing decayTopology.nmbDecayVertices()")

def dTTestnFP(): assert(decTo.nmbFsParticles() == 0)
do_test(dTTestnFP, "Testing decayTopology.nmbFsParticles()")

def dTTestnIFP():
	assert(decTo.nmbIndistFsParticles() == {})
	assert(consistentIsobarTopo.nmbIndistFsParticles() == {'pi+': 2, 'pi-': 3})
do_test(dTTestnIFP, "Testing decayTopology.nmbIndistFsParticles()")

def dTTestfPIP(): assert(decTo.fsParticlesIntrinsicParity() == 1)
do_test(dTTestfPIP, "Testing decayTopology.fsParticlesIntrinsicParity()")

def dTTestsIEV(): assert(consistentIsobarTopo.spaceInvEigenValue() == -1)
do_test(dTTestsIEV, "Testing decayTopology.spaceInvEigenValue()")

def dTTestrEV(): assert(consistentIsobarTopo.reflectionEigenValue() == 1)
do_test(dTTestrEV, "Testing decayTopology.reflectionEigenValue()")

def dTTestfP(): assert(decTo.fsParticles() == [])
do_test(dTTestfP, "Testing decayTopology.fsParticles()")

def dTTestdV():
	assert(decTo.decayVertices() == [])
	assert(len(consistentIsobarTopo.decayVertices()) == 4)
do_test(dTTestdV, "Testing decayTopology.decayVertices()")

def dTTestXpar(): assert(consistentIsobarTopo.XParticle().name == "X-")
do_test(dTTestXpar, "Testing decayTopology.XParticle()")

def dTTestprodVert(): assert(decTo.productionVertex() is None)
do_test(dTTestprodVert, "Testing decayTopology.productionVertex()")

def dTTestXDV():
	vtx = consistentIsobarTopo.XDecayVertex()
	assert(vtx.nmbInParticles == 1 and vtx.nmbOutParticles == 2)
do_test(dTTestXDV, "Testing decayTopology.XDecayVertex()")

def dTTesttransFSP():
	L = pyRootPwa.ROOT.TLorentzRotation(1, 1, 1)
	decTo.transformFsParticles(L)
do_test(dTTesttransFSP, "Testing decayTopology.transformFsParticles()")

def dTTesttransFSPUnsupT():
	print("\n")
	decTo.transformFsParticles([1, 2, 3])
	print
do_test(dTTesttransFSPUnsupT, "Testing decayTopology.transformFsParticles(UNSUPPORTED TYPE)")

def dTTestisbla():
	assert(not consistentIsobarTopo.isProductionVertex(isobDecVtx))
	assert(not consistentIsobarTopo.isDecayVertex(isobDecVtx))
	assert(not consistentIsobarTopo.isFsVertex(isobDecVtx))
	assert(not consistentIsobarTopo.isFsParticle(part))
	assert(consistentIsobarTopo.isDecayVertex(consistentIsobarTopo.XDecayVertex()))
	assert(consistentIsobarTopo.isProductionVertex(consistentIsobarTopo.productionVertex()))
	assert(not consistentIsobarTopo.isFsVertex(consistentIsobarTopo.decayVertices()[0]))
	assert(consistentIsobarTopo.isFsParticle(consistentIsobarTopo.fsParticles()[0]))
do_test(dTTestisbla, "Testing decayTopology.is{ProductionVertex/DecayVertex/FsVertex/FsParticle}()")

def dTTestdVInd():
	assert(consistentIsobarTopo.decayVertexIndex(consistentIsobarTopo.decayVertices()[3]) == 3)
do_test(dTTestdVInd, "Testing decayTopology.decayVertexIndex()")

def dTTestfPInd():
	assert(decTo.fsParticlesIndex(part))
do_test(dTTestfPInd, "Testing decayTopology.fsParticlesIndex()")

def tTTestCheckCons(): assert(decTo.checkConsistency())
do_test(tTTestCheckCons, "Testing decayTopology.checkConsistency()")

def tTTestAdDec():
	waste = consistentIsobarTopo.clone()
	add0r = pyRootPwa.core.isobarDecayTopology()
	waste.addDecay(add0r)
do_test(tTTestAdDec, "Testing decayTopology.addDecay()")

def tTTestSPV(): decTo.setProductionVertex(diDiVtx)
do_test(tTTestSPV, "Testing decayTopology.setProductionVertex()")

def tTTestiKD():
	tCA = pyRootPwa.ROOT.TClonesArray("TObjString", 1)
	tCA2 = pyRootPwa.ROOT.TClonesArray("TObjString", 0)
	tCA[0] = pyRootPwa.ROOT.TObjString("Kstar2(1430)2+")
	assert(decTo.initKinematicsData(tCA, tCA2))
do_test(tTTestiKD, "Testing decayTopology.initKinematicsData()")

def tTTestiKDUnsuT():
	print("\n")
	decTo.initKinematicsData("12", [])
	print
do_test(tTTestiKDUnsuT, "Testing decayTopology.initKinematicsData(UNSUPPORTED TYPE)")

def tTTestrKD():
	tCA = pyRootPwa.ROOT.TClonesArray("TVector3", 1)
	tCA[0] = pyRootPwa.ROOT.TVector3(1., 1., 1.)
	tCA2= pyRootPwa.ROOT.TClonesArray("TVector3", 0)
	assert(decTo.readKinematicsData(tCA, tCA2))
do_test(tTTestrKD, "Testing decayTopology.readKinematicsData()")

def tTTestrKDUnsT():
	print("\n")
	decTo.readKinematicsData("123", [])
	print
do_test(tTTestrKDUnsT, "Testing decayTopology.readKinematicsData(UNSUPPORTED TYPE)")

def tTTestFKD(): decTo.fillKinematicsDataCache
do_test(tTTestFKD, "Testing decayTopology.fillKinematicsDataCache()")

def tTTestRM():
	print("\n")
	assert(decTo.revertMomenta())
	assert(not decTo.revertMomenta([1,3,2]))
	try:
		decTo.revertMomenta(1)
	except TypeError:
		pass
	else:
		raise Exception("That shouldn't work.")
	print
do_test(tTTestRM, "Testing decayTopology.revertMomenta{()/(list)/(UNSUPPORTED TYPE)}")

def tTTestDebug():
	old_debug = decTo.debugDecayTopology
	decTo.debugDecayTopology = (not old_debug)
	assert(decTo.debugDecayTopology == (not old_debug))
	decTo.debugDecayTopology = old_debug
do_test(tTTestDebug, "Testing decayTopology debug flag")

def tTTestClear(): decTo.clear()
do_test(tTTestClear, "Testing decayTopology.clear()")

def tTTestToVertex():
	vertex = consistentIsobarTopo.toVertex(consistentIsobarTopo.XParticle())
	testVertex = consistentIsobarTopo.decayVertices()[0]
	assert(vertex.name() == testVertex.name())
do_test(tTTestToVertex, "Testing decayTopology.toVertex()")

def tTTestFromVertex():
	fsParticle = consistentIsobarTopo.fsParticles()[3]
	vertex = consistentIsobarTopo.fromVertex(fsParticle)
	testVertex = consistentIsobarTopo.decayVertices()[1]
	assert(vertex.name() == testVertex.name())
do_test(tTTestFromVertex, "Testing decayTopology.fromVertex()")

print
print("########################################################################")
print

# ---------------------------------------------------------
#
#	isobarDecayTopology
#
# ---------------------------------------------------------

def iDTTestConsts():
	t = pyRootPwa.core.isobarDecayTopology()
	t2 = pyRootPwa.core.isobarDecayTopology(t)
	t2 = pyRootPwa.core.isobarDecayTopology(decTo)
#	t = pyRootPwa.core.isobarDecayTopology(diDiVtx, [isobDecVtx], [part])
# Need a consistent topology for that
	return t
isoDecTop = do_test(iDTTestConsts, "Testing isobarDecayTopology constructors")

def iDTTestPrint(): print("\n\n" + str(isoDecTop) + "\n")
do_test(iDTTestPrint, "Testing print(isobarDecayTopology")

def iDTTestClone():
	t = isoDecTop.clone()
	t = isoDecTop.clone(True)
	t = isoDecTop.clone(True, True)
do_test(iDTTestClone, "Testing isobarDecayTopology.clone()")

def iDTTestiDV(): assert(isoDecTop.isobarDecayVertices() == [])
do_test(iDTTestiDV, "Testing isobarDecayTopology.isobarDecayVertices()")

def iDTTestiXDV(): assert(consistentIsobarTopo.XIsobarDecayVertex().name() == 'isobarDecayVertex')
do_test(iDTTestiXDV, "Testing isobarDecayTopology.XIsobarDecayVertex()")

def iDTTestcC(): assert(consistentIsobarTopo.checkTopology())
do_test(iDTTestcC, "Testing isobareDecayTopology.checkTopology()")

def iDTTestcCon(): assert(isoDecTop.checkConsistency())
do_test(iDTTestcCon, "Testing isobarDecayTopology.checkConsistency()")

def iDTTestSubDecay(): return(consistentIsobarTopo.subDecay(consistentIsobarTopo.XIsobarDecayVertex()))
do_test(iDTTestSubDecay, "Testing isobarDecayTopology.subDecay()")

def iDTTestJDD():
	waste = consistentIsobarTopo.clone()
	dec1 = waste.subDecay(waste.decayVertices()[3])
	dec2 = waste.subDecay(waste.decayVertices()[3])
	return pyRootPwa.core.isobarDecayTopology.joinDaughterDecays(isobDecVtx, dec1, dec2)
do_test(iDTTestJDD, "Testing isobarDecayTopology.joinDaughterDecays()")

def iDTTestCILV():
	zero = consistentIsobarTopo.calcIsobarLzVec().T() - 0.6978509
	assert(zero < 1e-5)
do_test(iDTTestCILV, "Testing isobarDecayTopology.calcIsobarLzVec()")

def iDTTestCIC(): consistentIsobarTopo.calcIsobarCharges()
do_test(iDTTestCIC, "Testing isobarDecayTopology.calcIsobarCharges()")

def iDTTestCIBN(): consistentIsobarTopo.calcIsobarBaryonNmbs()
do_test(iDTTestCIBN, "Testing isobarDecayTopology.calcIsobarBaryonNmbs()")

def iDTTestWGV():
	output = consistentIsobarTopo.writeGraphViz()
	test = 'digraph "\\"X-[1-(1+)0+]\\"" {\nnode [\nfillcolor=white, style=filled];\n0[label=diffractiveDissVertex, shape=box];\n1[label="isobarDecayVertex: L = 1, S = 0"];\n2[label="isobarDecayVertex: L = 2, S = 1"];\n3[label="isobarDecayVertex: L = 2, S = 2"];\n4[label="isobarDecayVertex: L = 2, S = 3"];\n5[label=fsVertex, shape=diamond, style="dashed,filled"];\n6[label=fsVertex, shape=diamond, style="dashed,filled"];\n7[label=fsVertex, shape=diamond, style="dashed,filled"];\n8[label=fsVertex, shape=diamond, style="dashed,filled"];\n9[label=fsVertex, shape=diamond, style="dashed,filled"];\n2 -> 1[label="rho(770)0[1+(1--)]"];\n3 -> 2[label="a2(1320)-[1-(2+)]"];\n0 -> 4[label="X-[1-(1+)0+]"];\n4 -> 3[label="rho3(1690)0[1+(3--)]"];\n4 -> 5[label="pi-[1-(0-)]"];\n3 -> 6[label="pi+[1-(0-)]"];\n2 -> 7[label="pi-[1-(0-)]"];\n1 -> 8[label="pi+[1-(0-)]"];\n1 -> 9[label="pi-[1-(0-)]"];\n}\n'
	assert(output == test)
do_test(iDTTestWGV, "Testing isobarDecayTopology.writeGraphViz()")

def iDTTestWGVtoFile():
	result = consistentIsobarTopo.writeGraphViz("test.out")
	if os.path.isfile("test.out"):
		os.system("rm -f test.out")
	else:
		assert(False)
	assert(result)
do_test(iDTTestWGVtoFile, "Testing isobarDecayTopology.writeGraphViz(outFileNmae)")

def iDTTestgIsoCGProd():
	assert(consistentIsobarTopo.getIsospinClebschGordanProduct() == -0.25000000000000006)
	assert(consistentIsobarTopo.getIsospinClebschGordanProduct(consistentIsobarTopo.isobarDecayVertices()[2]) == 0.50000000000000011)
do_test(iDTTestgIsoCGProd, "Testing isobarDecayTopology.getIsospinClebschGordanProduct()")

def iDTTestgetBoseSym():
	print
	testval = [{'fsPartPermMap': [0, 1, 2, 3, 4], 'factor': (0.2886751345948129+0j)}, {'fsPartPermMap': [0, 1, 4, 3, 2], 'factor': (0.2886751345948129+0j)}, {'fsPartPermMap': [0, 2, 1, 3, 4], 'factor': (0.2886751345948129+0j)}, {'fsPartPermMap': [0, 2, 4, 3, 1], 'factor': (0.2886751345948129+0j)}, {'fsPartPermMap': [0, 4, 1, 3, 2], 'factor': (0.2886751345948129+0j)}, {'fsPartPermMap': [0, 4, 2, 3, 1], 'factor': (0.2886751345948129+0j)}, {'fsPartPermMap': [3, 1, 2, 0, 4], 'factor': (0.2886751345948129+0j)}, {'fsPartPermMap': [3, 1, 4, 0, 2], 'factor': (0.2886751345948129+0j)}, {'fsPartPermMap': [3, 2, 1, 0, 4], 'factor': (0.2886751345948129+0j)}, {'fsPartPermMap': [3, 2, 4, 0, 1], 'factor': (0.2886751345948129+0j)}, {'fsPartPermMap': [3, 4, 1, 0, 2], 'factor': (0.2886751345948129+0j)}, {'fsPartPermMap': [3, 4, 2, 0, 1], 'factor': (0.2886751345948129+0j)}]
	retval = consistentIsobarTopo.getBoseSymmetrization()
	for itemlist in [retval, testval]:
		for item in itemlist:
			item['factor'] = complex(round(item['factor'].real, DIGITS_TO_ROUND_TO), round(item['factor'].imag, DIGITS_TO_ROUND_TO))
	assert(retval == testval)
do_test(iDTTestgetBoseSym, "Testing isobarDecayTopology.getBoseSymmetrization()")

def iDTTestgetIsoSym():
	testval = [{'fsPartPermMap': [0, 1, 2, 3, 4], 'factor': (-0.7071067811865475+0j)}, {'fsPartPermMap': [0, 1, 3, 2, 4], 'factor': (-0.7071067811865475+0j)}]
	retval = consistentIsobarTopo.getIsospinSymmetrization()
	for itemlist in [retval, testval]:
		for item in itemlist:
			item['factor'] = complex(round(item['factor'].real, DIGITS_TO_ROUND_TO), round(item['factor'].imag, DIGITS_TO_ROUND_TO))
	assert(retval == testval)
do_test(iDTTestgetIsoSym, "Testing isobarDecayTopology.getIsospinSymmetrization()")

def iDTTestisoAffPerm():
	testvals = [[False, False, False], [False, True, True], [False, False, False], [False, False, False]]
	i = 0
	for vertex in consistentIsobarTopo.isobarDecayVertices():
		assert(consistentIsobarTopo.isobarIsAffectedByPermutation(vertex, [0,1,2,3,4]) == testvals[i][0])
		assert(consistentIsobarTopo.isobarIsAffectedByPermutation(vertex, [0,1,2,4,3]) == testvals[i][1])
		assert(consistentIsobarTopo.isobarIsAffectedByPermutation(vertex, [1,0,2,4,3]) == testvals[i][2])
		i += 1
do_test(iDTTestisoAffPerm, "Testing isobarDecayTopology.isobarIsAffectedByPermutation()")

def iDTTestdausAffPerm():
	testvals = [[False, True, True], [False, True, True], [False, False, False], [False, False, False]]
	i = 0
	for vertex in consistentIsobarTopo.isobarDecayVertices():
		assert(consistentIsobarTopo.daughtersAreAffectedByPermutation(vertex, [0,1,2,3,4]) == testvals[i][0])
		assert(consistentIsobarTopo.daughtersAreAffectedByPermutation(vertex, [0,1,2,4,3]) == testvals[i][1])
		assert(consistentIsobarTopo.daughtersAreAffectedByPermutation(vertex, [1,0,2,4,3]) == testvals[i][2])
		i += 1
do_test(iDTTestdausAffPerm, "Testing isobarDecayTopology.daughtersAreAffectedByPermutation()")

def iDTTestgetFsPCTV():
	testvals = [[0,1,2,3,4], [0,1,2,3], [0,1,2], [0,1]]
	i = 0
	for vertex in consistentIsobarTopo.isobarDecayVertices():
		assert(consistentIsobarTopo.getFsPartIndicesConnectedToVertex(vertex) == testvals[i])
		i += 1
do_test(iDTTestgetFsPCTV, "Testing isobarDecayTopology.getFsPartIndicesConnectedToVertex()")

def iDTestfIsoBoseSymVerts():
	assert(consistentIsobarTopo.findIsobarBoseSymVertices() == [])
do_test(iDTestfIsoBoseSymVerts, "Testing isobarDecayTopology.findIsobarBoseSymVertices()")

def iDTTestDebug():
	old_debug = isoDecTop.debugIsobarDecayTopology
	isoDecTop.debugIsobarDecayTopology = (not old_debug)
	assert(isoDecTop.debugIsobarDecayTopology == (not old_debug))
	isoDecTop.debugIsobarDecayTopology = old_debug
do_test(iDTTestDebug, "Testing isobarDecayTopology debug flag")

def iDTTestClear(): isoDecTop.clear()
do_test(iDTTestClear, "Testing isobarDecayTopology.clear()")

print
print("########################################################################")
print

# ---------------------------------------------------------
#
#	isobarAmplitude
#
# ---------------------------------------------------------

def iATestgjTrans():
	lzVec1 = pyRootPwa.ROOT.TLorentzVector(-0.000905,-0.082895,192.513945,192.514013)
	lzVec2 = pyRootPwa.ROOT.TLorentzVector(-0.855221,0.115472,192.091726,192.100313)
	rot = pyRootPwa.core.isobarAmplitude.gjTransform(lzVec2, lzVec1)
	testval = round(-0.97198203396081984, DIGITS_TO_ROUND_GJTRAFO_TO)
	retval = round(rot.XX(), DIGITS_TO_ROUND_GJTRAFO_TO)
	assert(retval == testval)
do_test(iATestgjTrans, "Testing isobarAmplitude::gjTransform")

print
print("########################################################################")
print

# ---------------------------------------------------------
#
#	isobarCanonicalAmplitude
#
# ---------------------------------------------------------

def iCATestConst():
	t = pyRootPwa.core.isobarCanonicalAmplitude()
	# t2 = pyRootPwa.core.isobarCanonicalAmplitude(isoDecTop)
	# needs consistent topology
	return t
iCA = do_test(iCATestConst, "Testing isobarCanonicalAmplitude constructors")

def iCATestName(): assert(iCA.name() == "isobarCanonicalAmplitude")
do_test(iCATestName, "Testing isobarCanonicalAmplitude.name()")

def iCATestDebug():
	old_debug = iCA.debugIsobarCanonicalAmplitude
	iCA.debugIsobarCanonicalAmplitude = (not old_debug)
	assert(iCA.debugIsobarCanonicalAmplitude == (not old_debug))
	iCA.debugIsobarCanonicalAmplitude = old_debug
do_test(iCATestDebug, "Testing isobarCanonicalAmplitude debug flag")

print
print("########################################################################")
print

# ---------------------------------------------------------
#
#	isobarHelicityAmplitude
#
# ---------------------------------------------------------

def iHATestConst():
	t = pyRootPwa.core.isobarHelicityAmplitude()
	# t2 = pyRootPwa.core.isobarHelicityAmplitude(isoDecTop)
	# needs consistent topology
	return t
iHA = do_test(iHATestConst, "Testing isobarHelicityAmplitude constructors")

def iHATestName(): assert(iHA.name() == "isobarHelicityAmplitude")
do_test(iHATestName, "Testing isobarHelicityAmplitude.name()")

def iHATTesthfTransform():
	vec = pyRootPwa.ROOT.TLorentzVector(0., 1., 2., 3.)
	rot1 = pyRootPwa.core.isobarHelicityAmplitude.hfTransform(vec)
	assert(rot1.XY() == 0.89442719099991586)
	rot2 = iHA.hfTransform(vec)
	assert(rot1 == rot2)
do_test(iHATTesthfTransform, "Testing isobarHelicityAmplitude::hfTransform")

def iHATestDebug():
	old_debug = iHA.debugIsobarHelicityAmplitude
	iHA.debugIsobarHelicityAmplitude = (not old_debug)
	assert(iHA.debugIsobarHelicityAmplitude == (not old_debug))
	iHA.debugIsobarHelicityAmplitude = old_debug
do_test(iHATestDebug, "Testing isobarHelicityAmplitude debug flag")

print
print("########################################################################")
print

# ---------------------------------------------------------
#
#	generatorParameters
#
# ---------------------------------------------------------

def beamTestConst():
	b = pyRootPwa.core.Beam()
	return b
beam = do_test(beamTestConst, "Testing Beam constructor")

def beamTestParticle():
	beam.particle = part
	assert(part == beam.particle)
do_test(beamTestParticle, "Testing Beam.particle data member")

def beamTestMomentum():
	beam.momentum = 123.5634
	assert(beam.momentum == 123.5634)
do_test(beamTestMomentum, "Testing Beam.momentum data member")

def beamTestMomentumSigma():
	beam.momentumSigma = 123.5634
	assert(beam.momentumSigma == 123.5634)
do_test(beamTestMomentumSigma, "Testing Beam.momentumSigma data member")

def beamTestDxDz():
	beam.DxDz = 123.5634
	assert(beam.DxDz == 123.5634)
do_test(beamTestDxDz, "Testing Beam.DxDz data member")

def beamTestDxDzSigma():
	beam.DxDzSigma = 123.5634
	assert(beam.DxDzSigma == 123.5634)
do_test(beamTestDxDzSigma, "Testing Beam.DxDzSigma data member")

def beamTestDyDz():
	beam.DyDz = 123.5634
	assert(beam.DyDz == 123.5634)
do_test(beamTestDyDz, "Testing Beam.DyDz data member")

def beamTestDyDzSigma():
	beam.DyDzSigma = 123.5634
	assert(beam.DyDzSigma == 123.5634)
do_test(beamTestDyDzSigma, "Testing Beam.DyDzSigma data member")

def beamTestPrint():
	print("\n")
	print(beam)
do_test(beamTestPrint, "Testing print(Beam)")

def targetTestConst():
	target = pyRootPwa.core.Target()
	return target
target = do_test(targetTestConst, "Testing Target constructor")

def targetTestTPart():
	target.targetParticle = part
	assert(target.targetParticle == part)
do_test(targetTestTPart, "Testing Target.targetParticle data member")

def targetTestRPart():
	target.recoilParticle = part
	assert(target.recoilParticle == part)
do_test(targetTestRPart, "Testing Target.recoilParticle data member")

def targetTestPos():
	vec = pyRootPwa.ROOT.TVector3(1., 2., 3.)
	target.position = vec
	assert(target.position == vec)
do_test(targetTestPos, "Testing Target.position data member")

def targetTestLength():
	target.length = 1434.2313
	assert(target.length == 1434.2313)
do_test(targetTestLength, "Testing Target.length data member")

def targetTestRadius():
	target.radius = 1434.2313
	assert(target.radius == 1434.2313)
do_test(targetTestRadius, "Testing Target.radius data member")

def targetTestIntLength():
	target.interactionLength = 1434.2313
	assert(target.interactionLength == 1434.2313)
do_test(targetTestIntLength, "Testing Target.interactionLength data member")

def targetTestPrint():
	print("\n")
	print(target)
do_test(targetTestPrint, "Testing print(Target)")

def finalStTestConst():
	fs = pyRootPwa.core.FinalState()
	return fs
finalState = do_test(finalStTestConst, "Testing FinalState constructor")

def finalStTestParts():
	l = [part, part, part]
	finalState.particles = l
do_test(finalStTestParts, "Testing FinalState.particles data member")

def finalStTestPrint():
	print("\n")
	print(finalState)
do_test(finalStTestPrint, "Testing print(FinalState)")

# ---------------------------------------------------------
#
#	Summary
#
# ---------------------------------------------------------

print
if errors == 0:
	print_green("All tests successful.")
elif errors == 1:
	print_red("There was " + str(errors) + " error.")
else:
	print_red("There were " + str(errors) + " errors.")
if skip > 0:
	if skip == 1:
		outstring = " test was skipped."
	else:
		outstring = " tests were skipped."
	print_yellow(str(skip) + outstring)

sys.exit(errors)
