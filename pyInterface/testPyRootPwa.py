
import math
import os
import sys

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

def getpDTinstance(): return pyRootPwa.core.particleDataTable.instance
particleTable = do_test( getpDTinstance, "Testing particleDataTable.instance")

def tPDTreadFile():
	print("\n")
	pyRootPwa.core.particleDataTable.readFile(os.environ['ROOTPWA'] + "/amplitude/particleDataTable.txt")
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
do_test_raw(partTestConsts, "Testing particle other constructors")

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
#	retval = pyRootPwa.core.isobarDecayVertex(part, part, part, 0, 0, pyRootPwa.core.flatMassDependence())
#	retval = pyRootPwa.core.isobarDecayVertex(part, part, part, 0, 0, [0, 1, 2])
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
do_test(flatMassDepTestDebug, "Testing flatMassDependence debug flag.")

def flatMassDepTestAmp(): assert(flatMassDep.amp(isobDecVtx) == (1+0j))
do_test(flatMassDepTestAmp, "Testing flatMassDependence.amp()")

def flatMassDepTestName(): assert(flatMassDep.name() == "flatMassDependence")
do_test(flatMassDepTestName, "Testing flatMassDependence.name()")

def relBreitWigTestConst(): return pyRootPwa.core.relativisticBreitWigner()
relBreitWig = do_test(relBreitWigTestConst, "Testing relativisticBreitWigner default constructor")

def relBreitWigTestAmp():
	amp = relBreitWig.amp(isobDecVtx)
	assert(math.isnan(amp.real) and math.isnan(amp.imag))
do_test(relBreitWigTestAmp, "Testing relativisticBreitWigner.amp()")

def relBreitWigTestName(): assert(relBreitWig.name() == "relativisticBreitWigner")
do_test(relBreitWigTestName, "Testing relativisticBreitWigner.name()")

def SAuMoPenMTestConst(): return pyRootPwa.core.piPiSWaveAuMorganPenningtonM()
SAuMoPenM = do_test(SAuMoPenMTestConst, "Testing piPiSWaveAuMorganPenningtonM default constructor")

def SAuMoPenMTestAmp():
	amp = SAuMoPenM.amp(isobDecVtx)
	zero = amp - (0.00349779419823-1.21769219865e-05j)
	assert(zero.real < 10e-15)
	assert(zero.imag < 10e-15)
do_test(SAuMoPenMTestAmp, "Testing piPiSWaveAuMorganPenningtonM.amp()")

def SAuMoPenMTestName(): assert(SAuMoPenM.name() == "piPiSWaveAuMorganPenningtonM")
do_test(SAuMoPenMTestName, "Testing piPiSWaveAuMorganPenningtonM.name()")

def SAuMoPenVesTestConst(): return pyRootPwa.core.piPiSWaveAuMorganPenningtonVes()
SAuMoPenVes = do_test(SAuMoPenVesTestConst, "Testing piPiSWaveAuMorganPenningtonVes default constructor")

def SAuMoPenVesTestAmp():
	amp = SAuMoPenVes.amp(isobDecVtx)
	zero = amp - (0.00349779421172-1.2174984393e-05j)
	assert(zero.real < 10e-15)
	assert(zero.imag < 10e-15)
do_test(SAuMoPenVesTestAmp, "Testing piPiSWaveAuMorganPenningtonVes.amp()")

def SAuMoPenVesTestName(): assert(SAuMoPenVes.name() == "piPiSWaveAuMorganPenningtonVes")
do_test(SAuMoPenVesTestName, "Testing piPiSWaveAuMorganPenningtonVes.name()")

def SAuMoPenKachaevTestConst(): return pyRootPwa.core.piPiSWaveAuMorganPenningtonKachaev()
SAuMoPenKachaev = do_test(SAuMoPenKachaevTestConst, "Testing piPiSWaveAuMorganPennigtonKachaev default constructor")

def SAuMoPenKachaevTestAmp():
	amp = SAuMoPenKachaev.amp(isobDecVtx)
	zero = amp - (-0.00448845210035-2.00514420305e-05j)
	assert(zero.real < 10e-15)
	assert(zero.imag < 10e-15)
do_test(SAuMoPenKachaevTestAmp, "Testing piPiSWaveAuMorganPennigtonKachaev.amp()")

def SAuMoPenKachaevTestName(): assert(SAuMoPenKachaev.name() == "piPiSWaveAuMorganPenningtonKachaev")
do_test(SAuMoPenKachaevTestName, "Testing piPiSWaveAuMorganPennigtonKachaev.name()")

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

def dTTestnIFP(): assert(decTo.nmbIndistFsParticles() == {})
do_test(dTTestnIFP, "Testing decayTopology.nmbIndistFsParticles()")

def dTTestfPIP(): assert(decTo.fsParticlesIntrinsicParity() == 1)
do_test(dTTestfPIP, "Testing decayTopology.fsParticlesIntrinsicParity()")

def dTTestsIEV(): print(decTo.spaceInvEigenValue())
do_test(dTTestsIEV, "Testing decayTopology.spaceInvEigenValue()", True)
# needs a consistent topology

def dTTestrEV(): print(decTo.reflectionEigenValue())
do_test(dTTestrEV, "Testing decayTopology.reflectionEigenValue()", True)
# needs a consistent topology

def dTTestfP(): assert(decTo.fsParticles() == [])
do_test(dTTestfP, "Testing decayTopology.fsParticles()")

def dTTestdV(): assert(decTo.decayVertices() == [])
do_test(dTTestdV, "Testing decayTopology.decayVertices()")

def dTTestXpar(): print(decTo.XParticle())
do_test(dTTestXpar, "Testing decayTopology.XParticle()", True)
# needs a consistent topology

def dTTestprodVert(): assert(decTo.productionVertex() is None)
do_test(dTTestprodVert, "Testing decayTopology.productionVertex()")

def dTTestXDV(): print(decTo.XDecayVertex())
do_test(dTTestXDV, "Testing decayTopology.XDecayVertex()", True)
# needs a consistent topology

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
	assert(not decTo.isProductionVertex(isobDecVtx))
	assert(decTo.isDecayVertex(isobDecVtx))
	assert(not decTo.isFsVertex(isobDecVtx))
	assert(not decTo.isFsParticle(part))
do_test(dTTestisbla, "Testing decayTopology.is{ProductionVertex/DecayVertex/FsVertex/FsParticle}()", True)
# needs a consistent topology

def dTTestfPInd():
	assert(decTo.fsParticlesIndex(part))
do_test(dTTestfPInd, "Testing decayTopology.fsParticlesIndex()")

def tTTestCheckbla(): print(not decTo.checkTopology())
do_test(tTTestCheckbla, "Testing decayTopology.checkTopology()", True)
# needs a consistent topology

def tTTestCheckCons(): assert(decTo.checkConsistency())
do_test(tTTestCheckCons, "Testing decayTopology.checkConsistency()")

def tTTestAdDec(): decTo.addDecay(decTo)
do_test(tTTestAdDec, "Testing decayTopology.addDecay()", True)
# needs a consistent topology

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
	assert(decTo.revertMomenta())
	assert(not decTo.revertMomenta([1,3,2]))
	try:
		decTo.revertMomenta(1)
	except TypeError:
		pass
	else:
		raise Exception("That shouldn't work.")
do_test(tTTestRM, "Testing decayTopology.revertMomenta{()/(list)/(UNSUPPORTED TYPE)}")

def tTTestDebug():
	old_debug = decTo.debugDecayTopology
	decTo.debugDecayTopology = (not old_debug)
	assert(decTo.debugDecayTopology == (not old_debug))
	decTo.debugDecayTopology = old_debug
do_test(tTTestDebug, "Testing decayTopology debug flag")

def tTTestClear(): decTo.clear()
do_test(tTTestClear, "Testing decayTopology.clear()")

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

def iDTTestiXDV(): print(isoDecTop.XIsobarDecayVertex())
do_test(iDTTestiXDV, "Testing isobarDecayTopology.XIsobarDecayVertex()", True)
# need a consistent topology for that

def iDTTestcC(): assert(isoDecTop.checkTopology())
do_test(iDTTestcC, "Testing isobareDecayTopology.checkTopology()", True)

def iDTTestcCon(): assert(isoDecTop.checkConsistency())
do_test(iDTTestcCon, "Testing isobarDecayTopology.checkConsistency()")

def iDTTestAddDec(): isoDecTop.addDecay(isoDecTop)
do_test(iDTTestAddDec, "Testing isobarDecayTopology.addDecay()", True)
# need a consistent topology for that

def iDTTestJDD(): pyRootPwa.core.isobarDecayTopology.joinDaughterDecays(isobDecVtx, isoDecTop, isoDecTop)
do_test(iDTTestJDD, "Testing isobarDecayTopology.joinDaughterDecays()", True)
# need a consistent topology for that

def iDTTestCILV(): print(isoDecTop.calcIsobarLzVec())
do_test(iDTTestCILV, "Tessting isobarDecayTopology.calcIsobarLzVec()", True)
# need a consistent topology for that

def iDTTestCIC(): print(isoDecTop.calcIsobarCharges())
do_test(iDTTestCIC, "Testing isobarDecayTopology.calcIsobarCharges()", True)
# need a consistent topology for that

def iDTTestCIBN(): print(isoDecTop.calcIsobarBaryonNmbs())
do_test(iDTTestCIC, "Testing isobarDecayTopology.calcIsobarBaryonNmbs()", True)
# need a consistent topology for that

def iDTTestWGV(): print(isoDecTop.writeGraphViz())
do_test(iDTTestWGV, "Testing isobarDecayTopology.writeGraphViz()", True)
# need a consistent topology for that

def iDTTestWGVtoFile():
	print(isoDecTop.writeGraphViz("test.out"))
	# Check if the file is there and if it is, remove it
do_test(iDTTestWGVtoFile, "Testing isobarDecayTopology.writeGraphViz(outFileNmae)", True)
# need a consistent topology for that

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
	rot2 = iHA.hfTransform(vec)
	assert(rot1 == rot2)
do_test_raw(iHATTesthfTransform, "Testing isobarHelicityAmplitude::hfTransform")

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
