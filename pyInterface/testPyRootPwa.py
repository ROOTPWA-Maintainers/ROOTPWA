
import sys

success = True

# Some functions
# ---------------------------------------------------------

def print_green(string):
	print('\033[92m' + string + '\033[0m')

def print_red(string):
	print('\033[91m' + string + '\033[0m')

def do_test(function, name):
	sys.stdout.write(name + "...")
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
#	particleDataTable
#
# ---------------------------------------------------------

def partPropTestConst(): return pyRootPwa.particleProperties()
partProp = do_test(partPropTestConst, "Testing particleProperties constructor")

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
