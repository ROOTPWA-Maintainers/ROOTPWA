
import sys

success = True

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

def impLib(): import pyRootPwa
do_test(impLib, "Importing pyRootPwa")

import pyRootPwa

def defConst(): return pyRootPwa.interactionVertex()
iV = do_test(defConst, "Testing default constructor")

def copyConst(): return pyRootPwa.interactionVertex(iV)
do_test(copyConst, "Testing copy constructor")

def tClone():
	iV2 = iV.clone()
	iV2 = iV.clone(True)
	iV2 = iV.clone(True, True)
do_test(tClone, "Testing interactionVertex::clone")

def tPrint(): print("\n\n" + str(iV))
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
if success:
	print_green("All tests successful.")
else:
	print_red("There were errors.")
