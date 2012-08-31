
print
print("Importing pyRootPwa")
print("---------------------------------")
print

import pyRootPwa

ROOT = pyRootPwa.ROOT

print
print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
print

print
print("filling particleDataTable")
print("---------------------------------")
print

pyRootPwa.particleDataTable.readFile("../../amplitude/particleDataTable.txt")

print
print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
print

print
print("Instantiated particle")
print("---------------------------------")
print
p = pyRootPwa.particle()

print
print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
print

print
print("p.read(\"\")")
print("---------------------------------")
print

p.read("Delta(1910)+     Deltabar(1910)-  1.91        0.25        +1           3    0    0    0     0     1    +1     0")

print
print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
print

print
print("print(p)")
print("---------------------------------")
print

print(p)

print
print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
print

print
print("Testing copy constructor -> creating p2")
print("---------------------------------")
print

p2 = pyRootPwa.particle(p)

print
print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
print

print
print("print(p2)")
print("---------------------------------")
print

print(p2)

print
print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
print

print
print("Testing constructors")
print("---------------------------------")
print

pP = pyRootPwa.particleProperties()
pP2 = pyRootPwa.particleProperties(pP)
pP.read("Delta(1910)+     Deltabar(1910)-  1.91        0.25        +1           3    0    0    0     0     1    +1     0")
t = ROOT.TVector3(1., 1., 1.)

p3 = pyRootPwa.particle(pP)
p3 = pyRootPwa.particle(pP, -1)
p3 = pyRootPwa.particle(pP, -1, 0)
p3 = pyRootPwa.particle(pP, -1, 0, 0)
p3.momentum.Print()
p3 = pyRootPwa.particle(pP, -1, 0, 0, t)
p3.momentum.Print()

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

print
print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
print

print
print("p.clone()")
print("---------------------------------")
print

p2 = p.clone()
print(p2)
p2.isospin = 5
print(p2.isospin)
print(p.isospin)

print
print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
print


print
print("p.qnSummary()")
print("---------------------------------")
print

print(p.qnSummary())

print
print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
print

print
print("p.label()")
print("---------------------------------")
print

print(p.label())

print
print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
print

print
print("p.momentum")
print("---------------------------------")
print

print(p.momentum)
print(p.momentum.X())
print(p.momentum.Y())
print(p.momentum.Z())
p.momentum.Print()
p.momentum = ROOT.TVector3(10, 10, 10)
p.momentum.Print()
p.momentum = "uiaeiaeuiae"
p.momentum.Print()

print
print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
print

print
print("p.transform(TVector3)")
print("---------------------------------")
print

tv3 = ROOT.TVector3(1., 1., 1.)
print(p.transform(tv3))

print
print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
print

print
print("p.transform(TLorentzVector)")
print("---------------------------------")
print

t = ROOT.TLorentzRotation()
print(p.transform(t))

print
print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
print

print
print("p.transform(<UNSUPPORTED_TYPE>)")
print("vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv")
print

first = False
second = False

t = ROOT.TH1D()
if p.transform(t) is None:
	first = True

t = [0, 1, 2]
if p.transform(t) is None:
	second = True

if not (first and second):
	raise Exception("p.transform(<UNSUPPORTED_TYPE>) did not behave as expected!")

print
print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
print

