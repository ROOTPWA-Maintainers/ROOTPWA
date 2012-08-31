
import pyRootPwa

ROOT = pyRootPwa.ROOT

print("Imported pyRootPwa")
print("---------------------------------")
print("")

p = pyRootPwa.particle()

print("Instantiated particle")
print("---------------------------------")
print("")

p.read("Delta(1910)+     Deltabar(1910)-  1.91        0.25        +1           3    0    0    0     0     1    +1     0")

print("p.read(\"\") executed")
print("---------------------------------")
print("")

print(p)


print("print(p) executed")
print("---------------------------------")
print("")

print(p.qnSummary())

print("p.qnSummary() executed")
print("---------------------------------")
print("")

print(p.label())

print("p.label() executed")
print("---------------------------------")
print("")

print(p.momentum())
print(p.momentum().X())
print(p.momentum().Y())
print(p.momentum().Z())
print(p.momentum().Print())

print("p.momentum() executed")
print("---------------------------------")
print("")

tv3 = ROOT.TVector3(1., 1., 1.)
print(p.transform(tv3))

print("p.transform(TVector3) executed")
print("---------------------------------")
print("")

t = ROOT.TLorentzRotation()
print(p.transform(t))

print("p.transform(TLorentzVector) executed")
print("---------------------------------")
print("")

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

print("p.transform(<UNSUPPORTED_TYPE>) executed")
print("---------------------------------")
print("")

