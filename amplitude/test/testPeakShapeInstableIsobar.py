#!/usr/bin/env python

import math
import multiprocessing

import pyRootPwa
import pyRootPwa.utils
ROOT = pyRootPwa.ROOT

#initialize the printing functors
printingCounter = multiprocessing.Array('i', [0]*5)
pyRootPwa.utils.printErr = pyRootPwa.utils.printErrClass(printingCounter)
pyRootPwa.utils.printWarn = pyRootPwa.utils.printWarnClass(printingCounter)
pyRootPwa.utils.printSucc = pyRootPwa.utils.printSuccClass(printingCounter)
pyRootPwa.utils.printInfo = pyRootPwa.utils.printInfoClass(printingCounter)
pyRootPwa.utils.printDebug = pyRootPwa.utils.printDebugClass(printingCounter)


def getPhaseSpaceVolume(n, M, m):
	return (4.*math.pi)**(n-1) * (1. / math.factorial(n-2)) * ( M - n*m )**(n-2)


def integrate(subAmp, nEvents, parentMass):

	pyRootPwa.utils.printInfo("calculating integral for M=" + str(parentMass) + " with " + str(nEvents) + " events.")

	prodKinNames = ROOT.TClonesArray("TObjString", 1)
	decayKinNames = ROOT.TClonesArray("TObjString", 3)

	prodKinNames[0] = ROOT.TObjString("a2(1320)-")
	decayKinNames[0] = ROOT.TObjString("pi-")
	decayKinNames[1] = ROOT.TObjString("pi+")
	decayKinNames[2] = ROOT.TObjString("pi-")

	subAmp.decayTopology().initKinematicsData(prodKinNames, decayKinNames)
	subAmp.reflectivityBasis = False
	subAmp.init()

	#pionMass = pyRootPwa.core.particleDataTable.entry('pi+').mass
	#parentMass = 1.3

	pionMass = 0.13957
#	parentMass = 5.

	psGen = pyRootPwa.core.nBodyPhaseSpaceGen()
	psGen.setDecay([pionMass, pionMass, pionMass])
	psGen.setMaxWeight(1.01 * psGen.estimateMaxWeight(parentMass, 100000))

	parent = ROOT.TLorentzVector(0, 0, 0, parentMass)

	integral = 0.

	prog = pyRootPwa.utils.progressBar(0, nEvents)
	prog.start()

#	g2d = ROOT.TGraph2D(nEvents)
	g2d = ROOT.TH2D("2d_"+str(parentMass), "2d_"+str(parentMass), 750, 0, 5, 750, 0, 5)

	for i in range(nEvents):

		psGen.generateDecay(parent)

		prodKinM = ROOT.TClonesArray("TVector3", 1)
		decayKinM = ROOT.TClonesArray("TVector3", 3)

		decayKinM[0] = psGen.daughter(0).Vect()
		decayKinM[1] = psGen.daughter(1).Vect()
		decayKinM[2] = psGen.daughter(2).Vect()

		prodKinM[0] = (psGen.daughter(0) + psGen.daughter(1) + psGen.daughter(2)).Vect()


		subAmp.decayTopology().readKinematicsData(prodKinM, decayKinM)
		eventWeight = psGen.eventWeight()

		ampVal = abs(subAmp())
		integral += ampVal * ampVal * eventWeight
#		integral += eventWeight

#		g2d.SetPoint(i, (psGen.daughter(0)+psGen.daughter(1)).M2(), (psGen.daughter(1)+psGen.daughter(2)).M2(), ampVal*ampVal)
#		g2d.Fill((psGen.daughter(0)+psGen.daughter(1)).M2(), (psGen.daughter(1)+psGen.daughter(2)).M2(), ampVal*ampVal)

		prog.update(i)

	integral /= float(nEvents)
	integral *= getPhaseSpaceVolume(3, parentMass, pionMass)

	g2d.Write()

	error = 0.
	return (integral, error)


pyRootPwa.core.particleDataTable.readFile('/opt/rootpwa/amplitude/particleDataTable.txt')

keyfile = '/home/kbicker/analysis/data_horsing_around/keys/1-2-00+rho770_02_a21320-=rho770_21_pi-.key'

waveDesc = pyRootPwa.core.waveDescription()
waveDesc.parseKeyFile(keyfile)
(result, amplitude) = waveDesc.constructAmplitude()
amplitude.init()
topo = amplitude.decayTopology()

#print(amplitude)
#print(topo)

decayVertices = topo.decayVertices()

v = None

for vertex in decayVertices:
	if vertex.parent().label() == "a2(1320)-[1-(2+)]":
		v = vertex

#print(v)

v.setMassDependence(pyRootPwa.core.flatMassDependence())
subTopo = topo.subDecayConsistent(v)

#print(subTopo)

subAmp = pyRootPwa.core.isobarHelicityAmplitude(subTopo)

#print(subAmp)

outfile = ROOT.TFile.Open("integralTest.root", "RECREATE")

gr = ROOT.TGraphErrors(10)
gr.SetTitle("int")
gr.SetName("int")

g = ROOT.TGraphErrors(10)
g.SetTitle("dyn")
g.SetName("dyn")

gbw = ROOT.TGraphErrors(10)
gbw.SetTitle("bw")
gbw.SetName("bw")

nEvents = 20000

M0     = 1.35;
Gamma0 = 0.35;

(int0, error0) = integrate(subAmp, 100000, pyRootPwa.core.particleDataTable.entry('a2(1320)-').mass)
#(int0, error0) = integrate(subAmp, nEvents, 1.32)
pyRootPwa.utils.printSucc("calculated int0: " + str(int0))

pyRootPwa.core.randomNumberGenerator.instance.setSeed(987654321)

if True:

	for i in range(20):
#	for i in range(10):
#	for i in range(1):

		M = 0.5 + i * 0.175
#		M = 0.5 + i * 0.2
#		M = 2.3

		(integral, error) = integrate(subAmp, nEvents, M)
		pyRootPwa.utils.printSucc("calculated integral: " + str(integral))
		dyn = integral / int0
		gr.SetPoint(i, M, integral)
		g.SetPoint(i, M, dyn)
		gbw.SetPoint(i, M, abs((M0 * Gamma0) / ((M0*M0) - (M*M) - (0+1j)*M0*dyn*Gamma0)))

gr.Write()
g.Write()
gbw.Write()
outfile.Write()
outfile.Close()
