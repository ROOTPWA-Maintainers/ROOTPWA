#!/usr/bin/env python

import argparse
import sys
import multiprocessing

import pyRootPwa
import pyRootPwa.utils

if __name__ == "__main__":

    # parse command line arguments
    parser = argparse.ArgumentParser(
         description="converts ASCII *.evt files "
                     "to ROOT tree file "
    )

    parser.add_argument("-e", type=str, metavar="evt-file", dest="evtFileName", default='/nfs/mnemosyne/user/odrotleff/private/proj/rootpwaMaster/ROOTPWA/DATA/1120.1190/1120.1190.evt', help="path to evt file")
    parser.add_argument("-o", type=str, metavar="output-file", dest="outFileName", default="./test.root", help="path to output file")
    parser.add_argument("-n", type=int, metavar="event number", default=-1, dest="numEvents", help="")
    parser.add_argument("-a", type=str, metavar="treename", dest="treeName", default="rootPwaEvtTree", help="")
    parser.add_argument("-b", type=str, metavar="prodKinParticlesName", dest="prodKinParticlesName", default="prodKinParticlesName", help="")
    parser.add_argument("-c", type=str, metavar="prodKinMomentaLeafName", dest="prodKinMomentaLeafName", default="prodKinMomentaLeafName", help="")
    parser.add_argument("-d", type=str, metavar="decayKinParticlesName", dest="decayKinParticlesName", default="decayKinParticlesName", help="")
    parser.add_argument("-f", type=str, metavar="decayKinMomentaLeafName", dest="decayKinMomentaLeafName", default="decayKinMomentaLeafName", help="")
    parser.add_argument("-g", action='store_true', dest="debug", default=True)

    args = parser.parse_args()

    # print some info
    pyRootPwa.core.printCompilerInfo()
    pyRootPwa.core.printLibraryInfo()
    pyRootPwa.core.printGitHash()

    pyRootPwa.utils.stdoutisatty = sys.stdout.isatty()
    pyRootPwa.utils.stderrisatty = sys.stderr.isatty()

    # initialize the printing functors
    printingCounter = multiprocessing.Array('i', [0] * 5)
    pyRootPwa.utils.printErr = pyRootPwa.utils.printErrClass(printingCounter)
    pyRootPwa.utils.printWarn = pyRootPwa.utils.printWarnClass(printingCounter)
    pyRootPwa.utils.printSucc = pyRootPwa.utils.printSuccClass(printingCounter)
    pyRootPwa.utils.printInfo = pyRootPwa.utils.printInfoClass(printingCounter)
    pyRootPwa.utils.printDebug = pyRootPwa.utils.printDebugClass(printingCounter)

    printErr = pyRootPwa.utils.printErr
    printWarn = pyRootPwa.utils.printWarn
    printSucc = pyRootPwa.utils.printSucc
    printInfo = pyRootPwa.utils.printInfo
    printDebug = pyRootPwa.utils.printDebug

    outputFile = pyRootPwa.ROOT.TFile.Open(args.outFileName, "NEW")
    outTree = pyRootPwa.ROOT.TTree(args.treeName, args.treeName)
    debug = args.debug

    # create leaf variables
    prodKinMomenta = pyRootPwa.ROOT.TClonesArray("TVector3")
    decayKinMomenta = pyRootPwa.ROOT.TClonesArray("TVector3")

    # connect leaf variables to tree branches
    splitLevel = 99
    bufSize = 256000
    outTree.Branch(args.prodKinMomentaLeafName, "TClonesArray", prodKinMomenta, bufSize, splitLevel)
    outTree.Branch(args.decayKinMomentaLeafName, "TClonesArray", decayKinMomenta, bufSize, splitLevel)

    # loop over events and fill tree
    prodNames = []
    decayNames = []
    countEvents = 0
    countLines = 0

    with open(args.evtFileName) as evtFile:
        for line in evtFile:

            # read number of particles
            nmbParticles = 0
            columns = line.split()

            if len(columns) == 1:
                nmbParticles = columns[0]
            else:
                printErr("Should only be one column in first line of event")

            if debug:
                printDebug("# of particles = " + nmbParticles)

            # read production kinematics data (beam + fixed target)
            prodNames = []
            prodKinMomenta.Clear()

            line = next(evtFile)
            columns = line.split()

            (id, charge, momX, momY, momZ, E) = (columns)
            if debug:
                printDebug(id + " " + charge + " " + momX + " " + momY + " " + momZ + " " + E)

            partName = particleNameFromGeantId(id)
            prodNames.append(partName)

            if not checkParticleCharge(id, charge):
                success = false;

            prodKinMomenta.AddLast(pyRootPwa.ROOT.TVector3(momX, momY, momZ))

            # check consistency

            # read decay kinematics data
            decayNames = []
            decayKinMomenta.Clear()

            for i in range(0, nmbParticles - 2):
                line = next(evtFile)
                columns = line.split()

                (id, charge, momX, momY, momZ, E) = (columns)
                if debug:
                    printDebug(id + " " + charge + " " + momX + " " + momY + " " + momZ + " " + E)
    
                partName = particleNameFromGeantId(id)
                decayNames.append(partName)
                
                if not checkParticleCharge(id, charge):
                    success = false;
    
                decayKinMomenta.AddLast(pyRootPwa.ROOT.TVector3(momX, momY, momZ))

            # check consistency
            
            outTree.Fill()
        outTree.Write()
        outputFile.Close()













