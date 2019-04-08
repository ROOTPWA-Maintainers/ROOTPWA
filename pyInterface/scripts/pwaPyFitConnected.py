#!/usr/bin/env python
# pylint: disable=C0413,W0122
import os
# disable multithreading by default
os.environ['OPENBLAS_NUM_THREADS'] = "1"
import argparse
import sys
import imp

import pyRootPwa
import pyRootPwa.pyPartialWaveFit

def buildBinRange(args):

	config = pyRootPwa.rootPwaConfig()
	if not config.initialize(args.configFileName):
		pyRootPwa.utils.printErr("loading config file '" + args.configFileName + "' failed. Aborting...")
		sys.exit(1)
	fileManager = pyRootPwa.loadFileManager(config.fileManagerPath)
	if not fileManager:
		pyRootPwa.utils.printErr("loading the file manager failed. Aborting...")
		sys.exit(1)

	centralBin = fileManager.binList[args.integralBin]
	otherMassBins = sorted([ b for b in fileManager.binList if b in centralBin.getSubMultiBin(exception='mass')])
	iCentralBin = otherMassBins.index(centralBin)
	nOtherMassBins = len(otherMassBins)
	iBinStart = max(iCentralBin - int(args.integralBinRange)/2, 0)
	iBinEnd = iBinStart + int(args.integralBinRange)
	if iBinEnd > nOtherMassBins:
		if args.integralBinRange > nOtherMassBins:
			pyRootPwa.utils.printErr("Too large range")
			sys.exit(1)
		iBinEnd = nOtherMassBins
		iBinStart = iBinEnd - args.integralBinRange
	selectedBins = [ otherMassBins[i] for i in range(iBinStart, iBinEnd) ]
	selectedBinIndices = [ fileManager.binList.index(b) for b in selectedBins ]
	jCentralBin = selectedBins.index(centralBin)
	pyRootPwa.utils.printInfo("Including the following bins:")
	for i in selectedBinIndices:
		pyRootPwa.utils.printInfo("\t"+str(fileManager.binList[i]))

	pyRootPwa.utils.printInfo("The central bin is:" + str(selectedBins[jCentralBin]))

	return selectedBinIndices, jCentralBin


def main():
	parser = argparse.ArgumentParser(
	                                 description="pwa pyPartialWaveFit fit executable"
	                                )

	parser.add_argument("outputFileName", type=str, metavar="fileName", help="path to output file")
	parser.add_argument("-c", type=str, metavar="configFileName", dest="configFileName", default="./rootpwa.config", help="path to config file (default: './rootpwa.config')")
	parser.add_argument("-b", type=int, metavar="#", dest="integralBin", default=0, help="integral bin id of fit (default: 0)")
	parser.add_argument("-B", type=int, metavar="#", dest="integralBinRange", default=7, help="Number of integral bins around the given interbral bin to fit (default: 7)")
	parser.add_argument("-s", type=int, metavar="#", dest="seed", default=None, help="random seed (default: draw from system)")
	parser.add_argument("-N", type=int, metavar="#", dest="nAttempts", default=1, help="number of fit attempts to perform")
	parser.add_argument("-w", type=str, metavar="path", dest="waveListFileName", default="", help="path to wavelist file (default: none)")
	parser.add_argument("-r", type=int, metavar="#", dest="rank", default=1, help="rank of spin density matrix (default: 1)")
	parser.add_argument("-H", "--checkHessian", help="check analytical Hessian eigenvalues (default: false)", action="store_true")
	parser.add_argument("-z", "--saveSpace", help="save space by not saving integral and covariance matrices (default: false)", action="store_true")
	parser.add_argument("--saveIntegrals", help="save integrals even if 'saveSpace' is defined -> do not calculate hessian but store integrals (default: false)", action="store_true")
	parser.add_argument("--saveAll", action="store_true",
	                    help="saving integral and covariance matrices of all fit attempts, not only of the best and best converged one (default: false)")
	parser.add_argument("-v", "--verbose", help="verbose; print debug output (default: false)", action="store_true")
	parser.add_argument("-L", "--likelihood", metavar="classname", default=None,
	                    help="Name of the likelihood class to use. Classes are: " + ", ".join(pyRootPwa.pyPartialWaveFit.getLikelihoodClassNames()))
	parser.add_argument("--likelihoodParameters", metavar="parameter-string", default=None, help="Parameter string given to the likelihood.setParameters(<parameter-string>) function")
	parser.add_argument("--likelihoodModule", metavar="path-to-likelihood-model", default=None, help="Implement the likelihood class not from ROOTPWA but from the given module-file")
	parser.add_argument("--dataset", action='append', dest='datasets', default=None, help="Define data-set to fit via data-set label.")
	args = parser.parse_args()

	clsModel = pyRootPwa.pyPartialWaveFit.ModelConnected
	clsLikelihood = pyRootPwa.pyPartialWaveFit.Likelihood
	clsParameterMapping = pyRootPwa.pyPartialWaveFit.ParameterMappingRpwa
	clsFitter = pyRootPwa.pyPartialWaveFit.NLoptFitter

	likelihoodModule = pyRootPwa.pyPartialWaveFit
	if args.likelihoodModule is not None:
		if args.likelihood is None:
			pyRootPwa.utils.printErr("Cannot use --'likelihoodModule' without '--likelihood'")
			sys.exit(1)
		if not os.path.exists(args.likelihoodModule):
			pyRootPwa.utils.printErr("Cannot find likelihood module '{0}'".format(args.likelihoodModule))
			sys.exit(1)
		likelihoodModule = imp.load_source('likelihoodModule', args.likelihoodModule)

	if args.likelihood is not None:
		exec("clsLikelihood = likelihoodModule.{l}".format(l=args.likelihood))

	pyRootPwa.utils.printInfo("Using likelihood '{0}' from module '{1}'.".format(clsLikelihood.__name__, likelihoodModule.__file__))

	binIndices, jCentralBin = buildBinRange(args)

	model = clsModel(clsLikelihood, clsParameterMapping)
	model.initModelInBins(args.configFileName, binIndices, args.waveListFileName, args.rank, args.rank, args.datasets)

	if args.likelihoodParameters is not None:
		exec("model.likelihood.setParameters({p})".format(p=args.likelihoodParameters))

	checkLevel = 0 if args.saveSpace else 1
	if args.checkHessian:
		checkLevel = 2

	storageLevel = 0 if args.saveSpace else 1
	if args.saveAll:
		storageLevel = 2
	fitter = clsFitter(
	                model,
	                checkLevel,
	                storageLevel,
	                pyRootPwa.pyPartialWaveFit.StartParameterGeneratorUniform(model, args.seed)
	               )

	fitResults = fitter.fit(args.nAttempts)
	if not fitResults:
		pyRootPwa.utils.printErr("didn't get valid fit result(s). Aborting...")
		sys.exit(1)

	fitResultsCentralBin = []
	for result in fitResults:
		result = dict(result)
		result['parameters'] = model.parameterMapping.paraFitterOfBin(result['parameters'], jCentralBin)
		fitResultsCentralBin.append(result)

	if args.saveIntegrals and storageLevel < 1:
		storageLevel = 1
	pyRootPwa.pyPartialWaveFit.writeResultsRpwa(fitter.model.models[jCentralBin], fitResultsCentralBin, args.outputFileName, storageLevel)

	sys.exit(0)


if __name__ == "__main__":
	main()
