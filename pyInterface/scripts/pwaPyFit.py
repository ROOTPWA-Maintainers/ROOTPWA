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


def main():
	parser = argparse.ArgumentParser(
	                                 description="pwa pyPartialWaveFit fit executable"
	                                )

	parser.add_argument("outputFileName", type=str, metavar="fileName", help="path to output file")
	parser.add_argument("-c", type=str, metavar="configFileName", dest="configFileName", default="./rootpwa.config", help="path to config file (default: './rootpwa.config')")
	parser.add_argument("-b", type=int, metavar="#", dest="integralBin", default=0, help="integral bin id of fit (default: 0)")
	parser.add_argument("-s", type=int, metavar="#", dest="seed", default=None, help="random seed (default: draw from system)")
	parser.add_argument("-N", type=int, metavar="#", dest="nAttempts", default=1, help="number of fit attempts to perform")
	parser.add_argument("-w", type=str, metavar="path", dest="waveListFileName", default="", help="path to wavelist file (default: none)")
	parser.add_argument("-r", type=int, metavar="#", dest="rank", default=1, help="rank of spin density matrix (default: 1)")
	parser.add_argument("-H", "--checkHessian", help="check analytical Hessian eigenvalues (default: false)", action="store_true")
	parser.add_argument("-z", "--saveSpace", help="save space by not saving integral and covariance matrices (default: false)", action="store_true")
	parser.add_argument("--saveAll", action="store_true",
	                    help="saving integral and covariance matrices of all fit attempts, not only of the best and best converged one (default: false)")
	parser.add_argument("-v", "--verbose", help="verbose; print debug output (default: false)", action="store_true")
	parser.add_argument("-L", "--likelihood", metavar="classname", default=None,
	                    help="Name of the likelihood class to use. Classes are: " + ", ".join(pyRootPwa.pyPartialWaveFit.getLikelihoodClassNames()))
	parser.add_argument("--likelihoodParameters", metavar="parameter-string", default=None, help="Parameter string given to the likelihood.setParameters(<parameter-string>) function")
	parser.add_argument("--likelihoodModule", metavar="path-to-likelihood-model", default=None, help="Implement the likelihood class not from ROOTPWA but from the given module-file")
	parser.add_argument("--dataset", action='append', dest='datasets', default=None, help="Define data-set to fit via data-set label.")
	parser.add_argument("--drop-flatwave", action='store_true', dest='dropFlatwave', default=False, help="Do not include incoherent flat wave into the fit.")
	parser.add_argument("--noAcceptance", help="do not take acceptance into account (default: false)", action="store_true")
	parser.add_argument("--nEvents", type=int, metavar="#", dest="nEvents", default=None,
	                    help="Do not fit the complete data set but include only the given number of events (choosen randomly)")
	parser.add_argument("--ratioOfEvents", type=float, metavar="#", dest="ratioOfEvents", default=None,
	                    help="Do not fit the complete data set but include only the given ratio w.r.t. the total number of events (choosen randomly)")
	args = parser.parse_args()

	clsModel = pyRootPwa.pyPartialWaveFit.ModelRpwa
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

	model = clsModel(clsLikelihood, clsParameterMapping)
	model.initModelInBin(args.configFileName, args.integralBin, args.waveListFileName, args.rank, args.rank, args.datasets,
	                     addFlatWave=not args.dropFlatwave, noAcceptance=args.noAcceptance)

	if args.likelihoodParameters is not None:
		exec("model.likelihood.setParameters({p})".format(p=args.likelihoodParameters))

	startParameterGenerator = pyRootPwa.pyPartialWaveFit.StartParameterGeneratorRpwaUniform(model, args.seed)

	if args.ratioOfEvents is not None:
		if args.nEvents is not None:
			pyRootPwa.utils.printErr("Cannot use `--ratioOfEvents` and `--nEvents` at the same time!")
			sys.exit(1)
		args.nEvents = int(round(args.ratioOfEvents*model.likelihood.nmbEvents.sum()))

	if args.nEvents is not None:
		model.likelihood.useNEvents(args.nEvents, startParameterGenerator.generator)
		pyRootPwa.utils.printInfo("Use {0} events in the fit.".format(model.likelihood.nmbEvents.sum()))

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
	                startParameterGenerator
	               )

	fitResults = fitter.fit(args.nAttempts, verbosity=2 if args.verbose else 1)
	if not fitResults:
		pyRootPwa.utils.printErr("didn't get valid fit result(s). Aborting...")
		sys.exit(1)

	pyRootPwa.pyPartialWaveFit.writeResultsRpwa(fitter.model, fitResults, args.outputFileName, storageLevel)

	sys.exit(0)


if __name__ == "__main__":
	main()
