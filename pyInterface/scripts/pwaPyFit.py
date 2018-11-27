#!/usr/bin/env python

import argparse
import sys
import os

import pyRootPwa
import pyRootPwa.pyPartialWaveFit


def main(clsModel, clsLikelihood, clsParameterMapping, clsFitter):
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
# 	parser.add_argument("--noAcceptance", help="do not take acceptance into account (default: false)", action="store_true")
# 	parser.add_argument("-C", "--cauchyPriors", help="use half-Cauchy priors (default: false)", action="store_true")
# 	parser.add_argument("-P", "--cauchyPriorWidth", type=float, metavar ="WIDTH", default=0.5, help="width of half-Cauchy prior (default: 0.5)")
	args = parser.parse_args()

	# disable multithreading by default
	os.environ['OPENBLAS_NUM_THREADS'] = 1

	model = clsModel(clsLikelihood, clsParameterMapping)

	model.initModelInBin(args.configFileName, args.integralBin, args.waveListFileName, args.rank, args.rank)

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
	                pyRootPwa.pyPartialWaveFit.StartParameterGeneratorRpwaUniform(model, args.seed)
	               )

	fitResults = fitter.fit(args.nAttempts)
	if not fitResults:
		pyRootPwa.utils.printErr("didn't get valid fit result(s). Aborting...")
		sys.exit(1)

	fitter.writeResultsRpwa(fitResults, args.outputFileName)

	sys.exit(0)


if __name__ == "__main__":
	main(pyRootPwa.pyPartialWaveFit.ModelRpwa,
	     pyRootPwa.pyPartialWaveFit.Likelihood,
	     pyRootPwa.pyPartialWaveFit.ParameterMappingRpwa,
	     pyRootPwa.pyPartialWaveFit.NLoptFitter)
