#!/usr/bin/env python

import argparse
import math
import sys

import pyRootPwa
import pyRootPwa.core
ROOT = pyRootPwa.ROOT


if __name__ == "__main__":

	parser = argparse.ArgumentParser(
	                                 description="pwa fit executable"
	                                )

	parser.add_argument("inputFileName", type=str, metavar="inputFile", help="path to fit result")
	parser.add_argument("outputFileName", type=str, metavar="outputFile", help="path to output file")
	parser.add_argument("-c", type=str, metavar="configFileName", dest="configFileName", default="./rootpwa.config", help="path to config file (default: './rootpwa.config')")
	parser.add_argument("-b", type=int, metavar="#", dest="integralBin", default=0, help="integral bin id of fit (default: 0)")
	parser.add_argument("-s", action="store_true", dest="ranWithHesse", default=False, help="indicates that the fit result contains an analytically calculated covariance matrix")
	parser.add_argument("-C", "--cauchyPriors", help="use half-Cauchy priors (default: false)", action="store_true")
	parser.add_argument("-P", "--cauchyPriorWidth", type=float, metavar ="WIDTH", default=0.5, help="width of half-Cauchy prior (default: 0.5)")
	parser.add_argument("-A", type=int, metavar="#", dest="accEventsOverride", default=0,
	                    help="number of input events to normalize acceptance to (default: use number of events from normalization integral file)")
	parser.add_argument("--noAcceptance", help="do not take acceptance into account (default: false)", action="store_true")
	parser.add_argument("-v", "--verbose", help="verbose; print debug output (default: false)", action="store_true")
	args = parser.parse_args()

	config = pyRootPwa.rootPwaConfig()
	if not config.initialize(args.configFileName):
		pyRootPwa.utils.printErr("loading config file '" + args.configFileName + "' failed. Aborting...")
		sys.exit(1)
	pyRootPwa.core.particleDataTable.readFile(config.pdgFileName)
	fileManager = pyRootPwa.loadFileManager(config.fileManagerPath)
	if not fileManager:
		pyRootPwa.utils.printErr("loading the file manager failed. Aborting...")
		sys.exit(1)

	if args.integralBin < 0:
		pyRootPwa.utils.printErr("bin < 0 (" + str(args.integralBin) + "). Aborting...")
		sys.exit(1)
	elif args.integralBin >= len(fileManager.binList):
		pyRootPwa.utils.printErr("bin out of range (" + str(args.integralBin) + ">=" + str(len(fileManager.binList)) + "). Aborting...")
		sys.exit(1)
	multiBin = fileManager.binList[args.integralBin]
	eventAndAmpFileDict = fileManager.getEventAndAmplitudeFilePathsInBin(multiBin, pyRootPwa.core.eventMetadata.REAL)
	if not eventAndAmpFileDict:
		pyRootPwa.utils.printErr("could not retrieve valid amplitude file list. Aborting...")
		sys.exit(1)

	psIntegralPath  = fileManager.getIntegralFilePath(multiBin, pyRootPwa.core.eventMetadata.GENERATED)
	accIntegralPath = psIntegralPath
	if not args.noAcceptance:
		accIntegralPath = fileManager.getIntegralFilePath(multiBin, pyRootPwa.core.eventMetadata.ACCEPTED)
	elif args.accEventsOverride != 0:
		# for a fit without acceptance corrections the number of events
		# the acceptance matrix is normalized to needs to be equal to
		# the number of events in the normalization matrix
		intFile = pyRootPwa.ROOT.TFile.Open(psIntegralPath, "READ")
		intMeta = pyRootPwa.core.ampIntegralMatrixMetadata.readIntegralFile(intFile)
		intMatrix = intMeta.getAmpIntegralMatrix()
		if args.accEventsOverride != intMatrix.nmbEvents():
			pyRootPwa.utils.printErr("incorrect number of events for normalization of integral matrix for "
			                         "a fit without acceptance (got: {:d}, expected: {:d}). Aborting...".format(args.accEventsOverride, intMatrix.nmbEvents()))
			sys.exit(1)

	result = pyRootPwa.utils.getFitResultFromFile(fitResultFileName = args.inputFileName,
	                                              fitResultTreeName = config.fitResultTreeName,
	                                              fitResultBranchName = config.fitResultBranchName)
	if not result:
		pyRootPwa.utils.printErr("could not get fit result from file '" + args.inputFileName + "'. Aborting...")
		sys.exit(1)

	waveDescThres = pyRootPwa.utils.getWaveDescThresFromFitResult(result, fileManager.getWaveDescriptions())
	if not waveDescThres:
		pyRootPwa.utils.printErr("error while getting wave names, descriptions and thresholds. Aborting...")
		sys.exit(1)

	likelihood = pyRootPwa.initLikelihood(waveDescThres = waveDescThres,
	                                      massBinCenter = result.massBinCenter(),
	                                      eventAndAmpFileDict = eventAndAmpFileDict,
	                                      normIntegralFileName = psIntegralPath,
	                                      accIntegralFileName = accIntegralPath,
	                                      multiBin = multiBin,
	                                      accEventsOverride = args.accEventsOverride,
	                                      cauchy = args.cauchyPriors,
	                                      cauchyWidth = args.cauchyPriorWidth,
	                                      rank = result.rank(),
	                                      verbose = args.verbose)
	if not likelihood:
		pyRootPwa.utils.printErr("error while initializing likelihood. Aborting...")
		sys.exit(1)

	minimum = ROOT.TVectorD(likelihood.nmbPars())
	for i, parameter in enumerate(likelihood.parameters()):
		minimum[i] = result.fitParameter(parameter.parName())

	minimumLikelihood = likelihood.DoEval( [ minimum[i] for i in xrange(minimum.GetNrows()) ] )
	pyRootPwa.utils.printInfo("likelihood at minimum is {: .15e}.".format(minimumLikelihood))

	covMatrixMinuit = result.fitParCovMatrix()
	eigenVectorsMinuit = likelihood.HessianEigenVectors(covMatrixMinuit)

	if not args.ranWithHesse:
		covMatrixAna = likelihood.CovarianceMatrix( [ minimum[i] for i in xrange(minimum.GetNrows()) ] )
		eigenVectorsAna = likelihood.HessianEigenVectors(covMatrixAna)

		if len(eigenVectorsMinuit) != len(eigenVectorsAna):
			pyRootPwa.utils.printErr("different number of Eigenvalues for covariance matrix of Minuit vs. analytic covariance matrix.")
			sys.exit(1)

	outputFile = ROOT.TFile.Open(args.outputFileName, "NEW")
	if (not outputFile) or outputFile.IsZombie():
		pyRootPwa.utils.printErr("cannot open output file '" + args.outputFileName + "'. Aborting...")
		sys.exit(1)

	for par in xrange(len(eigenVectorsMinuit)):
		# The idea here is to scan the logLikelihood function from the
		# minimum in direction of the Eigenvectors.
		#
		# Around the minimum p the logLikelihood can be expanded as a
		# Taylor series:
		#
		# -logL(p+dp) = -logL(p) + dp^T * grad(-logL(x))|x=p + 1/2 * dp^T * hesse(-logL(x))|x=p * dp
		#
		# The second summand should be equal to zero at the minimum.
		# It is assumed that the Hessian matrix of -logL is the
		# inverse of the covariance matrix
		#
		# hesse(-logL(x))|x=p = C^-1
		#
		# Eigenvectors of the covariance matrix then are also
		# Eigenvectors of the Hessian matrix, albeit with inverse
		# Eigenvalues. If e is an Eigenvector of the covariance matrix
		# with Eigenvalue l
		#
		# C * e = l * e
		#
		# then
		#
		# hesse(-logL(x))|x=p * e = 1/l * e
		#
		# We now write dp = r * sqrt(l) * e where r in the following
		# is varied between -1 and +1. For a logLikelihood behaving
		# nicely we then should find the parabola
		#
		# -logL(p+dp) = -logL(p) + 1/2 * r^2
		#
		# This ideal parabola is drawn for a x between -sqrt(l) and
		# +sqrt(l) using x = r * sqrt(l).

		graphLikeli = ROOT.TGraph()
		for p in xrange(-50, 51):
			ratio = (p/50.0) * math.sqrt(eigenVectorsMinuit[par][1])
			# TODO: Fix this once root get's their operators in working order
			step = ROOT.TVectorD(eigenVectorsMinuit[par][0])
			step *= ratio
			pars = minimum + step
			likeli = likelihood.DoEval( [ pars[i] for i in xrange(pars.GetNrows()) ] ) - minimumLikelihood
			graphLikeli.SetPoint(p + 50, ratio, likeli)

		lowerLimit = -1.2 * math.sqrt(eigenVectorsMinuit[par][1])
		upperLimit =  1.2 * math.sqrt(eigenVectorsMinuit[par][1])

		# the magnitude of the eigenvector is 1. otherwise an
		# additional factor '|eigenvector|^2' would be required
		parabolaMinuit = ROOT.TF1("minuitParabola", "0.5 / {:.15e} * x*x".format(eigenVectorsMinuit[par][1]), lowerLimit, upperLimit)
		if not args.ranWithHesse:
			parabolaAna = ROOT.TF1("analyticParabola", "0.5 / {:.15e} * x*x".format(eigenVectorsAna[par][1]), lowerLimit, upperLimit)

		canvas = ROOT.TCanvas("eigenvectorSlice{:d}".format(par))
		canvas.cd()
		graphLikeli.Draw("A*")
		parabolaMinuit.Draw("Lsame")
		parabolaMinuit.SetLineColor(ROOT.kBlue)
		if not args.ranWithHesse:
			parabolaAna.Draw("Lsame")
			parabolaAna.SetLineColor(ROOT.kRed)
		canvas.Write()

	outputFile.Close()
	pyRootPwa.utils.printSucc("slices successfully written to file '{}'.".format(args.outputFileName))
	pyRootPwa.utils.printInfo("setting Minuit line color to blue, analytical solution to red.")
