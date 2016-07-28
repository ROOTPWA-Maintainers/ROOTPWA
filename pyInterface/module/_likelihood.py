import pyRootPwa.utils
ROOT = pyRootPwa.utils.ROOT

def initLikelihood(waveDescThres,
                   massBinCenter,
                   eventAndAmpFileDict,
                   normIntegralFileName,
                   accIntegralFileName,
                   multiBin,
                   accEventsOverride = 0,
                   cauchy = False,
                   cauchyWidth = 0.5,
                   rank = 1,
                   verbose = False
                  ):
	likelihood = pyRootPwa.core.pwaLikelihood()
	likelihood.useNormalizedAmps(True)
	if not verbose:
		likelihood.setQuiet()
	if cauchy:
		likelihood.setPriorType(pyRootPwa.core.pwaLikelihood.HALF_CAUCHY)
		likelihood.setCauchyWidth(cauchyWidth)
	if (not likelihood.init(waveDescThres,
	                        rank,
	                        massBinCenter)):
		pyRootPwa.utils.printErr("could not initialize likelihood. Aborting...")
		return None

	normIntFile = ROOT.TFile.Open(normIntegralFileName, "READ")
	normIntMeta = pyRootPwa.core.ampIntegralMatrixMetadata.readIntegralFile(normIntFile)
	normIntMatrix = normIntMeta.getAmpIntegralMatrix()
	if not likelihood.addNormIntegral(normIntMatrix):
		pyRootPwa.utils.printErr("could not add normalization integral. Aborting...")
		return None
	normIntFile.Close()
	accIntFile = ROOT.TFile.Open(accIntegralFileName, "READ")
	accIntMeta = pyRootPwa.core.ampIntegralMatrixMetadata.readIntegralFile(accIntFile)
	accIntMatrix = accIntMeta.getAmpIntegralMatrix()
	if not likelihood.addAccIntegral(accIntMatrix, accEventsOverride):
		pyRootPwa.utils.printErr("could not add acceptance integral. Aborting...")
		return None
	accIntFile.Close()

	eventMetas = []
	for eventFileName in eventAndAmpFileDict.keys():
		eventFile, eventMeta = pyRootPwa.utils.openEventFile(eventFileName)
		if not eventFile or not eventMeta:
			pyRootPwa.utils.printErr("could not open event file '" + eventFileName + "'. Aborting...")
			return None
		eventMetas.append(eventMeta)
	likelihood.setOnTheFlyBinning(multiBin.boundaries, eventMetas)
	for (waveName, _, _) in waveDescThres:
		ampMetas = []
		for eventFileName in eventAndAmpFileDict:
			ampFileName = eventAndAmpFileDict[eventFileName][waveName]
			ampFile = ROOT.TFile.Open(ampFileName, "READ")
			if not ampFile:
				pyRootPwa.utils.printErr("could not open amplitude file '" + ampFileName + "'.")
				return None
			meta = pyRootPwa.core.amplitudeMetadata.readAmplitudeFile(ampFile, waveName)
			if not meta:
				pyRootPwa.utils.printErr("could not get metadata for waveName '" + waveName + "'.")
				return None
			ampMetas.append(meta)
		if not likelihood.addAmplitude(ampMetas):
			pyRootPwa.utils.printErr("could not add amplitude '" + waveName + "'. Aborting...")
			return None
	if not likelihood.finishInit():
		pyRootPwa.utils.printErr("could not finish initialization of likelihood. Aborting...")
		return None

	return likelihood
