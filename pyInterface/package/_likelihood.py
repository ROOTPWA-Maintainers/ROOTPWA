import pyRootPwa
ROOT = pyRootPwa.ROOT

def initLikelihood(waveDescThres,
                   massBinCenter,
                   ampFileList,
                   normIntegralFileName,
                   accIntegralFileName,
                   accEventsOverride = 0,
                   cauchy = False,
                   cauchyWidth = 0.5,
                   rank = 1,
                   verbose = False
                  ):
	likelihood = pyRootPwa.core.pwaLikelihood()
	likelihood.useNormalizedAmps(True)
	if (not verbose):
		likelihood.setQuiet()
	if cauchy:
		likelihood.setPriorType(pyRootPwa.core.HALF_CAUCHY)
		likelihood.setCauchyWidth(cauchyWidth)
	if (not likelihood.init(waveDescThres,
	                        rank,
	                        massBinCenter)):
		printErr("could not initialize likelihood. Aborting...")
		return None

	normIntFile = ROOT.TFile.Open(normIntegralFileName, "READ")
	if len(normIntFile.GetListOfKeys()) != 1:
		pyRootPwa.utils.printWarn("'" + normIntegralFileName + "' does not contain exactly one TKey.")
		return None
	normIntMatrix = normIntFile.Get(pyRootPwa.core.ampIntegralMatrix.integralObjectName)
	if (not likelihood.addNormIntegral(normIntMatrix)):
		pyRootPwa.utils.printErr("could not add normalization integral. Aborting...")
		return None
	normIntFile.Close()
	accIntFile = ROOT.TFile.Open(accIntegralFileName, "READ")
	if len(accIntFile.GetListOfKeys()) != 1:
		pyRootPwa.utils.printWarn("'" + accIntegralFileName + "' does not contain exactly one TKey.")
		return None
	accIntMatrix = accIntFile.Get(pyRootPwa.core.ampIntegralMatrix.integralObjectName)
	if (not likelihood.addAccIntegral(accIntMatrix, accEventsOverride)):
		pyRootPwa.utils.printErr("could not add acceptance integral. Aborting...")
		return None
	accIntFile.Close()

	for wave in waveDescThres:
		waveName = wave[0]
		ampFileName = ampFileList[waveName]
		ampFile = ROOT.TFile.Open(ampFileName, "READ")
		if not ampFile:
			pyRootPwa.utils.printErr("could not open amplitude file '" + ampFileName + "'.")
			return None
		meta = pyRootPwa.core.amplitudeMetadata.readAmplitudeFile(ampFile, waveName)
		if not meta:
			pyRootPwa.utils.printErr("could not get metadata for waveName '" + waveName + "'.")
			return None
		if (not likelihood.addAmplitude(meta)):
			pyRootPwa.utils.printErr("could not add amplitude '" + waveName + "'. Aborting...")
			return None
	if (not likelihood.finishInit()):
		pyRootPwa.utils.printErr("could not finish initialization of likelihood. Aborting...")
		return None

	return likelihood
