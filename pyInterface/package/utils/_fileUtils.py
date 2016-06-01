
import pyRootPwa.core

import _root
ROOT = _root.ROOT

def openEventFile(eventFileName):
	eventFile = ROOT.TFile.Open(eventFileName, "READ")
	if not eventFile:
		pyRootPwa.utils.printWarn("could not open event file '" + eventFileName + "'.")
		return (None, None)
	eventMeta = pyRootPwa.core.eventMetadata.readEventFile(eventFile)
	if not eventMeta:
		pyRootPwa.utils.printWarn("could not read metadata from event file '" + eventFileName + "'.")
		return (None, None)
	return (eventFile, eventMeta)
