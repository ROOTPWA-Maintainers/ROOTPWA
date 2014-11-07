from pyRootPwa.utils import printInfo
from pyRootPwa.utils import printErr

def binningMapFromArgList(argList):
	binningMap = {}
	if argList:
		for bin in argList:
			try:
				splitUp = bin.split(";")
				if len(splitUp)==3:
					binningMap[splitUp[0]] = (float(splitUp[1]), float(splitUp[2]))
					printInfo("adding to binning map: " + splitUp[0] + " -> (" + splitUp[1] + "," + splitUp[2] + ")")
				else:
					printErr("did not get the right amount of semicolon seperated values for " + splitUp[0] + "-bin.")
					return False
			except ValueError:
				printErr("could not convert binning map boundaries of " + splitUp[0] + "-bin to float. Aborting...")
				return False
	return binningMap