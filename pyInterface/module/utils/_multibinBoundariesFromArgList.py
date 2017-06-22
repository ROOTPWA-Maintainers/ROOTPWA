from _printingUtils import printInfo
from _printingUtils import printErr

def multibinBoundariesFromArgList(argList):
	multibinBoundaries = {}
	if argList:
		for binArg in argList:
			try:
				splitUp = binArg.split(";")
				if len(splitUp)==3:
					multibinBoundaries[splitUp[0]] = (float(splitUp[1]), float(splitUp[2]))
					printInfo("adding to multibin boundaries: " + splitUp[0] + " -> (" + splitUp[1] + "," + splitUp[2] + ")")
				else:
					printErr("did not get the right amount of semicolon seperated values for " + splitUp[0] + "-bin.")
					return False
			except ValueError:
				printErr("could not convert multibin boundaries of " + splitUp[0] + "-bin to float. Aborting...")
				return False
	return multibinBoundaries
