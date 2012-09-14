
import inspect
import multiprocessing
import sys

import pyRootPwa.utils

stdoutLock = multiprocessing.Lock()
stdoutisatty = None
stderrisatty = None

def printPrintingSummary():
	print("Number of errors printed:    " + str(pyRootPwa.utils.printErr.count))
	print("Number of warnings printed:  " + str(pyRootPwa.utils.printWarn.count))
	print("Number of infos printed:     " + str(pyRootPwa.utils.printInfo.count))
	print("Number of successes printed: " + str(pyRootPwa.utils.printSucc.count))
	print("Number of debugs printed:    " + str(pyRootPwa.utils.printDebug.count))

class printClass:

	count = 0

	_terminalColorStrings = {}
	_terminalColorStrings['normal']     = "\033[0m"
	_terminalColorStrings['bold']       = "\033[1m"
	_terminalColorStrings['underline']  = "\033[4m"
	_terminalColorStrings['blink']      = "\033[5m"
	_terminalColorStrings['inverse']    = "\033[7m"
	_terminalColorStrings['fgBlack']    = "\033[30m"
	_terminalColorStrings['fgRed']      = "\033[31m"
	_terminalColorStrings['fgGreen']    = "\033[32m"
	_terminalColorStrings['fgYellow']   = "\033[33m"
	_terminalColorStrings['fgBlue']     = "\033[34m"
	_terminalColorStrings['fgMangenta'] = "\033[35m"
	_terminalColorStrings['fgCyan']     = "\033[36m"
	_terminalColorStrings['fgWhite']    = "\033[37m"
	_terminalColorStrings['bgBlack']    = "\033[40m"
	_terminalColorStrings['bgRed']      = "\033[41m"
	_terminalColorStrings['bgGreen']    = "\033[42m"
	_terminalColorStrings['bgYellow']   = "\033[43m"
	_terminalColorStrings['bgBlue']     = "\033[44m"
	_terminalColorStrings['bgMangenta'] = "\033[45m"
	_terminalColorStrings['bgCyan']     = "\033[46m"
	_terminalColorStrings['bgWhite']    = "\033[47m"

	def __init__(self):
		self.count = 0

	def printFormatted(self, msg, level):
		frame = inspect.currentframe().f_back.f_back
		if frame is None:
			printErr("This method cannot be called directly.")
		(filename, lineno, function, code_contex, index) = inspect.getframeinfo(frame)
		if function == "<module>":
			function = "__main__"
		string = ""
		if level == "err":
			if pyRootPwa.utils.stderrisatty: string += self._terminalColorStrings['fgRed']
			string += "!!! "
			string += function + " [" + filename + ":" + str(lineno) + "]: "
			string += "error: "
			if pyRootPwa.utils.stderrisatty: string += self._terminalColorStrings['normal']
			string += msg
		elif level == "warn":
			if pyRootPwa.utils.stderrisatty: string += self._terminalColorStrings['fgYellow']
			string += "??? "
			string += function + " [" + filename + ":" + str(lineno) + "]: "
			string += "warning: "
			if pyRootPwa.utils.stderrisatty: string += self._terminalColorStrings['normal']
			string += msg
		elif level == "suc":
			if pyRootPwa.utils.stdoutisatty: string += self._terminalColorStrings['fgGreen']
			string += "*** "
			string += function + ": success: "
			if pyRootPwa.utils.stdoutisatty: string += self._terminalColorStrings['normal']
			string += msg
		elif level == "info":
			if pyRootPwa.utils.stdoutisatty: string += self._terminalColorStrings['bold']
			string += ">>> "
			string += function + ": info: "
			if pyRootPwa.utils.stdoutisatty: string += self._terminalColorStrings['normal']
			string += msg
		elif level == "debug":
			if pyRootPwa.utils.stdoutisatty: string += self._terminalColorStrings['fgMangenta']
			string += "+++ "
			string += function + ": debug: "
			string += self._terminalColorStrings['normal']
			string += msg
		else:
			printErr("Invalid level string.")
			raise Exception()
		if level == "err" or level == "warn":
			sys.stderr.write(string + "\n")
		else:
			sys.stdout.write(string + "\n")

class printErrClass(printClass):

	def __call__(self, msg):
		self.count += 1
		self.printFormatted(str(msg), "err")

class printWarnClass(printClass):

	def __call__(self, msg):
		self.count += 1
		self.printFormatted(str(msg), "warn")

class printSuccClass(printClass):

	def __call__(self, msg):
		self.count += 1
		self.printFormatted(str(msg), "suc")

class printInfoClass(printClass):

	def __call__(self, msg):
		self.count += 1
		self.printFormatted(str(msg), "info")

class printDebugClass(printClass):

	def __call__(self, msg):
		self.count += 1
		self.printFormatted(str(msg), "debug")

printErr = printErrClass()
printWarn = printWarnClass()
printSucc = printSuccClass()
printInfo = printInfoClass()
printDebug = printDebugClass()

