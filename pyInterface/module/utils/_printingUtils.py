from __future__ import print_function
import inspect as _inspect
import multiprocessing as _multiprocessing
import sys as _sys

# making the output nice for multi-threading
stdoutisatty = _sys.stdout.isatty()
stderrisatty = _sys.stderr.isatty()

def printPrintingSummary(printingCounterVar):
	print('')
	print("Number of errors printed:    " + str(printingCounterVar[0]))
	print("Number of warnings printed:  " + str(printingCounterVar[1]))
	print("Number of infos printed:     " + str(printingCounterVar[3]))
	print("Number of successes printed: " + str(printingCounterVar[2]))
	print("Number of debugs printed:    " + str(printingCounterVar[4]))

class _printClass(object):

	def __init__(self):
		pass

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

	def printFormatted(self, msg, level):
		frame = _inspect.currentframe().f_back.f_back
		if frame is None:
			printErr("This method cannot be called directly.")
		(filename, lineno, function, _, _) = _inspect.getframeinfo(frame)
		if function == "<module>":
			function = "__main__"
		string = ""
		if level == "err":
			if stderrisatty: string += self._terminalColorStrings['fgRed']
			string += "!!! "
			string += function + " [" + filename + ":" + str(lineno) + "]: "
			string += "error: "
			if stderrisatty: string += self._terminalColorStrings['normal']
			string += msg
		elif level == "warn":
			if stderrisatty: string += self._terminalColorStrings['fgYellow']
			string += "??? "
			string += function + " [" + filename + ":" + str(lineno) + "]: "
			string += "warning: "
			if stderrisatty: string += self._terminalColorStrings['normal']
			string += msg
		elif level == "suc":
			if stdoutisatty: string += self._terminalColorStrings['fgGreen']
			string += "*** "
			string += function + ": success: "
			if stdoutisatty: string += self._terminalColorStrings['normal']
			string += msg
		elif level == "info":
			if stdoutisatty: string += self._terminalColorStrings['bold']
			string += ">>> "
			string += function + ": info: "
			if stdoutisatty: string += self._terminalColorStrings['normal']
			string += msg
		elif level == "debug":
			if stdoutisatty: string += self._terminalColorStrings['fgMangenta']
			string += "+++ "
			string += function + ": debug: "
			if stdoutisatty: string += self._terminalColorStrings['normal']
			string += msg
		else:
			printErr("Invalid level string.")
			raise Exception()
		if level == "err" or level == "warn":
			_sys.stderr.write(string + "\n")
		else:
			_sys.stdout.write(string + "\n")

class printErrClass(_printClass):

	def __init__(self, counter):
		_printClass.__init__(self)
		self.counter = counter

	def __call__(self, msg):
		self.counter[0] += 1
		self.printFormatted(str(msg), "err")

class printWarnClass(_printClass):

	def __init__(self, counter):
		_printClass.__init__(self)
		self.counter = counter

	def __call__(self, msg):
		self.counter[1] += 1
		self.printFormatted(str(msg), "warn")

class printSuccClass(_printClass):

	def __init__(self, counter):
		_printClass.__init__(self)
		self.counter = counter

	def __call__(self, msg):
		self.counter[2] += 1
		self.printFormatted(str(msg), "suc")

class printInfoClass(_printClass):

	def __init__(self, counter):
		_printClass.__init__(self)
		self.counter = counter

	def __call__(self, msg):
		self.counter[3] += 1
		self.printFormatted(str(msg), "info")

class printDebugClass(_printClass):

	def __init__(self, counter):
		_printClass.__init__(self)
		self.counter = counter

	def __call__(self, msg):
		self.counter[4] += 1
		self.printFormatted(str(msg), "debug")

try:
	printingCounter
except NameError:
	printingCounter = _multiprocessing.Array('i', [0]*5)
	printErr = printErrClass(printingCounter)
	printWarn = printWarnClass(printingCounter)
	printSucc = printSuccClass(printingCounter)
	printInfo = printInfoClass(printingCounter)
	printDebug = printDebugClass(printingCounter)
