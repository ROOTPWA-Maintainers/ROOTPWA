
import inspect
import sys

config = None

__terminalColorStrings = {}
__terminalColorStrings['normal']     = "\033[0m"
__terminalColorStrings['bold']       = "\033[1m"
__terminalColorStrings['underline']  = "\033[4m"
__terminalColorStrings['blink']      = "\033[5m"
__terminalColorStrings['inverse']    = "\033[7m"
__terminalColorStrings['fgBlack']    = "\033[30m"
__terminalColorStrings['fgRed']      = "\033[31m"
__terminalColorStrings['fgGreen']    = "\033[32m"
__terminalColorStrings['fgYellow']   = "\033[33m"
__terminalColorStrings['fgBlue']     = "\033[34m"
__terminalColorStrings['fgMangenta'] = "\033[35m"
__terminalColorStrings['fgCyan']     = "\033[36m"
__terminalColorStrings['fgWhite']    = "\033[37m"
__terminalColorStrings['bgBlack']    = "\033[40m"
__terminalColorStrings['bgRed']      = "\033[41m"
__terminalColorStrings['bgGreen']    = "\033[42m"
__terminalColorStrings['bgYellow']   = "\033[43m"
__terminalColorStrings['bgBlue']     = "\033[44m"
__terminalColorStrings['bgMangenta'] = "\033[45m"
__terminalColorStrings['bgCyan']     = "\033[46m"
__terminalColorStrings['bgWhite']    = "\033[47m"

def __printFormatted(msg, level):
	frame = inspect.currentframe().f_back.f_back
	if frame is None:
		printErr("This method cannot be called directly.")
	(filename, lineno, function, code_contex, index) = inspect.getframeinfo(frame)
	if function == "<module>":
		function = "__main__"
	string = ""
	if level == "err":
		if sys.stderr.isatty(): string += __terminalColorStrings['fgRed']
		string += "!!! "
		string += function + " [" + filename + ":" + str(lineno) + "]: "
		string += "error: "
		if sys.stderr.isatty(): string += __terminalColorStrings['normal']
		string += msg
	elif level == "warn":
		if sys.stderr.isatty(): string += __terminalColorStrings['fgYellow']
		string += "??? "
		string += function + " [" + filename + ":" + str(lineno) + "]: "
		string += "warning: "
		if sys.stderr.isatty(): string += __terminalColorStrings['normal']
		string += msg
	elif level == "suc":
		if sys.stdout.isatty(): string += __terminalColorStrings['fgGreen']
		string += "*** "
		string += function + ": success: "
		if sys.stdout.isatty(): string += __terminalColorStrings['normal']
		string += msg
	elif level == "info":
		if sys.stdout.isatty(): string += __terminalColorStrings['bold']
		string += ">>> "
		string += function + ": info: "
		if sys.stdout.isatty(): string += __terminalColorStrings['normal']
		string += msg
	elif level == "debug":
		if sys.stdout.isatty(): string += __terminalColorStrings['fgMangenta']
		string += "+++ "
		string += function + ": debug: "
		string += __terminalColorStrings['normal']
		string += msg
	else:
		printErr("Invalid level string.")
		raise Exception()
	if level == "err" or level == "warn":
		sys.stderr.write(string + "\n")
	else:
		sys.stdout.write(string + "\n")


def printErr(msg): __printFormatted(msg, "err")
def printWarn(msg): __printFormatted(msg, "warn")
def printSucc(msg): __printFormatted(msg, "suc")
def printInfo(msg): __printFormatted(msg, "info")
def printDebug(msg): __printFormatted(msg, "debug")

