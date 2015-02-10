#!/usr/bin/env python

import argparse
import sys

import pyRootPwa
import pyRootPwa.utils

if __name__ == "__main__":

	# parse command line arguments
	parser = argparse.ArgumentParser(
	                                 description="calculates decay amplitudes "
	                                             "for given wave for events in "
	                                             "input data files and "
	                                             "writes amplitudes to file"
	                                )

	parser.add_argument("-c", type=str, metavar="config-file", default="rootpwa.config", dest="configFileName", help="path to config file (default: ./rootpwa.config)")
	parser.add_argument("-n", type=int, metavar="#", default=-1, dest="maxNmbEvents",  help="maximum number of events to read (default: all)")
	parser.add_argument("-b", type=str, metavar="massBin(s)", default="all", dest="massBins", help="mass bins to be calculated (default: all)")
	parser.add_argument("-j", type=int, metavar=("jobs"), default=1, dest="nJobs", help="number of jobs (default: 1)")
	parser.add_argument("-f", "--no-progress-bar", action="store_true", dest="noProgressBar", help="disable progress bars (decreases computing time)")
	parser.add_argument("-k", "--keyfiles", type=str, metavar="keyfiles", dest="keyfiles", nargs="*", help="keyfiles to calculate amplitude for (overrides settings from the config file)")
#	parser.add_argument("-v", action="store_true", dest="debug", help="verbose; print debug output (default: false)")
	parser.add_argument("--profiler", type=str, metavar="profiler-output", dest="proFile", help="path to profiler output file")

	arguments = parser.parse_args()

	if arguments.proFile is None:

		args = {"configFileName": arguments.configFileName,
		        "massBins": arguments.massBins,
		        "nJobs": arguments.nJobs,
		        "progressBar": (not arguments.noProgressBar),
		        "maxNmbEvents": arguments.maxNmbEvents}

		if arguments.keyfiles is not None:
			args["keyfiles"] = arguments.keyfiles

		pyRootPwa.calcAmplitudes(**args)

	else:

		import profile
		import time
		print("Profiler activated.")
		print("WARNING: runtime will increase significantly, use this only for debug purposes!")
		print
		print("To view the generated profile, do:")
		print("$> python")
		print(">>> import pstats")
		print(">>> p = pstats.Stats('<profile file>')")
		print(">>> p.strip_dirs().sort_stats('time').print_stats()")
		print
		print("Sleeping for 10 seconds...")
		try:
			time.sleep(10)
		except KeyboardInterrupt:
			print('\nRecieved KeyboardInterrupt. Aborting...')
			sys.exit(1)
		print("Starting now...")
		print
		print
		profile.run('pyRootPwa.calcAmplitudes(configFileName=arguments.configFileName, massBins=arguments.massBins, nJobs=arguments.nJobs, proFile=arguments.proFile)', arguments.proFile)
