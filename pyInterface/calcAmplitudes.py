#!/usr/bin/python2.7

import argparse

import pyRootPwa.bin

if __name__ == "__main__":

	# parse command line arguments
	parser = argparse.ArgumentParser(
	                                 description="calculates decay amplitudes "
	                                             "for given wave for events in "
	                                             "input data files and "
	                                             "writes amplitudes to file",
	                                 formatter_class=argparse.RawTextHelpFormatter
	                                )

	parser.add_argument("-c", type=str, metavar="file", default="rootpwa.config", dest="configFileName", help="path ot config file (default: ./rootpwa.config)")
	parser.add_argument("-n", type=int, metavar="#", default=-1, dest="maxNmbEvents",  help="maximum number of events to read (default: all)")
	parser.add_argument("-mb", type=str, metavar="massBin(s)", default="all", dest="massBins", help="mass bins to be calculated (default: all)")
	parser.add_argument("-j", type=int, metavar=("jobs"), default=1, dest="nJobs", help="EXPERIMENTAL: number of jobs used (default: 1)")
	parser.add_argument("-v", action="store_true", dest="debug", help="verbose; print debug output (default: false)")

	arguments = parser.parse_args()

	pyRootPwa.bin.calcAmplitudes(configFileName=arguments.configFileName,
                                 massBins=arguments.massBins,
                                 nJobs=arguments.nJobs)

