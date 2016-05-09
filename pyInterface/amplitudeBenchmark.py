#!/usr/bin/env python
# pylint: disable=line-too-long
import sys
import urllib
import os
import glob
import argparse
import numpy

import compareAmplitudeFiles

import pyRootPwa
import pyRootPwa.core
ROOT = pyRootPwa.ROOT

def treeToArray(tree):
	arr = numpy.zeros(tree.GetEntries())
	i = 0
	for evt in tree:
		arr[i] = evt.diff
		i += 1
	arr = numpy.sort(arr)
	return arr

def createDirectory(path):
	if not os.path.exists(path):
		os.mkdir(path)
		pyRootPwa.utils.printInfo("created directory '" + path + "'.")

def createDirectories(path):
	createDirectory(path)
	createDirectory(path + "/keyfiles")
	createDirectory(path + "/datafiles")
	createDirectory(path + "/ampfiles")
	createDirectory(path + "/ampfilesBenchmark")
	createDirectory(path + "/diffs")
	createDirectory(path + "/diffsBenchmark")

FILE_LINKS = {
				"keyfiles/1-0-+0+f01500_00_pi-.key": "https://www.dropbox.com/s/e5ip52e3yh5rd1m/1-0-%2B0%2Bf01500_00_pi-.key?dl=1",
				"keyfiles/1-1-+1+f21270_22_pi-.key":"https://www.dropbox.com/s/j9wv0c1x3d4zqtw/1-1-%2B1%2Bf21270_22_pi-.key?dl=1",
				"keyfiles/1-1++1-rho770_01_pi0.key":"https://www.dropbox.com/s/ds1rsujcufgkign/1-1%2B%2B1-rho770_01_pi0.key?dl=1",
				"keyfiles/1-2-+0+f21270_02_pi-.key":"https://www.dropbox.com/s/acpdz1d3zda0drr/1-2-%2B0%2Bf21270_02_pi-.key?dl=1",
				"keyfiles/1-2-+1-rho31690_13_pi0.key":"https://www.dropbox.com/s/c8zxcniwrx090yn/1-2-%2B1-rho31690_13_pi0.key?dl=1",
				"keyfiles/1-2++0-rho31690_23_pi0.key":"https://www.dropbox.com/s/uvnc9jww14rq4i3/1-2%2B%2B0-rho31690_23_pi0.key?dl=1",
				"keyfiles/1-3-+0-rho770_31_pi0.key":"https://www.dropbox.com/s/untnvbxzv5gqs60/1-3-%2B0-rho770_31_pi0.key?dl=1",
				"keyfiles/1-3-+2+rho770_31_pi0.key":"https://www.dropbox.com/s/dzgitymwietqr3l/1-3-%2B2%2Brho770_31_pi0.key?dl=1",
				"keyfiles/1-3++1+rho31690_23_pi0.key":"https://www.dropbox.com/s/538cpu0wf37e9on/1-3%2B%2B1%2Brho31690_23_pi0.key?dl=1",
				"keyfiles/1-6-+2+sigma_60_pi-.key":"https://www.dropbox.com/s/5s0mb7i6qf3dk8z/1-6-%2B2%2Bsigma_60_pi-.key?dl=1",
				"ampfilesBenchmark/[1-,0-,0+]=[f0_1500_0=[pi0[0,0]pi0][0,0]pi-]_binID-0_0.root":"https://www.dropbox.com/s/1xhhxfm24yl5cxp/%5B1-%2C0-%2C0%2B%5D%3D%5Bf0_1500_0%3D%5Bpi0%5B0%2C0%5Dpi0%5D%5B0%2C0%5Dpi-%5D_binID-0_0.root?dl=1",
				"ampfilesBenchmark/[1-,1-,1+]=[f2_1270_0=[pi0[2,0]pi0][2,2]pi-]_binID-0_0.root":"https://www.dropbox.com/s/4nujd7x8ua4zaz4/%5B1-%2C1-%2C1%2B%5D%3D%5Bf2_1270_0%3D%5Bpi0%5B2%2C0%5Dpi0%5D%5B2%2C2%5Dpi-%5D_binID-0_0.root?dl=1",
				"ampfilesBenchmark/[1-,1+,1-]=[rho_770_-=[pi-[1,0]pi0][0,1]pi0]_binID-0_0.root":"https://www.dropbox.com/s/45fq2ofy0val6su/%5B1-%2C1%2B%2C1-%5D%3D%5Brho_770_-%3D%5Bpi-%5B1%2C0%5Dpi0%5D%5B0%2C1%5Dpi0%5D_binID-0_0.root?dl=1",
				"ampfilesBenchmark/[1-,2-,0+]=[f2_1270_0=[pi0[2,0]pi0][0,2]pi-]_binID-0_0.root":"https://www.dropbox.com/s/a9mm0he7v4laeqd/%5B1-%2C2-%2C0%2B%5D%3D%5Bf2_1270_0%3D%5Bpi0%5B2%2C0%5Dpi0%5D%5B0%2C2%5Dpi-%5D_binID-0_0.root?dl=1",
				"ampfilesBenchmark/[1-,2-,1-]=[rho3_1690_-=[pi-[3,0]pi0][1,3]pi0]_binID-0_0.root":"https://www.dropbox.com/s/jco1i2eopxjeqko/%5B1-%2C2-%2C1-%5D%3D%5Brho3_1690_-%3D%5Bpi-%5B3%2C0%5Dpi0%5D%5B1%2C3%5Dpi0%5D_binID-0_0.root?dl=1",
				"ampfilesBenchmark/[1-,2+,0-]=[rho3_1690_-=[pi-[3,0]pi0][2,3]pi0]_binID-0_0.root":"https://www.dropbox.com/s/a3m6u92wrdck9ga/%5B1-%2C2%2B%2C0-%5D%3D%5Brho3_1690_-%3D%5Bpi-%5B3%2C0%5Dpi0%5D%5B2%2C3%5Dpi0%5D_binID-0_0.root?dl=1",
				"ampfilesBenchmark/[1-,3-,0-]=[rho_770_-=[pi-[1,0]pi0][3,1]pi0]_binID-0_0.root":"https://www.dropbox.com/s/puro2171viwosic/%5B1-%2C3-%2C0-%5D%3D%5Brho_770_-%3D%5Bpi-%5B1%2C0%5Dpi0%5D%5B3%2C1%5Dpi0%5D_binID-0_0.root?dl=1",
				"ampfilesBenchmark/[1-,3-,2+]=[rho_770_-=[pi-[1,0]pi0][3,1]pi0]_binID-0_0.root":"https://www.dropbox.com/s/226kazrim3w2dbu/%5B1-%2C3-%2C2%2B%5D%3D%5Brho_770_-%3D%5Bpi-%5B1%2C0%5Dpi0%5D%5B3%2C1%5Dpi0%5D_binID-0_0.root?dl=1",
				"ampfilesBenchmark/[1-,3+,1+]=[rho3_1690_-=[pi-[3,0]pi0][2,3]pi0]_binID-0_0.root":"https://www.dropbox.com/s/vg9xaa2zmdenqtj/%5B1-%2C3%2B%2C1%2B%5D%3D%5Brho3_1690_-%3D%5Bpi-%5B3%2C0%5Dpi0%5D%5B2%2C3%5Dpi0%5D_binID-0_0.root?dl=1",
				"ampfilesBenchmark/[1-,6-,2+]=[sigma0=[pi0[0,0]pi0][6,0]pi-]_binID-0_0.root":"https://www.dropbox.com/s/6pijfga6olpml09/%5B1-%2C6-%2C2%2B%5D%3D%5Bsigma0%3D%5Bpi0%5B0%2C0%5Dpi0%5D%5B6%2C0%5Dpi-%5D_binID-0_0.root?dl=1",
				"rootpwa.config":"https://www.dropbox.com/s/6nenli5ubnd824g/rootpwa.config?dl=1",
				"datafiles/1660.1700.root":"https://www.dropbox.com/s/emg1ehf5felgnkd/1660.1700.root?dl=1",
				"diffsBenchmark/diff0_imag.root":"https://www.dropbox.com/s/v7zih3jt70s2oik/diff0_imag.root?dl=1",
				"diffsBenchmark/diff0_real.root":"https://www.dropbox.com/s/pan3v1xi4vg0t6l/diff0_real.root?dl=1",
				"diffsBenchmark/diff1_imag.root":"https://www.dropbox.com/s/mgfc9a4mrwvh8ra/diff1_imag.root?dl=1",
				"diffsBenchmark/diff1_real.root":"https://www.dropbox.com/s/idoymz1e7m1tp9s/diff1_real.root?dl=1",
				"diffsBenchmark/diff2_imag.root":"https://www.dropbox.com/s/2bhutkbc3o6nxxc/diff2_imag.root?dl=1",
				"diffsBenchmark/diff2_real.root":"https://www.dropbox.com/s/047d7tloholh2by/diff2_real.root?dl=1",
				"diffsBenchmark/diff3_imag.root":"https://www.dropbox.com/s/00x8dt6ipkjwk83/diff3_imag.root?dl=1",
				"diffsBenchmark/diff3_real.root":"https://www.dropbox.com/s/7ef36cso71biosj/diff3_real.root?dl=1",
				"diffsBenchmark/diff4_imag.root":"https://www.dropbox.com/s/jvib64xzt5p4oj7/diff4_imag.root?dl=1",
				"diffsBenchmark/diff4_real.root":"https://www.dropbox.com/s/a19nwwp84pcbxvr/diff4_real.root?dl=1",
				"diffsBenchmark/diff5_imag.root":"https://www.dropbox.com/s/cujdsq4vyrnnsyg/diff5_imag.root?dl=1",
				"diffsBenchmark/diff5_real.root":"https://www.dropbox.com/s/uf65yejjrkmpkpv/diff5_real.root?dl=1",
				"diffsBenchmark/diff6_imag.root":"https://www.dropbox.com/s/nzhpg5n4b8ynvtu/diff6_imag.root?dl=1",
				"diffsBenchmark/diff6_real.root":"https://www.dropbox.com/s/qpz3ojbjxdhdj9p/diff6_real.root?dl=1",
				"diffsBenchmark/diff7_imag.root":"https://www.dropbox.com/s/u4mzmatalynlk1l/diff7_imag.root?dl=1",
				"diffsBenchmark/diff7_real.root":"https://www.dropbox.com/s/kxp3lat7kxh9cay/diff7_real.root?dl=1",
				"diffsBenchmark/diff8_imag.root":"https://www.dropbox.com/s/xklhxlh6tcbvr24/diff8_imag.root?dl=1",
				"diffsBenchmark/diff8_real.root":"https://www.dropbox.com/s/wibdhorphui4t4e/diff8_real.root?dl=1",
				"diffsBenchmark/diff9_imag.root":"https://www.dropbox.com/s/qvidhbtzm0lig7b/diff9_imag.root?dl=1",
				"diffsBenchmark/diff9_real.root":"https://www.dropbox.com/s/ic2zvcvnga17kbr/diff9_real.root?dl=1",
				}

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description="amplitude benchmark")

	parser.add_argument("path", type=str, metavar="<path>", help="path where benchmark should be run (will be created if it doesn't exist)")
	args = parser.parse_args()

	createDirectories(args.path)

	pyRootPwa.utils.printInfo("downloading test files...")
	progressBar = pyRootPwa.utils.progressBar(0, len(FILE_LINKS)+1, sys.stdout)
	progressBar.start()
	currentFile = 0
	for keyFileName in FILE_LINKS:
		urllib.urlretrieve (FILE_LINKS[keyFileName], args.path + "/" + keyFileName)
		currentFile += 1
		progressBar.update(currentFile)

	rootPwaPath = os.path.expandvars("$ROOTPWA")
	configPath = args.path + "/rootpwa.config"
	pyRootPwa.utils.printInfo("creating file manager from config file: '" + configPath + "'.")
	os.system(rootPwaPath + "/pyInterface/userInterface/createFileManager.py -c " + configPath)

	pyRootPwa.utils.printInfo("calculating amplitudes")
	os.system(rootPwaPath + "/pyInterface/userInterface/calcAmplitudes.py -c " + configPath)

	config = pyRootPwa.rootPwaConfig()
	if not config.initialize(configPath):
		pyRootPwa.utils.printErr("loading config file '" + configPath + "' failed. Aborting...")
		sys.exit(1)
	fileManager = pyRootPwa.loadFileManager(config.fileManagerPath)
	n = 0
	for binID in fileManager.getBinIDList():
		for eventsType in [ pyRootPwa.core.eventMetadata.DATA, pyRootPwa.core.eventMetadata.GENERATED, pyRootPwa.core.eventMetadata.ACCEPTED ]:
			for ampFilePathCreated in fileManager.getAmplitudeFilePaths(binID, eventsType):
				ampFilePathBenchmark = args.path + "/ampfilesBenchmark/" + os.path.basename(ampFilePathCreated)
				real1, imag1 = compareAmplitudeFiles.extractRealImagListsFromAmpFile(ampFilePathCreated)
				real2, imag2 = compareAmplitudeFiles.extractRealImagListsFromAmpFile(ampFilePathBenchmark)
				diffReal = compareAmplitudeFiles.getDiffListAbs(real1, real2)
				diffImag = compareAmplitudeFiles.getDiffListAbs(imag1, imag2)

		compareAmplitudeFiles.saveTree(args.path + "/diffs/diff" + str(n) + "_real.root", "tree", "diffs", diffReal)
		compareAmplitudeFiles.saveTree(args.path + "/diffs/diff" + str(n) + "_imag.root", "tree", "diffs", diffImag)
		n += 1

	success = True
	for diffFileName in sorted(glob.glob(args.path + "/diffsBenchmark/*.root")):
		diffFileBenchmark = ROOT.TFile(diffFileName, "READ")
		diffCreatedFileName = args.path + "/diffs/" + os.path.basename(diffFileName)
		diffFileCreated = ROOT.TFile(diffCreatedFileName, "READ")

		keyListBenchmark = diffFileBenchmark.GetListOfKeys()
		if not len(keyListBenchmark) == 1:
			pyRootPwa.utils.printErr("not a valid histogram file: '" + diffFileName + "'. Aborting...")
			sys.exit(1)
		treeBenchmark = diffFileBenchmark.Get(keyListBenchmark[0].GetName())
		arrBenchmark = treeToArray(treeBenchmark)

		keyListCreated = diffFileCreated.GetListOfKeys()
		treeCreated = diffFileCreated.Get(keyListCreated[0].GetName())
		arrCreated = treeToArray(treeCreated)

		if not treeCreated.GetEntries() == treeBenchmark.GetEntries():
			pyRootPwa.utils.printErr("lalala")
			sys.exit(1)
		entryCount = treeCreated.GetEntries()
		result = ROOT.TMath.KolmogorovTest(entryCount, arrCreated, entryCount, arrBenchmark, "")
		pyRootPwa.utils.printInfo("probability: " + str(result*100) + "% ('" + diffFileName + "' vs. '" + diffCreatedFileName + "'.")
		if result < .955:
			success = False

	if success:
		pyRootPwa.utils.printSucc("TEST SUCCEEDED!!")
	else:
		pyRootPwa.utils.printWarn("TEST FAILED!!")
