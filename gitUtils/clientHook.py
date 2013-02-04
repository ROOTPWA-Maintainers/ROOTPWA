#!/usr/bin/env python

import sys
import os


if __name__ == "__main__":

	msgFileName = sys.argv[1]

	with open(msgFileName, 'r') as msgFile:
		prelimCommitMessage = msgFile.readlines()

	commitMessage = []
	for line in prelimCommitMessage:
		if line.startswith('#'):
			continue
		commitMessage.append(line.rstrip('\n'))

	error = False
	if len(commitMessage) < 3:
		sys.stderr.write("[POLICY] Commit message has less than 3 lines. Aborting commit...\n")
		error = True
	if not error and (len(commitMessage[0]) > 80 or len(commitMessage[0]) == 0):
		sys.stderr.write("[POLICY] Commit message's first line is longer than 80 characters. Aborting commit...\n")
		error = True
	if not error and (len(commitMessage[1]) != 0):
		sys.stderr.write("[POLICY] Commit message's second line is not empty. Aborting commit...\n")
		error = True
	if not error and (len(commitMessage[2]) == 0):
		sys.stderr.write("[POLICY] Commit message's third line must not be empty. Aborting commit...\n")
		error = True
	if not error:
		for line in commitMessage[2:]:
			if len(line) > 80:
				sys.stderr.write("[POLICY] Commit message has a line longer than 80 characters. Aborting commit...\n")
				error = True

	if error:
		i = 1
		commitMessageFileName = "git-commit.tmp"
		while os.path.isfile(commitMessageFileName):
			commitMessageFileName = "git-commit_" + str(i) + ".tmp"
			i += 1
		os.system("cp " + msgFileName + " " + commitMessageFileName)
		print("Could not commit, your commit message has been saved to '" + commitMessageFileName + "'")
		sys.exit(1)
