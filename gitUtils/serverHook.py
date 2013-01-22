#!/usr/bin/env python

import sys
import subprocess

refname = sys.argv[1]
oldrev = sys.argv[2]
newrev = sys.argv[3]

if oldrev == "0000000000000000000000000000000000000000":
	allCommits = [newrev]
else:
	allCommits = subprocess.check_output("git rev-list " + oldrev + ".." + newrev, shell=True).split('\n')

while '' in allCommits:
	allCommits.remove('')

for sha in allCommits:
	commitInfo = subprocess.check_output("git cat-file commit " + sha, shell=True).split('\n')
	commitMessage = []
	for i in range(len(commitInfo)):
		if commitInfo[i] == '':
			commitMessage = commitInfo[i+1:]
			break
		i += 1

	if len(commitMessage) < 3:
		sys.stderr.write("[POLICY] Commit [" + sha[:15] + "]: Commit message has less than 3 lines. Aborting push...\n")
		sys.exit(1)
	if len(commitMessage[0]) > 80 or len(commitMessage[0]) == 0:
		sys.stderr.write("[POLICY] Commit [" + sha[:15] + "]: Commit message's first line is longer than 80 characters. Aborting push...\n")
		sys.exit(1)
	if len(commitMessage[1]) != 0:
		sys.stderr.write("[POLICY] Commit [" + sha[:15] + "]: Commit message's second line is not empty. Aborting push...\n")
		sys.exit(1)
	if len(commitMessage[2]) == 0:
		sys.stderr.write("[POLICY] Commit [" + sha[:15] + "]: Commit message's third line must not be empty. Aborting push...\n")
		sys.exit(1)
	for line in commitMessage[2:]:
		if len(line) > 80:
			sys.stderr.write("[POLICY] Commit [" + sha[:15] + "]: Commit message has a line longer than 80 characters. Aborting push...\n")
			sys.exit(1)

