#!/usr/bin/env python2.6

import os
import sys
import subprocess

refname = sys.argv[1]
oldrev = sys.argv[2]
newrev = sys.argv[3]

def runCommand(command):
	process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	retval = process.communicate()[0].split('\n')
	if process.returncode != 0:
		print("Something went wrong in the update hook. Please contact a ROOTPWA administrator...")
		print("Add the following information to your report:")
		print("Command: \"" + command + "\"")
		print("Output:\n" + str(retval))
		sys.exit(255)
	return retval


allowedUsers = [ ]

branchOrTagName = refname.rsplit('/', 1)[-1]

if newrev == "0000000000000000000000000000000000000000":
	if branchOrTagName.startswith('_') or branchOrTagName == "master":
		sys.stderr.write("Branch '" + refname + "' cannot be deleted. Aborting push...\n")
		sys.exit(1)
	allCommits = [oldrev]
elif oldrev == "0000000000000000000000000000000000000000":
	if branchOrTagName.startswith('_') or branchOrTagName == "master":
		if os.environ['USER'] not in allowedUsers:
			sys.stderr.write("Branch '" + refname + "' cannot be created. Aborting push...\n")
			sys.exit(1)
	allCommits = [newrev]
else:
	command = "git rev-list " + oldrev + ".." + newrev
	allCommits = runCommand(command)

while '' in allCommits:
	allCommits.remove('')

for sha in allCommits:
	command = "git cat-file commit " + sha
	commitInfo = runCommand(command)
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

