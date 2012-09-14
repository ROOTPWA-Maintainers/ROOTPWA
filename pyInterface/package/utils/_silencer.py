
import os
import sys
import tempfile

class Silencer():

	_outputStream = None
	_save = None

	output = ""

	def __init__(self):
		pass

	def __enter__(self):
		sys.stdout.flush()
		self._outputStream = tempfile.SpooledTemporaryFile(1000000, 'rw')
		self._save = os.dup(1), os.dup(2)
		os.dup2(self._outputStream.fileno(), 1)
		os.dup2(self._outputStream.fileno(), 2)
		return self

	def __exit__(self, *args):
		os.dup2(self._save[0], 1)
		os.dup2(self._save[1], 2)
		self.output = ""
		self._outputStream.seek(0)
		for line in self._outputStream.readlines():
			self.output += line
		self._outputStream.close()
