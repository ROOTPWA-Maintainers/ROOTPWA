
import os
import sys
import tempfile

class silencer():

	output = ""

	_outputStream = None
	_save = None
	_silence = True

	def __init__(self, silence = True):
		self._silence = silence

	def __enter__(self):
		if not self._silence:
			return self
		sys.stdout.flush()
		self._outputStream = tempfile.TemporaryFile().__enter__()
		self._save = os.dup(1), os.dup(2)
		os.dup2(self._outputStream.fileno(), 1)
		os.dup2(self._outputStream.fileno(), 2)
		return self

	def __exit__(self, *args):
		if not self._silence:
			self.output = ""
			return
		os.dup2(self._save[0], 1)
		os.dup2(self._save[1], 2)
		os.close(self._save[0])
		os.close(self._save[1])
		self.output = ""
		self._outputStream.seek(0)
		for line in self._outputStream.readlines():
			self.output += line
		self._outputStream.__exit__()
