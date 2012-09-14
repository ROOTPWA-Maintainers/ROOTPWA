
import multiprocessing

class Keyfile():

	lock = multiprocessing.Lock()
	fileName = ""

	def __init__(self, fileName):
		self.fileName = fileName

	def __str__(self):
		return self.fileName
