
import sys

class progressBar:

	minimum = 0.0
	maximum = 100.0
	stars = 0
	fp = None
	full = False
	effRange = 0.

	def __init__(self, minimum = 0, maximum = 100, fp=sys.stdout):
		self.fp = fp
		self.reset(minimum, maximum - 1)

	def reset(self, minimum = 0, maximum = 100):
		self.full = False
		self.minimum = float(minimum)
		self.maximum = float(maximum - 1)
		try:
			self.effRange = 51. / (self.maximum - self.minimum)
		except ZeroDivisionError:
			self.effRange = 1.
			self.full = True
		self.minimum = self.minimum - (0.5 / self.effRange)
		self.stars = 0

	def start(self):
		self.fp.write("0%   10   20   30   40   50   60   70   80   90   100%\n")
		self.fp.write("|----|----|----|----|----|----|----|----|----|----|\n")
		if self.full:
			self.fp.write(51 * "*" + "\n")
			self.fp.flush()

	def cancel(self):
		self.fp.write("<!\n")

	def update(self, prog):
		if not self.full:
			stars = int((prog - self.minimum) * self.effRange)
			diff = stars - self.stars
			if diff > 0:
				self.fp.write(diff * "*")
				self.fp.flush()
				self.stars += diff
				if self.stars == 51:
					self.fp.write('\n')
					self.full = True
