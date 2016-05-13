from _printingUtils import printDebug

class multiBin(object):

	def __init__(self, boundaries):
		self.boundaries = boundaries # { "binningVariable": (lowerBound, upperBound) }


	def __lt__(self, other):
		keys = sorted(self.boundaries.keys())
		if keys != sorted(self.boundaries.keys()):
			return False
		for key in keys:
			selfCenter = (self.boundaries[key][0] + self.boundaries[key][1]) / 2.
			otherCenter = (other.boundaries[key][0] + other.boundaries[key][1]) / 2.
			if selfCenter > otherCenter:
				return False
			elif selfCenter < otherCenter:
				return True
		return False


	def __le__(self, other):
		return self < other or self == other


	def __eq__(self, other):
		return self.boundaries == other.boundaries


	def __ne__(self, other):
		return not self == other


	def __gt__(self, other):
		return not self <= other


	def __ge__(self, other):
		return not self < other


	def __str__(self):
		retval = ""
		retval += "multiBin: { "
		for key in sorted(self.boundaries.keys()):
			retval += "\"" + key + "\": (" + str(self.boundaries[key][0]) + ", " + str(self.boundaries[key][1]) + "), "
		retval = retval[:-2] + " }"
		return retval


	def overlap(self, other):
		keys = sorted(self.boundaries.keys())
		if keys != sorted(self.boundaries.keys()):
			return False
		for key in keys:
			if self.boundaries[key][1] <= other.boundaries[key][0]:
				return False
			if self.boundaries[key][0] >= other.boundaries[key][1]:
				return False
		return True


def _testMultiBin():
	bins = []
	for i in xrange(3):
		for j in xrange(3):
			bins.append(multiBin( { "mass": (float(i), float(i+1)), "tPrime": (float(j), float(j+1))} ))
	for mBin in bins:
		printDebug(str(mBin) + " <(" + str(mBin < bins[4]) + ") >(" + str(mBin > bins[4]) + ") ol(" + str(mBin.overlap(bins[4])) + ")")
