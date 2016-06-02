from _printingUtils import printDebug


class multiBin(object):


	def __init__(self, boundaries): # boundaries = { "binningVariable": (lowerBound, upperBound) }
		if not isinstance(boundaries, dict):
			raise TypeError("Boundaries is not of type 'dict'.")
		if not boundaries:
			raise ValueError("Bin boundaries are empty.")
		for key in boundaries.keys():
			if not isinstance(key, str):
				raise TypeError("Binning variable name is not of type 'str'.")
			binRange = boundaries[key]
			if not isinstance(binRange, tuple):
				raise TypeError("Binning range is not of type 'tuple' for binning variable '" + key + "'.")
			if len(binRange) != 2:
				raise ValueError("Binning range does not have two entries for binning variable '" + key + "'.")
			if not (isinstance(binRange[0], float) or isinstance(binRange[0], int)):
				raise TypeError("Lower bound of bin range is not a number for binning variable '" + key + "'.")
			if not (isinstance(binRange[1], float) or isinstance(binRange[1], int)):
				raise TypeError("Upper bound of bin range is not a number for binning variable '" + key + "'.")
			if binRange[0] >= binRange[1]:
				raise ValueError("Lower bound of bin range (" + str(binRange[0]) + ") is larger or equal to upper bound of bin range (" + str(binRange[1]) + ").")
			if isinstance(binRange[0], int) or isinstance(binRange[1], int):
				boundaries[key] = (float(binRange[0]), float(binRange[1]))
		self.boundaries = boundaries


	def __lt__(self, other):
		keys = sorted(self.boundaries.keys())
		if keys != sorted(other.boundaries.keys()):
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


	def overlap(self, other, strict=True):
		def comparator(left, right, direction):
			if direction == ">":
				right, left = left, right
			if strict:
				return left < right
			else:
				return left <= right

		keys = sorted(self.boundaries.keys())
		if keys != sorted(other.boundaries.keys()):
			return False
		for key in keys:
			if comparator(self.boundaries[key][1], other.boundaries[key][0], "<"):
				return False
			if comparator(self.boundaries[key][0], other.boundaries[key][1], ">"):
				return False
		return True


	def sameBinningVariables(self, other):
		return sorted(self.boundaries.keys()) == sorted(other.boundaries.keys())


	def inBin(self, binningInfo):
		# binningInfo = { "variableName": value }
		if sorted(binningInfo.keys()) != sorted(self.boundaries.keys()):
			return False
		for variableName, value in binningInfo.iteritems():
			if value < self.boundaries[variableName][0] or value >= self.boundaries[variableName][1]:
				return False
		return True


def _testMultiBin():
	bins = []
	for i in xrange(3):
		for j in xrange(3):
			bins.append(multiBin( { "mass": (float(i), float(i+1)), "tPrime": (float(j), float(j+1))} ))
	for mBin in bins:
		printDebug(str(mBin) + " <(" + str(mBin < bins[4]) + ") >(" + str(mBin > bins[4]) + ") ol(" + str(mBin.overlap(bins[4])) + ")")
