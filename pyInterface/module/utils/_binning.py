import copy
from _printingUtils import printDebug, printErr


class multiBin(object):


	def __init__(self, boundaries): # boundaries = { "binningVariable": (lowerBound, upperBound) }
		if not isinstance(boundaries, dict):
			raise TypeError("Boundaries is not of type 'dict'.")
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
			if binRange[0] > binRange[1]:
				raise ValueError("Lower bound of bin range (" + str(binRange[0]) + ") is larger than upper bound of bin range (" + str(binRange[1]) + ").")
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
		retval = "multiBin: " + self.__repr__()
		return retval

	def __repr__(self):
		retVariables = []
		for key in sorted(self.boundaries.keys()):
			retVariables.append( "\"" + key + "\": (" + str(self.boundaries[key][0]) + ", " + str(self.boundaries[key][1]) + ")")
		retval = "{ "
		retval += ", ".join(retVariables)
		retval += " }"
		return retval

	def __hash__(self):
		keys = tuple(self.variables())
		values = tuple([ self.boundaries[k] for k in keys])
		return hash((keys, values))

	def __contains__(self, other):
		'''
		The other multiBin is contained in this multiBin if the intervals of all variables
		of this multiBin are larger than the corresponding intervals of the other multiBin.
		If this multiBin is binned in a variable in which the other multiBin is not binned,
		the other multiBin is not contained in this multiBin.
		If the other multiBin is binned in a variable in which this multiBin is not binned,
		the other multiBin is contained in this multiBin.
		'''
		if isinstance(other, dict):
			other = multiBin(other)
		if not isinstance(other, multiBin):
			msg = "Cannot check whether object of type {0} is contained in a multiBin".format(type(other))
			printErr(msg)
			raise ValueError(msg)

		for key in self.boundaries.keys():
			if key in other.boundaries.keys():
				if not self.boundaries[key][0] <= other.boundaries[key][0]:
					return False
				if not self.boundaries[key][1] >= other.boundaries[key][1]:
					return False
			else:
				return False
		return True

	def variables(self):
		'''
		@return: sorted list of variables defined in the multibin
		'''
		return sorted(self.boundaries.keys())

	def uniqueStr(self):
		'''
		@return: string, unique for the multibin (without whitespaces)
		'''
		out = []
		for variable in self.variables():
			variableStr = "{0}_{1!r}_{2!r}".format(variable, self.boundaries[variable][0], self.boundaries[variable][1])
			out.append(variableStr)
		return "__".join(out)

	@classmethod
	def fromUniqueStr(cls, strIn):
		boundaries = {}
		if strIn:
			variables = strIn.split('__')
			for variable in variables:
				if variable.count('_') == 2:
					name, lower, upper = variable.split("_")
					boundaries[name] = (float(lower), float(upper))
				else:
					printErr("Cannot get multiBin form string '{0}'".format(strIn))
		return multiBin(boundaries)

	def overlap(self, other, strict=True):
		'''
		@return: True if this and the other multiBin overlap in all variables, False otherwise.
		'''
		keys = self.variables()
		if keys != other.variables():
			return False
		for key in keys:
			if not self.overlapInVariable(other, key, strict):
				return False
		return True


	def overlapInVariable(self, other, key, strict=True):
		'''
		@return: True if this and the other multiBin overlap in variable key, False otherwise.
		'''
		def comparator(left, right, direction):
			if direction == ">":
				right, left = left, right
			if strict:
				return left < right
			else:
				return left <= right

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

	def getSubMultiBin(self, exception=None):
		'''
		@param exception: Str or list of Str with the variables which should be removed from the returned multiBin object
		@return: multiBin in all binning variables except those given in exception
		@rtype: multiBin
		'''
		if exception is None:
			return copy.deepcopy(self)

		if isinstance(exception, str):
			exception = [exception]
		boundaries = {}
		for k,values in self.boundaries.iteritems():
			if k not in exception:
				boundaries[k] = values
		return multiBin(boundaries)

	def getBinCenters(self):
		return {k: 0.5*(v[0]+v[1]) for k,v in self.boundaries.iteritems()}


def _testMultiBin():
	bins = []
	for i in xrange(3):
		for j in xrange(3):
			bins.append(multiBin( { "mass": (float(i), float(i+1)), "tPrime": (float(j), float(j+1))} ))
	for mBin in bins:
		printDebug(str(mBin) + " <(" + str(mBin < bins[4]) + ") >(" + str(mBin > bins[4]) + ") ol(" + str(mBin.overlap(bins[4])) + ")")
