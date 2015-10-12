
import ROOT

from _config import rootPwaConfig
from _fileManager import fileManager
from _fileManager import saveFileManager
from _fileManager import loadFileManager
from amplitude import calcAmplitude
from fit import pwaFit
from fit import pwaNloptFit
from integrals import calcIntegrals
from likelihood import initLikelihood

del _config
del _fileManager
del amplitude
del integrals
del fit

config = None
