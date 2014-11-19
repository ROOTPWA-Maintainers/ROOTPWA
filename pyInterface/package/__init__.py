
import ROOT

from _config import rootPwaConfig
from _fileManager import fileManager
from _fileManager import saveFileManager
from _fileManager import loadFileManager
from amplitude import calcAmplitude
from integrals import calcIntegrals

del _config
del _fileManager
del amplitude
del integrals

config = None
