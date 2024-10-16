from pathlib import Path
import sys
CodePath = Path('/Users/charlesh/Documents/Codes/');sys.path.insert(1,CodePath)
ProjectPath = Path(CodePath/'OBS_Methods/NOISE');sys.path.insert(1,ProjectPath)
MethodSector,CompSector = (ProjectPath/'METHODS',ProjectPath/'COMPS')
sys.path.append(str(MethodSector));sys.path.insert(0, str(CompSector))
from OrientPy import *
# sys.path.append('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/METHODS')
k = 1