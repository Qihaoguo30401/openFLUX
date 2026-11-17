"""openFLUX: 3D nonlinear magnetostatic FEA for axial-flux permanent-magnet machines."""
from .materials import BH, add_bh, energy_cf, mu0
from .geometry import build, P as DEFAULT_PARAMS
from .solver import AxialFluxModel

__version__ = '0.1.0'
