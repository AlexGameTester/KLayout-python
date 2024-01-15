# import built-ins
from typing import List
from importlib import reload
import os
from dataclasses import dataclass

# import good 3rd party
import numpy as np

# import project specific 3rd party
import pya
from pya import DVector

# import self-made API
import classLib

reload(classLib)
from classLib.josJ import AsymSquidParams
from classLib.chipTemplates import CHIP_14x14_20pads

# chip external couplings and dimensions
CHIP = CHIP_14x14_20pads

# positioning of md open-end lines relative to qubit
# represents shift from qubit center for CPW's central conductor
# open-ended center
VERT_ARR_SHIFT = DVector(-50e3, -150e3)

# directory to output files
PROJECT_DIR = os.path.dirname(__file__)

# print("global definitions loaded")
BC_dx = 2.5e3 * np.sqrt(2) + 1e3
SQUID_PARS = AsymSquidParams(
    squid_dx=14.2e3, squid_dy=10e3,
    TC_dx=2.5e3 * np.sqrt(2) + 1e3,
    TC_dy=5e3 * np.sqrt(2) / 2 + 1e3,
    BC_dx=(BC_dx, BC_dx),
    BC_dy=5e3 * np.sqrt(2) / 2 + 1e3,
    TCW_dx=0.5e3, TCW_dy=5.232e3, BCW_dy=0e3
)
