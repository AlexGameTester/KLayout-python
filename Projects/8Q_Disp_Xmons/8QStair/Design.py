__version__ = "8QS_0.0.0.1"

'''
Changes log

8QS_0.0.0.1
Based on 8Q_0.0.0.1
'''

# import built-ins
from typing import List
import os
PROJECT_DIR = os.path.dirname(__file__)
import sys
sys.path.append(PROJECT_DIR)
from importlib import reload

# import good 3rd party
import numpy as np

# import project specific 3rd party
import pya
from pya import Point, Vector, DPoint, DVector, DEdge, \
    DSimplePolygon, \
    SimplePolygon, DPolygon, DBox, Polygon, Region

from pya import DCplxTrans

# import project lib
import classLib
reload(classLib)
from classLib.coplanars import CPW, CPWParameters
from classLib.chipDesign import ChipDesign
from classLib.marks import MarkBolgar
from classLib.contactPads import ContactPad
from classLib.helpers import fill_holes, split_polygons, extended_region

# import sonnet simulation self-made package
import sonnetSim
reload(sonnetSim)
from sonnetSim import SonnetLab, SonnetPort, SimulationBox

# import local dependencies in case project file is too big
from classLib.baseClasses import ComplexBase
from classLib.helpers import FABRICATION
import globalDefinitions
reload(globalDefinitions)
from globalDefinitions import CHIP, VERT_ARR_SHIFT, SQUID_PARS
import qubitGeomPars
reload(qubitGeomPars)
from qubitGeomPars import QubitParams, Qubit, QubitsGrid
from qubitGeomPars import DiskConn8Pars


class Design8QStair(ChipDesign):
    def __init__(self, cell_name):
        super().__init__(cell_name)

        ''' DEFINE LYTOGRAPHY LAYERS AND REGIONS '''
        # `Region` objects are used as intermediate buffers for actual
        # `cell.layer`structures. This optimizes logical operations
        # between different layers.

        # for DC contact deposition
        dc_bandage_layer_i = pya.LayerInfo(3, 0)
        self.dc_bandage_reg = Region()
        self.dc_bandage_layer = self.layout.layer(dc_bandage_layer_i)

        info_bridges1 = pya.LayerInfo(4, 0)  # bridge photo layer 1
        self.region_bridges1 = Region()
        self.layer_bridges1 = self.layout.layer(info_bridges1)

        info_bridges2 = pya.LayerInfo(5, 0)  # bridge photo layer 2
        self.region_bridges2 = Region()
        self.layer_bridges2 = self.layout.layer(info_bridges2)

        # layer with rectangles that will be converted to the CABL format
        info_el_protection = pya.LayerInfo(6, 0)
        self.region_el_protection = Region()
        self.layer_el_protection = self.layout.layer(info_el_protection)

        # has to call it once more to add new layers
        # defined in this child of `ChipDesign`
        self.lv.add_missing_layers()

        ''' CHIP PARAMETERS '''
        self.chip = CHIP
        self.chip_box: pya.DBox = self.chip.box
        # Z = 49.5656 E_eff = 6.30782 (E = 11.45)
        self.z_md_fl: CPWParameters = CPWParameters(10e3, 5.7e3)

        self.contact_pads: List[ContactPad] = self.chip.get_contact_pads(
            [self.z_md_fl] * 16 + [self.chip.chip_Z] * 4,
            back_metal_gap=200e3,
            back_metal_width=0e3,
            pad_length=700e3,
            transition_len=250e3
        )

        ''' QUBITS GRID PARAMETERS '''
        self.qubits_grid: QubitsGrid = QubitsGrid()
        self.qubits: List[Qubit] = []

        ''' QUBIT DRAWING PARAMETERS '''
        ''' SQUID POSITIONING AND PARAMETERS SECTION STARAT '''
        self.squid_vertical_shifts_list = []


    def draw(self):
        """

        Parameters
        ----------

        Returns
        -------
        None
        """
        self.draw_chip()
        self.draw_qubits_array()
        self.draw_readout_lines()

    def draw_chip(self):
        self.region_bridges2.insert(self.chip_box)

        self.region_ph.insert(self.chip_box)
        for contact_pad in self.contact_pads:
            contact_pad.place(self.region_ph)

    def draw_qubits_array(self):
        for pt_i in range(len(self.qubits_grid.pts_grid)):
            pt = self.qubits_grid.get_pt(pt_i)
            qubit_pars = QubitParams(squid_params=SQUID_PARS,
                                     qubit_cap_params=DiskConn8Pars())
            qubit = Qubit(origin=pt, qubit_params=qubit_pars)
            self.qubits.append(qubit)
            qubit.place(self.region_ph, region_id="ph")
            qubit.place(self.region_el, region_id="el")

    def draw_readout_lines(self):
        pass

    def _transfer_regs2cell(self):
        '''

        Returns
        -------

        '''
        # this method is used after all drawing
        # have placed their objects on related regions.
        # This design avoids extensive copy/pasting of polygons to/from
        # `cell.shapes` structure.
        self.cell.shapes(self.layer_ph).insert(self.region_ph)
        self.cell.shapes(self.layer_el).insert(self.region_el)
        self.cell.shapes(self.dc_bandage_layer).insert(self.dc_bandage_reg)
        self.cell.shapes(self.layer_bridges1).insert(self.region_bridges1)
        self.cell.shapes(self.layer_bridges2).insert(self.region_bridges2)
        self.cell.shapes(self.layer_el_protection).insert(
            self.region_el_protection)
        self.lv.zoom_fit()


if __name__ == "__main__":
    ''' draw and show design for manual design evaluation '''
    FABRICATION.OVERETCHING = 0.0e3
    design = Design8QStair("testScript")
    design.draw()
    design.show()

    # design.save_as_gds2(
    #     os.path.join(
    #         PROJECT_DIR,
    #         "Dmon_" + __version__ + "_overetching_0um.gds"
    #     )
    # )

    # FABRICATION.OVERETCHING = 0.5e3
    # design = Design8QStair("testScript")
    # design.draw()
    # design.show()
    # design.save_as_gds2(
    #     os.path.join(
    #         PROJECT_DIR,
    #         "Dmon_" + __version__ + "_overetching_0um5.gds"
    #     )
    # )

    ''' C_qr sim '''
    # simulate_Cqr(resolution=(1e3, 1e3), mode="Cqr", pts=11, par_d=10e3)
    # simulate_Cqr(resolution=(1e3, 1e3), mode="Cq", pts=3, par_d=20e3)
    # simulate_Cqr(resolution=(1e3, 1e3), mode="Cqr")

    ''' Simulation of C_{q1,q2} in fF '''
    # simulate_Cqq(2, 3, resolution=(1e3, 1e3))

    ''' MD line C_qd for md1,..., md6 '''
    # for md_idx in [0,1]:
    #     for q_idx in range(2):
    #         simulate_md_Cg(md_idx=md_idx, q_idx=q_idx, resolution=(1e3, 1e3))

    ''' Resonators Q and f sim'''
    # simulate_resonators_f_and_Q(resolution=(2e3, 2e3))

    ''' Resonators Q and f when placed together'''
    # simulate_resonators_f_and_Q_together()
