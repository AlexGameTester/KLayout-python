__version__ = "8QS_0.0.0.1"

'''
Changes log

8QS_0.0.0.1
Based on 8Q_0.0.0.1
'''

# import built-ins
from typing import List
import os
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
from classLib.josJ import AsymSquidParams, AsymSquid
from classLib.chipDesign import ChipDesign
from classLib.marks import MarkBolgar
from classLib.contactPads import ContactPad
from classLib.helpers import fill_holes, split_polygons, extended_region

import sonnetSim

reload(sonnetSim)
from sonnetSim import SonnetLab, SonnetPort, SimulationBox

# import local dependencies in case project file is too big
from classLib.baseClasses import ComplexBase
from classLib.helpers import FABRICATION
from classLib.chipTemplates import CHIP_14x14_20pads

CHIP = CHIP_14x14_20pads

PROJECT_DIR = os.path.dirname(__file__)

# positioning of md open-end lines relative to qubit
# represents shift from qubit center for CPW's central conductor
# open-ended center
VERT_ARR_SHIFT = DVector(-50e3, -150e3)

from dataclasses import dataclass, field


@dataclass()
class QubitsGrid:
    # in fractions of chip dimensions
    origin: float = DVector(CHIP.dx/2, CHIP.dy/2)
    # step of 2D grid in `x` and `y` directions correspondingly
    dx: float = 2e6
    dy: float = 2e6
    pts_grid: np.ndarray = np.array(
        [
            # grid iterates from left to right, from bottom to top,
            # starting from bl corner.
            (0, 0), (1, 0), (2, 0),
            (0, 1), (1, 1), (2, 1),
            (0, 2), (1, 2),
        ],
        dtype=int
    )

    def __post_init__(self):
        self.__centralize_grid()

    def __centralize_grid(self):
        # grid is centralized such that bbox center of the grid has
        # origin (0,0)
        grid_center = np.array([1, 1], dtype=int)
        for i, pt_pos in enumerate(self.pts_grid):
            self.pts_grid[i] -= grid_center

    def get_pt(self, idx) -> DPoint:
        pt_pos = self.pts_grid[idx]
        pt_x = pt_pos[0] * self.dx
        pt_y = pt_pos[1] * self.dy
        origin = self.origin
        return origin + DVector(pt_x, pt_y)


@dataclass()
class DiskConn8Pars:
    disk_r = 0.5e6
    pimp_l = 100e3
    conn_width = 20e3
    conn_gap = 10e3


class DiskConn8(ComplexBase):
    """
    Single superconducting Island represents shunting ground capacitor
    for qubit
    """

    def __init__(self, origin,
                 pars: DiskConn8Pars = DiskConn8Pars(),
                 trans_in=None,
                 region_id="ph"):
        self.pars: DiskConn8Pars = pars
        self.disk: DiskConn8 = None
        self.conn8_list: List[CPW] = []
        super().__init__(
            origin=origin, trans_in=trans_in, region_id=region_id
        )

    def init_primitives(self):
        origin = DPoint(0, 0)

        # draw star-like connection flanges
        angles = np.linspace(0, 360, 8, endpoint=False)  # degree
        for i, angle in enumerate(angles):
            cpw_l = self.pars.disk_r + self.pars.pimp_l
            cpw = CPW(
                start=origin,
                end=origin + DVector(cpw_l, 0),
                width=self.pars.conn_width,
                gap=self.pars.conn_gap,
                trans_in=DCplxTrans(1, angle, False, 0, 0),
                region_id=self.region_id
            )
            self.conn8_list.append(cpw)
            self.primitives["conn" + str(i)] = cpw

        from classLib.shapes import Disk
        self.disk = Disk(
            center=origin, r=self.pars.disk_r,
            region_id=self.region_id
        )
        self.primitives["circle"] = self.disk


class QubitParams:
    def __init__(
            self,
            squid_params: AsymSquidParams = AsymSquidParams(),
            qubit_cap_params: DiskConn8Pars = DiskConn8Pars()
    ):
        self.squid_params: AsymSquidParams = squid_params
        self.qubit_cap_params: DiskConn8Pars = qubit_cap_params


class Qubit(ComplexBase):
    def __init__(
            self,
            origin: DPoint = DPoint(0, 0),
            qubit_params: QubitParams = QubitParams(),
            trans_in=None
    ):
        self.qubit_params = qubit_params
        self.squid: AsymSquid = None
        self.cap_shunt: DiskConn8 = None

        super().__init__(origin=origin, trans_in=trans_in,
                         region_id="default")

    def init_primitives(self):
        origin = DPoint(0,0)
        self.squid = AsymSquid(
            origin=origin,
            params=self.qubit_params.squid_params,
            region_id="el"
        )
        self.primitives["squid"] = self.squid

        self.cap_shunt = DiskConn8(
            origin=origin,
            pars=self.qubit_params.qubit_cap_params,
            region_id="ph"
        )
        self.primitives["cap_shunt"] = self.cap_shunt


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

    def draw_chip(self):
        self.region_bridges2.insert(self.chip_box)

        self.region_ph.insert(self.chip_box)
        for contact_pad in self.contact_pads:
            contact_pad.place(self.region_ph)

    def draw_qubits_array(self):
        # draw rectangles to check QubitsGrid() class
        # for i, _ in enumerate(self.qubits_grid.pts_grid):
        #     p_center = self.qubits_grid.get_pt(idx=i)
        #     dv = 1e6 * DVector(0.1, 0.1)
        #     box = DBox(p_center - dv, p_center + dv)
        #     self.region_el.insert(box)

        for pt_i in range(len(self.qubits_grid.pts_grid)):
            pt = self.qubits_grid.get_pt(pt_i)
            qubit = Qubit(origin=pt)
            self.qubits.append(qubit)
            qubit.place(self.region_ph, region_id="ph")
            qubit.place(self.region_el, region_id="el")


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
