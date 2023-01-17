__version__ = "8QS_0.0.0.1"

'''
Changes log

Design is based on schematics located on YandexDisk: 
https://disk.yandex.com/d/F1Uz4Qk79VytSA
Developing journal is also available on YandexDisk:
https://disk.yandex.com/i/KLZmyRAYXG4mGA
'''

# import built-ins
from typing import List
import os
import itertools

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
from classLib.coplanars import CPW, CPWParameters, DPathCPW
from classLib.chipDesign import ChipDesign
from classLib.marks import MarkBolgar
from classLib.contactPads import ContactPad
from classLib.helpers import fill_holes, split_polygons, extended_region
from classLib.helpers import simulate_cij, save_sim_results
from classLib.shapes import RingSector

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
import designElementsGeometry

reload(designElementsGeometry)
from designElementsGeometry import QubitParams, Qubit, QubitsGrid
from designElementsGeometry import DiskConn8Pars


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

        ''' QUBITS GRID '''
        self.qubits_grid: QubitsGrid = QubitsGrid()
        self.qubits_n = len(self.qubits_grid.pts_grid)
        self.qubits: List[Qubit] = []
        ''' QUBIT COUPLINGS '''
        self.q_couplings: np.array = np.empty(
            (self.qubits_n, self.qubits_n),
            dtype=object
        )
        self.circle_hull_d = 5e3

        ''' READOUT LINES '''
        self.ro_lines: List[DPathCPW] = [None, None]
        self.qCenter_roLine_distance = 1e6

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
        self.draw_qq_couplings()
        self.draw_readout_lines()

    def draw_postpone(self):
        """
        Placing elements that were scheduled for postpone drawing.
        Initially implemented to use for geometry parameters sweeps in
        order to be able to change geometry parameters after
        `ChipDesign` class instance is initialized.

        Returns
        -------

        """
        for qubit in self.qubits:
            qubit.place(self.region_ph, region_id="ph")
            qubit.place(self.region_el, region_id="el")

    def draw_chip(self):
        self.region_bridges2.insert(self.chip_box)

        self.region_ph.insert(self.chip_box)
        for contact_pad in self.contact_pads:
            contact_pad.place(self.region_ph)

    def draw_qubits_array(self):
        def get_squid_connector_idx(qubit_idx: int):
            """
            Returns squid connector idx based on qubit idx
            Parameters
            ----------
            qubit_idx: int
                from 0 to 7

            Returns
            -------
            int
                qubit disk's connector idx
            """
            if qubit_idx in [0, 1, 3, 6]:
                return 5
            else:
                return 1

        for qubit_idx in range(len(self.qubits_grid.pts_grid)):
            pt = self.qubits_grid.get_pt(qubit_idx)
            qubit_pars = QubitParams(
                squid_params=SQUID_PARS,
                qubit_cap_params=DiskConn8Pars(),
                squid_connector_idx=get_squid_connector_idx(qubit_idx)
            )
            qubit = Qubit(
                origin=pt,
                qubit_params=qubit_pars,
                postpone_drawing=False
            )
            self.qubits.append(qubit)
            # TODO: not working, qubit.squid.origin is wrong | partially
            #  dealt with. Check this.
            # shift squid to suit into scheme
            # qubit.squid.make_trans(DCplxTrans(1, 0, False, 0, -20e3))
            # q_origin = qubit.origin.dup()  # memorize origin
            # # transfer to origin
            # qubit.make_trans(DCplxTrans(1, 0, False, -q_origin))
            # # rotate depending on qubit group
            # if pt_i in [0, 1, 3]:
            #     qubit.make_trans(DCplxTrans(1, -45, False, 0, 0))
            # else:
            #     qubit.make_trans(DCplxTrans(1, -45, True, 0, 0))
            # qubit.make_trans(DCplxTrans(1, 0, False, q_origin))
            qubit.place(self.region_ph, region_id="ph")
            qubit.place(self.region_el, region_id="el")

    def draw_qq_couplings(self):
        it_1d = list(range(len(self.qubits_grid.pts_grid)))
        it_2d = itertools.product(it_1d, it_1d)
        for pt_i, pt_j in it_2d:
            if (pt_i < pt_j):
                continue
            row_i = pt_i // 3
            col_i = pt_i % 3
            row_j = pt_j // 3
            col_j = pt_j % 3
            if (abs(row_i - row_j) + abs(col_i - col_j)) == 1:
                pt_1 = self.qubits_grid.get_pt(pt_i)
                pt_2 = self.qubits_grid.get_pt(pt_j)
                dv = (pt_1 - pt_2)
                dv = dv/dv.abs()
                dv = dv*(
                    self.qubits[row_i*3 + col_i].cap_shunt.pars.disk_r
                    + self.circle_hull_d
                )
                pt_1 -= dv
                pt_2 += dv
                cpw = CPW(start=pt_1, end=pt_2, width=40e3, gap=10e3)
                cpw.place(self.region_ph)
                self.q_couplings[row_i, row_j] = cpw

    def draw_readout_lines(self):
        # left readout line
        p0_start = self.contact_pads[-1].end
        p0_end = self.contact_pads[8].end
        p1 = p0_start + DVector(0, -3e6)
        p2 = DPoint(self.qubits[3].origin.x - self.qCenter_roLine_distance, p1.y)
        p3 = DPoint(p2.x, self.qubits[0].origin.y - self.qCenter_roLine_distance)
        p4 = DPoint(self.qubits[2].origin.x + self.qCenter_roLine_distance, p3.x)
        p5 = DPoint(p0_end.x, p4.y)
        pts = [p0_start, p1, p2, p3, p4, p5, p0_end]
        self.ro_lines[0] = DPathCPW(
            points=pts,
            cpw_parameters=[CPWParameters(width=20e3, gap=10e3)],
            turn_radii=[60e3],
            trans_in=None,
            region_id="ph"
        )
        self.ro_lines[0].place(self.region_ph, region_id="ph")

        # right readout line
        

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


def simulate_Cqq(q1_idx, q2_idx=-1, resolution=(5e3, 5e3)):
    resolution_dx, resolution_dy = resolution
    dl_list = [-5e3, 0, 5e3]
    # x_distance_dx_list = [0]
    for dl in dl_list:
        ''' DRAWING SECTION START '''
        design = Design8QStair("testScript")
        design.draw()

        # design.N_coils = [1] * design.NQUBITS

        q1 = design.qubits[q1_idx]
        q1.qubit_params.qubit_cap_params.disk_r = 120e3 + dl
        # q2 = design.qubits[q2_idx]
        design.draw_postpone()

        design.show()
        design.lv.zoom_fit()
        design.layout.write(
            os.path.join(PROJECT_DIR, f"Cqq_{q1_idx}_{q2_idx}_"
                                      f"{dl:.3f}_.gds")
        )
        '''DRAWING SECTION END'''
        q1_reg = q1.metal_regions["ph"]
        C1, _, _ = simulate_cij(
            design, env_reg=None, layer=design.layer_ph,
            subregs=[q1_reg],
            resolution=resolution, print_values=True
        )

        '''SAVING REUSLTS SECTION START'''
        output_filepath = os.path.join(
            PROJECT_DIR,
            f"Xmon_Cqq_{q1_idx}_{q2_idx}_results.csv"
        )
        save_sim_results(
            output_filepath=output_filepath,
            design=design,
            additional_pars={"C1, fF": C1}
        )


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
    # simulate_Cqq(0, resolution=(2e3, 2e3))

    ''' MD line C_qd for md1,..., md6 '''
    # for md_idx in [0,1]:
    #     for q_idx in range(2):
    #         simulate_md_Cg(md_idx=md_idx, q_idx=q_idx, resolution=(1e3, 1e3))

    ''' Resonators Q and f sim'''
    # simulate_resonators_f_and_Q(resolution=(2e3, 2e3))

    ''' Resonators Q and f when placed together'''
    # simulate_resonators_f_and_Q_together()
