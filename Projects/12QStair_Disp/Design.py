__version__ = "8QS_0.0.0.1"

'''
Changes log

Design is based on schematics located on YandexDisk: 
https://disk.yandex.com/d/F1Uz4Qk79VytSA
Developing journal is also available on YandexDisk:
https://disk.yandex.com/i/KLZmyRAYXG4mGA

# TODO: change `self.origin` behaviour such it is always shows
objects local coordinate system origin.
Example: 
    "get rid of this shite, this is a plague. Major refactoring will be needed"
    def _refresh_named_connections(self):
        self.origin = self.connections[0]
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

from pya import DCplxTrans, Trans

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
from classLib.resonators import EMResonatorTL3QbitWormRLTailXmonFork

# import local dependencies in case project file is too big
from classLib.baseClasses import ComplexBase
from classLib.helpers import FABRICATION

import globalDefinitions

reload(globalDefinitions)
from globalDefinitions import CHIP, VERT_ARR_SHIFT, SQUID_PARS, PROJECT_DIR

import designElementsGeometry

reload(designElementsGeometry)
from designElementsGeometry import QubitParams, Qubit, QubitsGrid, ResonatorParams
from designElementsGeometry import DiskConn8Pars

import Cqq_couplings

reload(Cqq_couplings)
from Cqq_couplings import CqqCouplingType2, CqqCouplingParamsType2
from Cqq_couplings import CqqCouplingType1, CqqCouplingParamsType1


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
        self.qubits: List[Qubit] = [] * self.qubits_n
        ''' QUBIT COUPLINGS '''
        self.q_couplings: np.array = np.empty(
            (self.qubits_n, self.qubits_n),
            dtype=object
        )
        self.circle_hull_d = 5e3

        ''' READOUT RESONATORS '''
        self.resonators: List[EMResonatorTL3QbitWormRLTailXmonFork] = [None] * self.qubits_n

        ''' READOUT LINES '''
        self.ro_lines: List[DPathCPW] = [None, None]
        self.qCenter_roLine_distance = None

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
        self.draw_readout_resonators()
        # self.draw_readout_lines()

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

    def draw_qubits_array(self, new_disk_r=120e3):
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
            if qubit_idx in [0, 1, 2, 3, 4, 7, 10]:
                return 6
            else:
                return 2

        for qubit_idx in range(len(self.qubits_grid.pts_grid)):
            pt = self.qubits_grid.get_pt(qubit_idx)
            qubit_pars = QubitParams(
                squid_params=SQUID_PARS,
                qubit_cap_params=DiskConn8Pars(disk_r=new_disk_r),
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

    def draw_qq_couplings(self, donut_metal_width=40e3):
        # TODO: maybe transfer this datastructure to another file
        # incidence matrix for qubits graph
        # incidence matrix entries consists of 2 numbers - corresponding
        # qubits connectors idxs (see schematics for details)
        # if any of connectors idxs is equal to `-1` then qubit pair considered disconnected
        qq_coupling_connectors_map = np.zeros((12, 12, 2), dtype=int) - 1
        # TODO: fill structure automatically for more qubits
        # horizontal
        qq_coupling_connectors_map[0, 1] = np.array((0, 4))
        qq_coupling_connectors_map[1, 2] = np.array((0, 4))
        #
        qq_coupling_connectors_map[3, 4] = np.array((0, 4))
        qq_coupling_connectors_map[4, 5] = np.array((0, 4))
        qq_coupling_connectors_map[5, 6] = np.array((0, 4))
        #
        qq_coupling_connectors_map[7, 8] = np.array((0, 4))
        qq_coupling_connectors_map[8, 9] = np.array((0, 4))

        qq_coupling_connectors_map[10, 11] = np.array((0, 4))

        # vertical
        qq_coupling_connectors_map[0, 4] = np.array((2, 6))
        qq_coupling_connectors_map[1, 5] = np.array((2, 6))
        qq_coupling_connectors_map[2, 6] = np.array((2, 6))
        #
        qq_coupling_connectors_map[3, 7] = np.array((2, 6))
        qq_coupling_connectors_map[4, 8] = np.array((2, 6))
        qq_coupling_connectors_map[5, 9] = np.array((2, 6))
        #
        qq_coupling_connectors_map[7, 10] = np.array((2, 6))
        qq_coupling_connectors_map[8, 11] = np.array((2, 6))

        it_1d = list(enumerate(self.qubits_grid.pts_grid))
        it_2d = itertools.product(it_1d, it_1d)
        # TODO: refactor code:
        #  `it_2d` and `qq_coupling_connectors_map` has to be putted into the QubitsGrid
        #  datastructure
        for (pt1_1d_idx, pt1), (pt2_1d_idx, pt2) in it_2d:
            # cast (1,2) float arrays into points
            pt1 = DPoint(pt1[0], pt1[1])
            pt2 = DPoint(pt2[0], pt2[1])

            if pt1_1d_idx >= pt2_1d_idx:
                continue

            qq_coupling_connectors_idxs = qq_coupling_connectors_map[pt1_1d_idx, pt2_1d_idx]

            if all(
                    [
                        (pt1 - pt2).abs() == 1,  # Manhattan qubits-neighbours
                        (qq_coupling_connectors_idxs >= 0).all()  # if coupled
                    ]
            ):
                q1_connector_idx = qq_coupling_connectors_idxs[0]
                q2_connector_idx = qq_coupling_connectors_idxs[1]
                # qq_coupling = CqqCouplingType2(
                #     origin=DPoint(0, 0),
                #     params=CqqCouplingParamsType2(
                #         disk1=self.qubits[pt1_1d_idx].cap_shunt,
                #         disk2=self.qubits[pt2_1d_idx].cap_shunt,
                #         disk1_connector_idx=q1_connector_idx,
                #         disk2_connector_idx=q2_connector_idx
                #     ),
                #     region_id="ph"
                # )
                qq_coupling = CqqCouplingType1(
                    origin=DPoint(0, 0),
                    params=CqqCouplingParamsType1(
                        donut_metal_width=donut_metal_width,
                        disk1=self.qubits[pt1_1d_idx].cap_shunt,
                        disk2=self.qubits[pt2_1d_idx].cap_shunt,
                        disk1_connector_idx=q1_connector_idx,
                        disk2_connector_idx=q2_connector_idx
                    ),
                    region_id="ph"
                )
                qq_coupling.place(self.region_ph, region_id="ph")
                self.q_couplings[pt1_1d_idx, pt2_1d_idx] = qq_coupling

    def draw_readout_resonators(self):
        resonator_kw_args_list = list(
            map(
                ResonatorParams.get_resonator_params_by_qubit_idx, range(12)
            )
        )

        '''
        Resonators are placed at origin and then translated to their corresponding qubit.
        See 8QStair.drawio for details 
        '''

        # TODO: create separeta structure for this?
        qubit_res_finger_lengths_list = [30e3] * 12
        qubit_res_finger_gaps_list = [10e3] * 12
        qubit_res_finger_width_list = [45e3] * 12

        q_res_idxs_pairs = [[6, 0], [3, 1], [0, 2], [1, 3],
                            [2, 7], [5, 6], [4, 5], [7, 4]]
        resonator_rotation_angles = [90, 90, 90, 180,
                                     270, 270, 270 + 45, 0]

        for (q_idx, res_idx), res_trans_angle in zip(q_res_idxs_pairs, resonator_rotation_angles):
            qubit = self.qubits[q_idx]
            resonator_kw_args = resonator_kw_args_list[res_idx]
            trans_res_rotation = DCplxTrans(1, res_trans_angle, False, 0, 0)
            resonator_kw_args.update(
                {
                    "start"   : DPoint(0, 0),
                    "trans_in": trans_res_rotation
                }
            )
            res = EMResonatorTL3QbitWormRLTailXmonFork(**resonator_kw_args)

            # moving resonator to it's corresponding qubit
            qubit_res_finger_length = qubit_res_finger_lengths_list[q_idx]
            qubit_res_finger_gap = qubit_res_finger_gaps_list[q_idx]
            qubit_res_finger_width = qubit_res_finger_width_list[q_idx]
            qubit_res_d = 254e3
            dv = res.start - res.end + qubit.origin + \
                 trans_res_rotation * DVector(0, qubit_res_d)
            res.make_trans(DCplxTrans(1, 0, False, dv))
            res.place(self.region_ph)
            self.resonators[res_idx] = res

            # drawing qubit finger coupling
            dv_qubit_res = (res.fork_y_cpw1.end + res.fork_y_cpw2.end) / 2 - qubit.origin
            dv_qubit_res /= dv_qubit_res.abs()

            qubit_res_finger_bandage = CPW(
                start=qubit.origin,
                end=qubit.origin + (qubit_res_finger_length +
                                    qubit.qubit_params.qubit_cap_params.disk_r) * dv_qubit_res,
                width=qubit_res_finger_width, gap=0,
                open_end_gap=res.fork_gnd_gap
            )
            qubit_res_finger_bandage.place(self.region_ph)
            qubit_res_finger = CPW(
                start=qubit.origin + qubit.qubit_params.qubit_cap_params.disk_r * dv_qubit_res,
                end=qubit.origin + (qubit_res_finger_length +
                                    qubit.qubit_params.qubit_cap_params.disk_r) * dv_qubit_res,
                width=qubit_res_finger_width, gap=qubit_res_finger_gap,
                open_end_gap=res.fork_gnd_gap
            )
            qubit_res_finger.place(self.region_ph)

    def draw_readout_lines(self):
        # readout line is extended around qubit square in order to
        # fit readout resonators `L_couplings` and left a bit more space
        # for consistent and easy simulation of notch port resonator
        self.qCenter_roLine_distance = abs((self.qubits[6].origin - self.resonators[0].start).x) + \
                                       ResonatorParams.to_line_list[6]
        ro_line_extension = self.qCenter_roLine_distance / 2
        turn_radii = ro_line_extension / 4

        # left readout line
        p0_start = self.contact_pads[-1].end
        p0_end = self.contact_pads[8].end
        p1 = DPoint(p0_start.x, self.qubits[3].origin.y + ro_line_extension)
        p2 = DPoint(self.qubits[3].origin.x - self.qCenter_roLine_distance, p1.y)
        p3 = DPoint(p2.x, self.qubits[0].origin.y - self.qCenter_roLine_distance)
        p4 = DPoint(self.qubits[2].origin.x + ro_line_extension, p3.y)
        p5 = p4 + DVector(0, -self.qCenter_roLine_distance)
        p6 = DPoint(p0_end.x, p5.y)
        pts = [p0_start, p1, p2, p3, p4, p5, p6, p0_end]
        self.ro_lines[0] = DPathCPW(
            points=pts,
            cpw_parameters=[CPWParameters(width=20e3, gap=10e3)],
            turn_radii=[turn_radii],
            trans_in=None,
            region_id="ph"
        )
        self.ro_lines[0].place(self.region_ph, region_id="ph")

        # right readout line
        p1_start = self.contact_pads[-2].end
        p1_end = self.contact_pads[9].end
        p1 = p1_start + DVector(0, -1e6)
        p2 = DPoint(self.qubits[6].origin.x - ro_line_extension, p1.y)
        p3 = DPoint(p2.x, self.qubits[6].origin.y + self.qCenter_roLine_distance)
        p4 = DPoint(self.qubits[2].origin.x + self.qCenter_roLine_distance, p3.y)
        p5 = DPoint(p4.x, self.qubits[2].origin.y - ro_line_extension)
        p6 = DPoint(p1_end.x, p5.y)
        pts = [p1_start, p1, p2, p3, p4, p5, p6, p1_end]
        self.ro_lines[1] = DPathCPW(
            points=pts,
            cpw_parameters=[CPWParameters(width=20e3, gap=10e3)],
            turn_radii=[turn_radii],
            trans_in=None,
            region_id="ph"
        )
        self.ro_lines[1].place(self.region_ph, region_id="ph")

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
            self.region_el_protection
        )
        self.lv.zoom_fit()


def simulate_Cqq(q1_idx, q2_idx=None, resolution=(5e3, 5e3)):
    resolution_dx, resolution_dy = resolution
    dl_list = [0, 10e3, 20e3]
    # x_distance_dx_list = [0]
    for dl in dl_list:
        ''' DRAWING SECTION START '''
        design = Design8QStair("testScript")
        design.draw_chip()
        design.draw_qubits_array(new_disk_r=DiskConn8Pars().disk_r)
        design.draw_qq_couplings(donut_metal_width=CqqCouplingParamsType1().donut_metal_width + dl)

        design.show()
        design.lv.zoom_fit()
        design.layout.write(
            os.path.join(
                PROJECT_DIR, f"Cqq_{q1_idx}_{q2_idx}_"
                             f"{dl:.3f}_.gds"
            )
        )
        '''DRAWING SECTION END'''

        ''' SIMULATION SECTION START '''
        if q2_idx is None:  # 1 terminal calculation: self-capacitance calculation
            q1_reg = design.qubits[q1_idx].metal_regions["ph"]
            C1, _, _ = simulate_cij(
                design, env_reg=None, layer=design.layer_ph,
                subregs=[q1_reg],
                resolution=resolution, print_values=True
            )
        else:  # 2 terminal calculation: C1, C2 and C12 capacitances calculated
            q1_reg = design.qubits[q1_idx].metal_regions["ph"]
            q2_reg = design.qubits[q2_idx].metal_regions["ph"]
            C1, C12, C2 = simulate_cij(
                design, env_reg=None, layer=design.layer_ph,
                subregs=[q1_reg, q2_reg],
                resolution=resolution, print_values=True
            )
        ''' SIMULATION SECTION END '''

        '''SAVING REUSLTS SECTION START'''
        output_filepath = os.path.join(
            PROJECT_DIR,
            f"Xmon_Cqq_{q1_idx}_{q2_idx}_results.csv"
        )

        if q2_idx is None:
            save_sim_results(
                output_filepath=output_filepath,
                design=design,
                additional_pars={"C1, fF": C1}
            )
        else:
            save_sim_results(
                output_filepath=output_filepath,
                design=design,
                additional_pars={"C1, fF": C1, "C2, fF": C2, "C12, fF": C12}
            )


class Cqq_type2(ChipDesign):
    def __init__(self, cell_name):
        super().__init__(cell_name)
        self.disks_d = 1000e3

    def draw(self):
        self.region_ph.insert(
            pya.DBox(DPoint(-2e6, -1e6), DPoint(2e6, 1e6))
        )

        origin = DPoint(0, 0)
        # `q1.cap_shunt = None` if `postpone_drawing=True`
        q1 = Qubit(origin=origin + DVector(0, -self.disks_d / 2), postpone_drawing=False)
        q2 = Qubit(origin=origin + DVector(0, self.disks_d / 2), postpone_drawing=False)
        q1.place(self.region_ph, region_id="ph")
        q2.place(self.region_ph, region_id="ph")

        coupling = CqqCouplingType2(
            origin=origin,
            params=CqqCouplingParamsType2(
                disk1=q1.cap_shunt, disk2=q2.cap_shunt,
                disk1_connector_idx=3, disk2_connector_idx=5
            ),
            postpone_drawing=False, region_id="ph", region_ids=["ph"]
        )
        coupling.place(self.region_ph, region_id="ph")


if __name__ == "__main__":
    ''' draw and show design for manual design evaluation '''
    FABRICATION.OVERETCHING = 0.0e3
    design = Design8QStair("testScript")
    design.draw()
    design.show()
    # test = Cqq_type2("cellName")
    # test.draw()
    # test.show()

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
    # simulate_Cqq(q1_idx=5, q2_idx=6, resolution=(2e3, 2e3))

    ''' MD line C_qd for md1,..., md6 '''
    # for md_idx in [0,1]:
    #     for q_idx in range(2):
    #         simulate_md_Cg(md_idx=md_idx, q_idx=q_idx, resolution=(1e3, 1e3))

    ''' Resonators Q and f sim'''
    # simulate_resonators_f_and_Q(resolution=(2e3, 2e3))

    ''' Resonators Q and f when placed together'''
    # simulate_resonators_f_and_Q_together()
