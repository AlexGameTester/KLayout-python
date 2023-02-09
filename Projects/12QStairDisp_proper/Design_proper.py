__version__ = "12QStair_0.0.0.1"

'''
Changes log

'''

# import built-ins
from typing import List
import os
import itertools
import copy

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

from pya import DCplxTrans, Trans, ICplxTrans

# import project lib
import classLib

reload(classLib)
from classLib.coplanars import CPW, CPW2CPW, CPWParameters, DPathCPW, Bridge1
from classLib.chipDesign import ChipDesign
from classLib.chipTemplates import CHIP_14x14_20pads
from classLib.marks import MarkBolgar
from classLib.contactPads import ContactPad
from classLib.helpers import fill_holes, split_polygons, extended_region
from classLib.helpers import simulate_cij, save_sim_results, rotate_around
from classLib.shapes import Donut
from classLib.resonators import EMResonatorTL3QbitWormRLTail
from classLib.josJ import AsymSquid
from classLib.testPads import TestStructurePadsSquare
from classLib.shapes import Rectangle

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
from designElementsGeometry import ROResonator, ROResonatorParams
from designElementsGeometry import CqrCouplingParamsType1

import Cqq_couplings

reload(Cqq_couplings)
from Cqq_couplings import CqqCouplingType2, CqqCouplingParamsType2
from Cqq_couplings import CqqCouplingType1, CqqCouplingParamsType1


class Design12QStair(ChipDesign):
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
        self.squids: List[AsymSquid] = [] * self.qubits_n
        ''' QUBIT COUPLINGS '''
        self.q_couplings: np.array = np.empty(
            (self.qubits_n, self.qubits_n),
            dtype=object
        )
        self.circle_hull_d = 5e3

        ''' READOUT RESONATORS '''
        self.resonators: List[EMResonatorTL3QbitWormRLTail] = [None] * self.qubits_n

        ''' READOUT LINES '''
        self.ro_lines: List[DPathCPW] = [None]*3
        self.qCenter_roLine_distance = None

        ''' Microwave control lines '''
        # Z = 49.0538 E_eff = 6.25103 (E = 11.45)
        self.z_md1: CPWParameters = CPWParameters(30e3, 15e3)
        self.z_md2: CPWParameters = CPWParameters(4e3, 4e3)
        self.cpw_md_lines: List[DPathCPW] = [None]*self.qubits_n
        # length of the smoothing part between normal thick and end-thin cpw for md line
        self.md_line_cpw12_smoothhing = 10e3

        # length ofa thin part of the microwave drive line end
        self.md_line_cpw2_len = 550e3

        ''' Flux control lines '''
        self.cpw_fl_lines: List[DPathCPW] = [None]*self.qubits_n
        # flux line widths at the end of flux line
        self.flux2ground_left_width = 2e3
        self.flux2ground_right_width = 4e3
        # Z = 49.0538 E_eff = 6.25103 (E = 11.45)
        self.z_fl1: CPWParameters = CPWParameters(30e3, 15e3)
        # Z = 50.136  E_eff = 6.28826 (E = 11.45)
        self.z_fl2: CPWParameters = CPWParameters(10e3, 5.7e3)

        ''' Test structures '''
        self.test_squids_pads: List[TestStructurePadsSquare] = []
        self.test_squids: List[AsymSquid] = []

        ''' Squid bandages '''
        self.bandages_regs_list: List = []
        self.bandage_width = 2.5e3 * np.sqrt(2)
        self.bandage_height = 5e3 * np.sqrt(2)
        self.bandage_r_outer = 2e3
        self.bandage_r_inner = 2e3
        self.bandage_curve_pts_n = 40

        ''' Litography alignment marks '''
        self.marks: List[MarkBolgar] = []

    def draw(self, design_params=None):
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
        self.draw_readout_lines()
        self.draw_microwave_drvie_lines()
        self.draw_flux_control_lines()

        self.draw_test_structures()
        self.draw_express_test_structures_pads()
        self.draw_bandages()

        self.add_chip_marking(text_bl=DPoint(7.5e6, 2e6), chip_name="8Q_0.0.0.1")

        self.draw_litography_alignment_marks()
        # self.draw_bridges()
        # self.draw_pinning_holes()
        # # 4Q_Disp_Xmon v.0.3.0.8 p.12 - ensure that contact pads has no holes
        # for contact_pad in self.contact_pads:
        #     contact_pad.place(self.region_ph)
        # self.extend_photo_overetching()
        # self.inverse_destination(self.region_ph)
        # # convert to gds acceptable polygons (without inner holes)
        # self.draw_cut_marks()
        # self.region_ph.merge()
        # self.resolve_holes()
        # # convert to litograph readable format. Litograph can't handle
        # # polygons with more than 200 vertices.
        # self.split_polygons_in_layers(max_pts=180)

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

    def draw_qubits_array(self, new_disk_r=DiskConn8Pars().disk_r):
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
            if qubit_idx in [0, 1, 2, 3, 4, 7]:
                return 5
            else:
                return 1

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

    def draw_qq_couplings(self, donut_metal_width=CqqCouplingParamsType1().donut_metal_width):
        # TODO: maybe transfer this datastructure to another file
        # incidence matrix for qubits graph
        # incidence matrix entries consists of 2 numbers - corresponding
        # qubits connectors idxs (see schematics for details)
        # if any of connectors idxs is equal to `-1` then qubit pair considered disconnected
        qq_coupling_connectors_map = np.zeros((12, 12, 2), dtype=int) - 1
        # TODO: fill structure automatically for more qubits
        # horizontal
        qq_coupling_connectors_map[0, 1] = np.array((0, 3))
        qq_coupling_connectors_map[1, 2] = np.array((0, 4))
        #
        qq_coupling_connectors_map[3, 4] = np.array((0, 3))
        qq_coupling_connectors_map[4, 5] = np.array((0, 4))
        qq_coupling_connectors_map[5, 6] = np.array((7, 4))
        #
        qq_coupling_connectors_map[7, 8] = np.array((0, 4))
        qq_coupling_connectors_map[8, 9] = np.array((7, 4))

        qq_coupling_connectors_map[10, 11] = np.array((0, 4))

        # vertical
        qq_coupling_connectors_map[0, 4] = np.array((2, 7))
        qq_coupling_connectors_map[1, 5] = np.array((2, 6))
        qq_coupling_connectors_map[2, 6] = np.array((2, 6))
        #
        qq_coupling_connectors_map[3, 7] = np.array((2, 7))
        qq_coupling_connectors_map[4, 8] = np.array((2, 6))
        qq_coupling_connectors_map[5, 9] = np.array((3, 6))
        #
        qq_coupling_connectors_map[7, 10] = np.array((2, 7))
        qq_coupling_connectors_map[8, 11] = np.array((3, 6))

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
                ROResonatorParams.get_resonator_params_by_qubit_idx, range(12)
            )
        )

        '''
        Resonators are placed at origin and then translated to their corresponding qubit.
        '''

        q_res_connector_idxs_pairs = [
            [10, 0, 4], [7, 1, 4], [3, 2, 4], [4, 3, 4],
            [11, 3, 0], [8, 2, 0], [9, 1, 0], [5, 0, 0],
            [6, 3, 0], [2, 2, 0], [1, 1, 4], [0, 0, 4]
        ]
        resonator_rotation_angles = [90, 90, 90, 90 + 45,
                                     0, -45, -45, -45,
                                     270, 270, 180, 180]

        for (q_idx, res_idx, q_res_connector_idx), res_trans_angle in zip(
                q_res_connector_idxs_pairs,
                resonator_rotation_angles
        ):
            qubit = self.qubits[q_idx]
            resonator_kw_args = resonator_kw_args_list[res_idx]
            trans_res_rotation = DCplxTrans(1, res_trans_angle, False, 0, 0)
            resonator_kw_args.update(
                {
                    "start": DPoint(0, 0),
                    "trans_in": trans_res_rotation,
                }
            )
            res = EMResonatorTL3QbitWormRLTail(**resonator_kw_args)

            # moving resonator to it's corresponding qubit
            qubit_res_d = 500e3
            dv = res.start - res.end + qubit.origin + \
                 trans_res_rotation * DVector(0, qubit_res_d)
            res.make_trans(DCplxTrans(1, 0, False, dv))

            # move resonators a bit more
            if q_idx in [0, 1, 4]:
                res.make_trans(
                    DCplxTrans(
                        1, 0, False,
                        DVector(
                            -CqrCouplingParamsType1().bendings_disk_center_d,
                            0
                        )
                    )
                )
                if q_idx == 4:
                    res.make_trans(
                        DCplxTrans(
                            1, 0, False,
                            DVector(
                                -self.qubits_grid.dx/2,
                                -self.qubits_grid.dy/2
                            )
                        )
                    )
            elif q_idx in [5, 8, 9]:
                res.make_trans(
                    DCplxTrans(
                        1, 0, False,
                        DVector(
                            CqrCouplingParamsType1().bendings_disk_center_d,
                            0
                        )
                    )
                )
                if q_idx in [5, 8]:
                    res.make_trans(
                        DCplxTrans(
                            1, 0, False,
                            DVector(
                                self.qubits_grid.dx/2,
                                self.qubits_grid.dy/2
                            )
                        )
                    )
                    if q_idx == 8:
                        res.make_trans(
                            DCplxTrans(
                                1, 0, False,
                                DVector(
                                    -self.qubits_grid.dx / 5,
                                    self.qubits_grid.dy / 5
                                )
                            )
                        )
                    if q_idx == 5:
                        res.make_trans(
                            DCplxTrans(
                                1, 0, False,
                                DVector(
                                    self.qubits_grid.dx / 5,
                                    -self.qubits_grid.dy / 5
                                )
                            )
                        )
            elif q_idx in [11]:
                res.make_trans(
                    DCplxTrans(
                        1, 0, False,
                        DVector(
                            CqrCouplingParamsType1().bendings_disk_center_d,
                            0
                        )
                    )
                )
            res.place(self.region_ph)
            self.resonators[res_idx] = res

            # draw res-q coupling
            coupling_pars = CqrCouplingParamsType1(
                disk1=qubit.cap_shunt, disk1_connector_idx=q_res_connector_idx,

            )
            angle1 = coupling_pars.disk1.angle_connections[
                         coupling_pars.disk1_connector_idx
                     ] / (2 * np.pi) * 360

            arc_coupler = Donut(
                origin=coupling_pars.disk1.origin,
                inner_r=coupling_pars.donut_disk_d + coupling_pars.disk1.pars.disk_r,
                ring_width=coupling_pars.donut_metal_width,
                alpha_start=-coupling_pars.donut_delta_alpha_deg / 2,
                alpha_end=coupling_pars.donut_delta_alpha_deg / 2,
            )
            rotate_around(arc_coupler, arc_coupler.origin, angle1)

            d_alpha_deg = 360 / 2 / np.pi * coupling_pars.donut_gnd_gap / (
                    (arc_coupler.inner_r +
                     arc_coupler.outer_r)
                    / 2)
            arc_coupler_empty = Donut(
                origin=arc_coupler.origin,
                inner_r=arc_coupler.inner_r - coupling_pars.donut_disk_d,
                outer_r=arc_coupler.outer_r + coupling_pars.donut_gnd_gap,
                alpha_start=arc_coupler.alpha_start - d_alpha_deg,
                alpha_end=arc_coupler.alpha_end + d_alpha_deg,
                inverse=True
            )
            rotate_around(arc_coupler_empty, arc_coupler_empty.origin, angle1)
            # first clear metal, then fill with donut
            arc_coupler_empty.place(self.region_ph)
            arc_coupler.place(self.region_ph)

            # empty first, filling metal later

            resonator_end = res.end
            connector_dv_n = DVector(
                np.cos(2 * np.pi * (angle1 / 360)), np.sin(2 * np.pi * (angle1 / 360))
                )
            disk_far_bending_point = coupling_pars.disk1.origin + \
                                     coupling_pars.bendings_disk_center_d * connector_dv_n
            res_end_dv_n = DVector(np.cos(res.alpha_end), np.sin(res.alpha_end))
            resonator_first_bending = resonator_end + resonator_kw_args["r"] * res_end_dv_n
            res_donut_cpw_path = DPathCPW(
                points=[resonator_end, resonator_first_bending, disk_far_bending_point,
                        arc_coupler.outer_arc_center],
                cpw_parameters=[resonator_kw_args["Z0"]],
                turn_radii=[resonator_kw_args["r"]]
            )
            res_donut_cpw_path.place(self.region_ph)

            coupling_bandage = CPW(
                start=res_donut_cpw_path.end,
                end=(arc_coupler.inner_arc_center + arc_coupler.outer_arc_center) / 2,
                width=resonator_kw_args["Z0"].width,
                gap=0
            )
            coupling_bandage.place(self.region_ph)

    def draw_readout_lines(self):
        # readout line is extended around qubit square in order to
        # fit readout resonators `L_couplings` and left a bit more space
        # for consistent and easy simulation of notch port resonator
        self.qCenter_roLine_distance = abs((self.qubits[10].origin - self.resonators[0].start).x)\
                                       + \
                                       ResonatorParams.to_line_list[10]
        ro_line_extension = self.qCenter_roLine_distance / 2
        turn_radii = ro_line_extension / 4

        # readout line 1
        p0_start = self.contact_pads[0].end
        p0_end = self.contact_pads[7].end
        p1 = p0_start + DPoint(1e6, 0)
        p2 = DPoint(3.7e6, 9.1e6)
        p3 = DPoint(p2.x, 6.0e6)
        p4 = DPoint(4.13e6, 5.134e6)
        p5 = p4 + 2*DVector(ro_line_extension, -ro_line_extension)
        p7 = p0_end + DPoint(0, 0.5e6)
        p6 = DPoint(p5.x, p7.y)
        pts = [p0_start, p1, p2, p3, p4, p5, p6, p7, p0_end]
        self.ro_lines[0] = DPathCPW(
            points=pts,
            cpw_parameters=[CPWParameters(width=20e3, gap=10e3)],
            turn_radii=[turn_radii],
            trans_in=None,
            region_id="ph"
        )
        self.ro_lines[0].place(self.region_ph, region_id="ph")

        # readout line 2
        p1_start = self.contact_pads[15].end
        p1_end = self.contact_pads[8].end
        p1 = p1_start + DVector(0, -0.5e6)
        p2 = DPoint(6.05e6, 10.2e6)
        p3 = DPoint(7.50e6, 10.2e6)
        p4 = DPoint(8.146e6, 9.6e6)
        p5 = p4 + 6*DVector(ro_line_extension, -ro_line_extension)
        p7 = p1_end + DVector(0, 0.5e6)
        p6 = DPoint(p5.x, p6.y)
        pts = [p1_start, p1, p2, p3, p4, p5, p6, p7, p1_end]
        self.ro_lines[1] = DPathCPW(
            points=pts,
            cpw_parameters=[CPWParameters(width=20e3, gap=10e3)],
            turn_radii=[turn_radii],
            trans_in=None,
            region_id="ph"
        )
        self.ro_lines[1].place(self.region_ph, region_id="ph")

    def draw_microwave_drvie_lines(self):
        r_turn = 100e3
        ''' for qubit group №1 '''
        # place caplanar line 0md
        trans_len = 50e3

        md_origin_x_disp = 200e3
        md_origin_y_disp = 200e3
        md_origin_disp = DVector(-md_origin_x_disp, md_origin_y_disp)
        p_start = self.contact_pads[19].end
        p1 = p_start + DVector(0, -0.25e6)
        p2 = DPoint(3.1e6, 12.1e6)
        p_end = self.qubits[10].origin + md_origin_disp
        cpwrl_md = DPathCPW(
            points=[p_start, p1, p2, p_end],
            cpw_parameters=[self.z_md1],
            turn_radii=[r_turn]
        )
        self.cpw_md_lines[10] = cpwrl_md

        p_start = self.contact_pads[17].end
        p1 = p_start + DVector(0, -0.25e6)
        p2 = DPoint(5.8e6, 12.1e6)
        p_end = self.qubits[11].origin + md_origin_disp
        cpwrl_md = DPathCPW(
            points=[p_start, p1, p2, p_end],
            cpw_parameters=[self.z_md1],
            turn_radii=[r_turn]
        )
        self.cpw_md_lines[11] = cpwrl_md

        md_origin_disp = DVector(md_origin_x_disp, md_origin_y_disp)
        p_start = self.contact_pads[13].end
        p1 = p_start + DVector(-0.25e6, 0)
        p2 = DPoint(7.934e6, 10.05e6)
        p3 = DPoint(7.475e6, 9.777e6)
        p4 = DPoint(7.27e6, 8.704e6)
        p_end = self.qubits[8].origin + md_origin_disp
        cpwrl_md = DPathCPW(
            points=[p_start, p1, p2, p3, p4, p_end],
            cpw_parameters=[self.z_md1],
            turn_radii=[r_turn]
        )
        self.cpw_md_lines[8] = cpwrl_md

        md_origin_disp = DVector(0, md_origin_y_disp)
        p_start = self.contact_pads[12].end
        p1 = p_start + DVector(-0.25e6, 0)
        p2 = DPoint(12.51e6, 9.15e6)
        p3 = DPoint(8.650e6, 9.254e6)
        p4 = DPoint(7.9731e6, 8.571e6)
        p_end = self.qubits[9].origin + md_origin_disp
        cpwrl_md = DPathCPW(
            points=[p_start, p1, p2, p3, p4, p_end],
            cpw_parameters=[self.z_md1],
            turn_radii=[r_turn]
        )
        self.cpw_md_lines[9] = cpwrl_md

        md_origin_disp = DVector(0, md_origin_y_disp)
        p_start = self.contact_pads[10].end
        p1 = p_start + DVector(-0.25e6, 0)
        p2 = DPoint(11.4e6, 8.46e6)
        p3 = DPoint(9.351e6, 8.577e6)
        p4 = DPoint(8.092e6, 7.394e6)
        p_end = self.qubits[5].origin + md_origin_disp
        cpwrl_md = DPathCPW(
            points=[p_start, p1, p2, p3, p4, p_end],
            cpw_parameters=[self.z_md1],
            turn_radii=[r_turn]
        )
        self.cpw_md_lines[5] = cpwrl_md

        q_idx = 7
        md_origin_disp = DVector(0, -md_origin_y_disp)
        p_start = self.contact_pads[2].end
        p1 = p_start + DVector(0.25e6, 0)
        p2 = DPoint(3.426e6, 6.905e6)
        p3 = DPoint(5.010e6, 6.905e6)
        p_end = self.qubits[q_idx].origin + md_origin_disp
        cpwrl_md = DPathCPW(
            points=[p_start, p1, p2, p3, p_end],
            cpw_parameters=[self.z_md1],
            turn_radii=[r_turn]
        )
        self.cpw_md_lines[q_idx] = cpwrl_md

        q_idx = 3
        md_origin_disp = DVector(0, -md_origin_y_disp)
        p_start = self.contact_pads[4].end
        p1 = p_start + DVector(0.25e6, 0)
        p2 = DPoint(3.87e6, 5.21e6)
        p3 = DPoint(4.89e6, 5.79e6)
        p_end = self.qubits[q_idx].origin + md_origin_disp
        cpwrl_md = DPathCPW(
            points=[p_start, p1, p2, p3, p_end],
            cpw_parameters=[self.z_md1],
            turn_radii=[r_turn]
        )
        self.cpw_md_lines[q_idx] = cpwrl_md

        q_idx = 4
        md_origin_disp = DVector(0, -md_origin_y_disp)
        p_start = self.contact_pads[6].end
        p1 = p_start + DVector(0, 0.25e6)
        p2 = DPoint(4.82e6, 4.23e6)
        p3 = DPoint(5.84e6, 5.23e6)
        p_end = self.qubits[q_idx].origin + md_origin_disp
        cpwrl_md = DPathCPW(
            points=[p_start, p1, p2, p3, p_end],
            cpw_parameters=[self.z_md1],
            turn_radii=[r_turn]
        )
        self.cpw_md_lines[q_idx] = cpwrl_md

        for i, cpw_md_line in enumerate(self.cpw_md_lines):
            if cpw_md_line is not None:
                self.modify_md_line_end_and_place(
                    cpw_md_line, mod_length=self.md_line_cpw2_len, smoothing=self.md_line_cpw12_smoothhing
                )

    def modify_md_line_end_and_place(self, md_line: DPathCPW,
                                     mod_length=100e3, smoothing=20e3):
        """
        Changes coplanar for `mod_length` length from the end of `md_line`.
        Transition region length along the `md_line` is controlled by passing `smoothing` value.

        Notes
        -------
        Unhandled behaviour if transition (smoothing) region overlaps with any number of `CPWArc`s.

        Parameters
        ----------
        md_line : DPathCPW
            line object to modify
        mod_length : float
            length counting from the end of line to be modified
        smoothing : float
            length of smoothing for CPW2CPW transition between normal-thick and end-thin parts

        Returns
        -------
        None
        """
        # make flux line wider to have a less chance to misalign
        # bandage-eBeam-photo layers at the qubit bandage region.
        last_lines = {}
        total_length = 0
        for key, primitive in reversed(list(md_line.primitives.items())):
            last_lines[key] = (primitive)
            total_length += primitive.length()
            if total_length > mod_length:
                break

        # iterate from the end of the line
        # TODO: this is not working, objects already drawn in their corresponding coordinate system
        for key, primitive in list(last_lines.items())[:-1]:
            primitive.width = self.z_md2.width
            primitive.gap = self.z_md2.gap

        transition_line = list(last_lines.values())[-1]
        length_to_mod_left = mod_length - sum(
            [primitive.length() for primitive in
             list(last_lines.values())[:-1]])
        # divide line into 3 sections with proportions `alpha_i`
        beta_2 = smoothing
        beta_3 = length_to_mod_left
        beta_1 = transition_line.length() - (beta_2 + beta_3)
        transition_line_dv = transition_line.end - transition_line.start
        # tangent vector
        transition_line_dv_s = transition_line_dv / transition_line_dv.abs()
        p1 = transition_line.start
        p2 = p1 + beta_1 * transition_line_dv_s
        p3 = p2 + beta_2 * transition_line_dv_s
        p4 = p3 + beta_3 * transition_line_dv_s
        cpw_transition_line1 = CPW(
            start=p1,
            end=p2,
            cpw_params=self.z_md1
        )
        cpw_transition = CPW2CPW(
            Z0=self.z_md1,
            Z1=self.z_md2,
            start=p2, end=p3
        )
        cpw_transition_line2 = CPW(
            start=p3 - 2 * transition_line_dv_s,
            # rounding error correction
            end=p4,
            cpw_params=self.z_md2
        )

        transition_line_name = list(last_lines.keys())[-1]
        new_primitives = {}
        for key, primitive in md_line.primitives.items():
            if key != transition_line_name:
                new_primitives[key] = primitive
            elif key == transition_line_name:
                new_primitives[
                    transition_line_name + "_1"] = cpw_transition_line1
                new_primitives[
                    transition_line_name + "_2"] = cpw_transition
                new_primitives[
                    transition_line_name + "_3"] = cpw_transition_line2

        # make open-circuit end for md line end
        cpw_end = list(last_lines.values())[0]
        end_dv = cpw_end.end - cpw_end.start
        end_dv_s = end_dv / end_dv.abs()
        p1 = cpw_end.end
        p2 = p1 + cpw_end.b / 2 * end_dv_s
        md_open_end = CPW(
            start=p1, end=p2,
            width=0, gap=cpw_transition_line2.b / 2
        )
        new_primitives["md_open_end"] = md_open_end
        md_line.primitives = new_primitives
        md_line.connections[1] = list(md_line.primitives.values())[-1].end
        md_line._refresh_named_connections()
        md_line.place(self.region_ph)

    def draw_flux_control_lines(self):
        r_turn = 100e3

        p_start = self.contact_pads[18].end
        p1 = p_start + DVector(0, -0.25e6)
        p_end = self.qubits[10].origin + DVector(
            0,
            self.qubits[10].cap_shunt.pars.disk_r + self.qubits[10].cap_shunt.pars.disk_gap
        ) + DVector(8.0169e3, 0)
        p_tr_start = DPoint(p_end.x, 9.7e6)
        fl_dpath = DPathCPW(
            points=[p_start, p1, p_tr_start, p_end],
            cpw_parameters=[self.z_fl1],
            turn_radii=[r_turn]
        )
        self.cpw_fl_lines[10] = fl_dpath

        p_start = self.contact_pads[16].end
        p1 = p_start + DVector(0, -0.25e6)
        p_end = self.qubits[11].origin + DVector(
            0,
            self.qubits[11].cap_shunt.pars.disk_r + self.qubits[11].cap_shunt.pars.disk_gap
        ) + DVector(8.0169e3, 0)
        p2 = DPoint(6.374e6, 11.471e6)
        p_tr_start = DPoint(p_end.x, 9.7e6)
        fl_dpath = DPathCPW(
            points=[p_start, p1, p2, p_tr_start, p_end],
            cpw_parameters=[self.z_fl1],
            turn_radii=[r_turn]
        )
        self.cpw_fl_lines[11] = fl_dpath

        p_start = self.contact_pads[14].end
        p1 = p_start + DVector(-0.25e6, 0)
        p_end = self.qubits[8].origin + DVector(
            0,
            self.qubits[8].cap_shunt.pars.disk_r + self.qubits[8].cap_shunt.pars.disk_gap
        ) + DVector(8.0169e3, 0)
        p2 = DPoint(7.703e6, 10.215e6)
        p3 = DPoint(7.2e6, 9.8e6)
        p4 = DPoint(7.1e6, 8.7e6)
        p_tr_start = DPoint(p_end.x, 7.90e6)
        fl_dpath = DPathCPW(
            points=[p_start, p1, p2, p3, p4, p_tr_start, p_end],
            cpw_parameters=[self.z_fl1],
            turn_radii=[r_turn]
        )
        self.cpw_fl_lines[8] = fl_dpath

        p_start = self.contact_pads[11].end
        p1 = p_start + DVector(-0.35e6, 0)
        p_end = self.qubits[9].origin + 1/np.sqrt(2)*DVector(
            self.qubits[9].cap_shunt.pars.disk_r + self.qubits[9].cap_shunt.pars.disk_gap,
            self.qubits[9].cap_shunt.pars.disk_r + self.qubits[9].cap_shunt.pars.disk_gap
        ) + 1/np.sqrt(2)*DVector(8.0169e3, -8.0169e3)
        p2 = DPoint(11.6e6, 8.92e6)
        p3 = DPoint(8.842e6, 9.081e6)
        p4 = DPoint(7.902e6, 8.177e6)
        p5 = DPoint(7.968e6, 7.925e6)
        p_tr_start = p_end + \
                     1/np.sqrt(2)*DVector(
            CqqCouplingParamsType1().bendings_disk_center_d,
            CqqCouplingParamsType1().bendings_disk_center_d
        )
        fl_dpath = DPathCPW(
            points=[p_start, p1, p2, p3, p4, p5, p_tr_start, p_end],
            cpw_parameters=[self.z_fl1],
            turn_radii=[r_turn]
        )
        self.cpw_fl_lines[9] = fl_dpath

        p_start = self.contact_pads[9].end
        p1 = p_start + DVector(0, 0.25e6)
        p2 = p1 + DVector(-0.25e6, 0)
        p_end = self.qubits[5].origin + 1 / np.sqrt(2) * DVector(
            self.qubits[5].cap_shunt.pars.disk_r + self.qubits[5].cap_shunt.pars.disk_gap,
            self.qubits[5].cap_shunt.pars.disk_r + self.qubits[5].cap_shunt.pars.disk_gap
        ) + 1 / np.sqrt(2) * DVector(8.0169e3, -8.0169e3)
        p3 = DPoint(11.1e6, 8.31e6)
        p4 = DPoint(9.525e6, 8.414e6)
        p5 = DPoint(8.632e6, 7.466e6)
        p_tr_start = p_end + \
                     1 / np.sqrt(2) * DVector(
            CqqCouplingParamsType1().bendings_disk_center_d,
            CqqCouplingParamsType1().bendings_disk_center_d
        )
        fl_dpath = DPathCPW(
            points=[p_start, p1, p2, p3, p4, p5, p_tr_start, p_end],
            cpw_parameters=[self.z_fl1],
            turn_radii=[r_turn]
        )
        self.cpw_fl_lines[5] = fl_dpath

        q_idx = 7
        p_start = self.contact_pads[1].end
        p1 = p_start + DVector(0.25e6, 0)
        p_end = self.qubits[q_idx].origin - 1 / np.sqrt(2) * DVector(
            self.qubits[q_idx].cap_shunt.pars.disk_r + self.qubits[q_idx].cap_shunt.pars.disk_gap,
            self.qubits[q_idx].cap_shunt.pars.disk_r + self.qubits[q_idx].cap_shunt.pars.disk_gap
        ) + 1 / np.sqrt(2) * DVector(-8.0169e3, 8.0169e3)
        p2 = DPoint(3.426e6, 7.137e6)
        p3 = DPoint(4.010e6, 7.137e6)
        p_tr_start = p_end - \
                     1 / np.sqrt(2) * DVector(
            CqqCouplingParamsType1().bendings_disk_center_d,
            CqqCouplingParamsType1().bendings_disk_center_d
        )
        fl_dpath = DPathCPW(
            points=[p_start, p1, p2, p3, p_tr_start, p_end],
            cpw_parameters=[self.z_fl1],
            turn_radii=[r_turn]
        )
        self.cpw_fl_lines[q_idx] = fl_dpath

        q_idx = 3
        p_start = self.contact_pads[3].end
        p1 = p_start + DVector(0.25e6, 0)
        p_end = self.qubits[q_idx].origin - 1 / np.sqrt(2) * DVector(
            self.qubits[q_idx].cap_shunt.pars.disk_r + self.qubits[q_idx].cap_shunt.pars.disk_gap,
            self.qubits[q_idx].cap_shunt.pars.disk_r + self.qubits[q_idx].cap_shunt.pars.disk_gap
        ) + 1 / np.sqrt(2) * DVector(-8.0169e3, 8.0169e3)
        p2 = DPoint(3.72e6, 5.68e6)
        p3 = DPoint(4.56e6, 5.92e6)
        p_tr_start = p_end - \
                     1 / np.sqrt(2) * DVector(
            CqqCouplingParamsType1().bendings_disk_center_d,
            CqqCouplingParamsType1().bendings_disk_center_d
        )
        fl_dpath = DPathCPW(
            points=[p_start, p1, p2, p3, p_tr_start, p_end],
            cpw_parameters=[self.z_fl1],
            turn_radii=[r_turn]
        )
        self.cpw_fl_lines[q_idx] = fl_dpath

        q_idx = 4
        p_start = self.contact_pads[5].end
        p1 = p_start + DVector(0, 0.25e6)
        p_end = self.qubits[q_idx].origin - 1 / np.sqrt(2) * DVector(
            self.qubits[q_idx].cap_shunt.pars.disk_r + self.qubits[q_idx].cap_shunt.pars.disk_gap,
            self.qubits[q_idx].cap_shunt.pars.disk_r + self.qubits[q_idx].cap_shunt.pars.disk_gap
        ) + 1 / np.sqrt(2) * DVector(-8.0169e3, 8.0169e3)
        p2 = DPoint(4.70e6, 4.43e6)
        p3 = DPoint(5.57e6, 5.36e6)
        p_tr_start = p_end - \
                     1 / np.sqrt(2) * DVector(
            CqqCouplingParamsType1().bendings_disk_center_d,
            CqqCouplingParamsType1().bendings_disk_center_d
        )
        fl_dpath = DPathCPW(
            points=[p_start, p1, p2, p3, p_tr_start, p_end],
            cpw_parameters=[self.z_fl1],
            turn_radii=[r_turn]
        )
        self.cpw_fl_lines[q_idx] = fl_dpath

        for flux_line in self.cpw_fl_lines:
            if flux_line is not None:
                self.modify_flux_line_end_and_place(flux_line)
        self.region_ph.merge()

    def modify_flux_line_end_and_place(self, flux_line: DPathCPW):
        # make flux line wider to have a less chance to misalign
        # bandage-eBeam-photo layers at the qubit bandage region.
        last_line = list(flux_line.primitives.values())[-1]
        last_line_name = list(flux_line.primitives.keys())[-1]

        # divide line into 3 sections with proportions `alpha_i`
        alpha_1 = 0.6
        alpha_2 = 0.3
        alpha_3 = 1 - (alpha_1 + alpha_2)
        last_line_dv = last_line.end - last_line.start
        last_line_dv_s = last_line_dv / last_line_dv.abs()
        p1 = last_line.start
        p2 = p1 + alpha_1 * last_line_dv
        p3 = p2 + alpha_2 * last_line_dv
        p4 = p3 + alpha_3 * last_line_dv
        cpw_normal = CPW(
            start=p1,
            end=p2,
            width=self.z_fl1.width,
            gap=self.z_fl1.gap
        )
        cpw_transition = CPW2CPW(
            Z0=self.z_fl1,
            Z1=self.z_fl2,
            start=p2, end=p3
        )
        cpw_thick = CPW(
            start=p3 - 2 * last_line_dv_s,  # rounding error correction
            end=p4,
            width=self.z_fl2.width,
            gap=self.z_fl2.gap
        )

        del flux_line.primitives[last_line_name]
        flux_line.primitives[last_line_name + "_1"] = cpw_normal
        flux_line.primitives[last_line_name + "_2"] = cpw_transition
        flux_line.primitives[last_line_name + "_3"] = cpw_thick

        # tangent vector
        last_line_dv_s = last_line_dv / last_line_dv.abs()
        # normal vector (90deg counterclockwise rotated tangent vector)
        last_line_dv_n = DVector(-last_line_dv_s.y, last_line_dv_s.x)

        # draw mutual inductance for squid, shunted to ground to the right
        # (view from flux line start towards squid)
        p1 = cpw_thick.end - cpw_thick.width / 2 * last_line_dv_n - \
             self.flux2ground_right_width / 2 * last_line_dv_s
        p2 = p1 - cpw_thick.gap * last_line_dv_n
        inductance_cpw = CPW(
            start=p1, end=p2,
            width=self.flux2ground_right_width, gap=0
        )
        flux_line.primitives["fl_end_ind_right"] = inductance_cpw

        # connect central conductor to ground to the left (view from
        # flux line start towards squid)
        p1 = cpw_thick.end + cpw_thick.width / 2 * last_line_dv_n - \
             self.flux2ground_left_width / 2 * last_line_dv_s
        p2 = p1 + cpw_thick.gap * last_line_dv_n
        flux_gnd_cpw = CPW(
            start=p1, end=p2,
            width=self.flux2ground_left_width, gap=0
        )
        flux_line.primitives["fl_end_ind_left"] = flux_gnd_cpw
        flux_line.connections[1] = list(flux_line.primitives.values())[
            -3].end
        flux_line._refresh_named_connections()
        flux_line.place(self.region_ph)

    def draw_test_structures(self):
        struct_centers = [
            DPoint(1.8e6, 6.0e6),
            DPoint(11e6, 10.7e6),
            DPoint(8.5e6, 4e6)
            ]
        for struct_center in struct_centers:
            ## JJ test structures ##
            dx = SQUID_PARS.SQB_dx / 2 - SQUID_PARS.SQLBT_dx / 2

            # test structure with big critical current (#1)
            test_struct1 = TestStructurePadsSquare(
                struct_center,
                # gnd gap in test structure is now equal to
                # the same of first xmon cross, where polygon is placed
                gnd_gap=20e3,
                pads_gap=self.qubits[0].cap_shunt.pars.disk_gap
            )
            self.test_squids_pads.append(test_struct1)
            test_struct1.place(self.region_ph)

            text_reg = pya.TextGenerator.default_generator().text(
                "56 nA", 0.001, 25, False, 0, 0)
            text_bl = test_struct1.empty_rectangle.p1 - DVector(0, 20e3)
            text_reg.transform(
                ICplxTrans(1.0, 0, False, text_bl.x, text_bl.y))
            self.region_ph -= text_reg

            pars_local = copy.deepcopy(SQUID_PARS)
            pars_local.SQRBT_dx = 0
            pars_local.SQRBJJ_dy = 0
            pars_local.bot_wire_x = [-dx]

            squid_center = test_struct1.center
            test_jj = AsymSquid(
                squid_center + DVector(0, -8.0001234e3),
                pars_local
            )
            self.test_squids.append(test_jj)
            test_jj.place(self.region_el)

            # test structure with low critical current (#2)
            test_struct2 = TestStructurePadsSquare(
                struct_center + DPoint(0.3e6, 0),
                gnd_gap=20e3,
                pads_gap=self.qubits[0].cap_shunt.pars.disk_gap
            )
            self.test_squids_pads.append(test_struct2)
            test_struct2.place(self.region_ph)

            text_reg = pya.TextGenerator.default_generator().text(
                "11 nA", 0.001, 25, False, 0, 0)
            text_bl = test_struct2.empty_rectangle.p1 - DVector(0, 20e3)
            text_reg.transform(
                ICplxTrans(1.0, 0, False, text_bl.x, text_bl.y))
            self.region_ph -= text_reg

            pars_local = copy.deepcopy(SQUID_PARS)
            pars_local.SQLBT_dx = 0
            pars_local.SQLBJJ_dy = 0
            pars_local.bot_wire_x = [dx]

            squid_center = test_struct2.center
            test_jj = AsymSquid(
                squid_center + DVector(0, -8.0001234e3),
                pars_local
            )
            self.test_squids.append(test_jj)
            test_jj.place(self.region_el)

            # test structure for bridge DC contact (#3)
            test_struct3 = TestStructurePadsSquare(
                struct_center + DPoint(0.6e6, 0))
            test_struct3.place(self.region_ph)
            text_reg = pya.TextGenerator.default_generator().text(
                "3xBrg 100um", 0.001, 25, False, 0, 0
            )
            text_bl = test_struct3.empty_rectangle.p1 - DVector(0, 20e3)
            text_reg.transform(
                ICplxTrans(1.0, 0, False, text_bl.x, text_bl.y)
            )
            self.region_ph -= text_reg

            test_bridges = []
            for i in range(3):
                bridge = Bridge1(
                    test_struct3.center + DPoint(50e3 * (i - 1), 0),
                    gnd2gnd_dy=100e3,
                    gnd_touch_dx=20e3
                )
                test_bridges.append(bridge)
                bridge.place(self.region_bridges1, region_id="bridges_1")
                bridge.place(self.region_bridges2, region_id="bridges_2")

        # bandages test structures
        test_dc_el2_centers = [
            DPoint(1.8e6, 7.8e6),
            DPoint(12.5e6, 10.7e6),
            DPoint(8.5e6, 3e6)
        ]
        for struct_center in test_dc_el2_centers:
            test_struct1 = TestStructurePadsSquare(struct_center)
            test_struct1.place(self.region_ph)
            text_reg = pya.TextGenerator.default_generator().text(
                "Bandage", 0.001, 40, False, 0, 0)
            text_bl = test_struct1.empty_rectangle.origin + DPoint(
                test_struct1.gnd_gap, -4 * test_struct1.gnd_gap
            )
            text_reg.transform(
                ICplxTrans(1.0, 0, False, text_bl.x, text_bl.y))
            self.region_ph -= text_reg

            rec_width = 10e3
            rec_height = test_struct1.pads_gap + 2 * rec_width
            p1 = struct_center - DVector(rec_width / 2, rec_height / 2)
            dc_rec = Rectangle(p1, rec_width, rec_height)
            dc_rec.place(self.dc_bandage_reg)

    def draw_express_test_structures_pads(self):
        el_pad_height = 30e3
        el_pad_width = 40e3
        # contact wire cpw parameters
        c_cpw = CPWParameters(width=1e3, gap=0)
        for squid, test_pad in zip(
                self.test_squids,
                self.test_squids_pads
        ):
            if not squid.squid_params.SQLBJJ_dy == 0:  # only left JJ is present
                # test pad to the right
                p1 = DPoint(test_pad.top_rec.p2.x, test_pad.center.y)
                p2 = p1 + DVector(-el_pad_width, 0)
                express_test_cpw_right = CPW(
                    start=p1, end=p2,
                    width=el_pad_height, gap=0
                )
                self.region_ph -= express_test_cpw_right.metal_region.sized(10e3)
                express_test_cpw_right.place(self.region_el)

                # right connection
                p3 = squid.TCW.center()
                p4 = DPoint(express_test_cpw_right.center().x, p3.y)
                right_conn_cpw = CPW(
                    start=p3, end=p4,
                    cpw_params=c_cpw
                )
                right_conn_cpw.place(self.region_el)

                # test pad on the left
                p1 = DPoint(test_pad.top_rec.p1.x, test_pad.center.y)
                p2 = p1 + DVector(el_pad_width, 0)
                express_test_cpw_left = CPW(
                    start=p1, end=p2,
                    width=el_pad_height, gap=0
                )
                self.region_ph -= express_test_cpw_left.metal_region.sized(10e3)
                express_test_cpw_left.place(self.region_el)

                p3 = squid.SQLBT.center()
                p4 = DPoint(express_test_cpw_left.center().x, p3.y)
                left_conn_cpw = CPW(
                    start=p3, end=p4,
                    cpw_params=c_cpw
                )
                left_conn_cpw.place(self.region_el)

            elif squid.squid_params.SQLBJJ_dy == 0:  # only right leg is present
                # test pad expanded to the left
                p1 = DPoint(test_pad.top_rec.p1.x, test_pad.center.y)
                p2 = p1 + DVector(el_pad_width, 0)
                express_test_cpw_lleft = CPW(
                    start=p1, end=p2,
                    width=el_pad_height,
                    gap=0
                )
                express_test_cpw_lleft.place(self.region_el)
                self.region_ph -= express_test_cpw_lleft.metal_region.sized(20e3)

                p1 = squid.SQRBT.center()
                p2 = DPoint(express_test_cpw_lleft.center().x, p1.y)
                left_conn_cpw = CPW(
                    start=p1, end=p2,
                    cpw_params=c_cpw
                )
                left_conn_cpw.place(self.region_el)

                # test pad expanded to the right
                p1 = DPoint(test_pad.top_rec.p2.x, test_pad.center.y)
                p2 = p1 + DVector(-el_pad_width, 0)
                express_test_cpw_right = CPW(
                    start=p1, end=p2,
                    width=el_pad_height,
                    gap=0
                )
                express_test_cpw_right.place(self.region_el)
                self.region_ph -= express_test_cpw_right.metal_region.size(20e3)

                p1 = squid.TCW.center()
                p2 = DPoint(express_test_cpw_right.center().x, p1.y)
                right_conn_cpw = CPW(
                    start=p1, end=p2,
                    cpw_params=c_cpw
                )
                right_conn_cpw.place(self.region_el)
                # cut_reg = etc4.metal_region.dup()
                # cut_reg.transform(Trans(Vector(1e3,0))).size(self.photo_recess_d)
                # self.region_ph -= cut_reg
                # self.region_el -= cut_reg

    def draw_bandages(self):
        """
        Returns
        -------

        """
        from itertools import chain
        for squid in chain(
                self.squids,
                self.test_squids
        ):
            # dc contact pad has to be completely
            # inside union of both  e-beam and photo deposed
            # metal regions.
            # `self.dc_cont_clearance` represents minimum distance
            # from dc contact pad`s perimeter to the perimeter of the
            # e-beam and photo-deposed metal perimeter.
            self.bandages_regs_list += self._draw_squid_bandage(
                squid,
                shift2sq_center=0
            )
            # collect all bottom contacts

    def _draw_squid_bandage(self, squid: AsymSquid = None,
                           shift2sq_center=0):
        # squid direction from bottom to top
        squid_BT_dv = squid.TC.start - squid.TC.end
        squid_BT_dv_s = squid_BT_dv / squid_BT_dv.abs()  # normalized

        bandages_regs_list: List[Region] = []

        # top bandage
        top_bandage_reg = self._get_bandage_reg(
            center=squid.TC.start,
            shift=-shift2sq_center * squid_BT_dv_s
        )
        bandages_regs_list.append(top_bandage_reg)
        self.dc_bandage_reg += top_bandage_reg

        # bottom contacts
        for i, _ in enumerate(squid.squid_params.bot_wire_x):
            BC = squid.BC_list[i]
            bot_bandage_reg = self._get_bandage_reg(
                center=BC.end,
                shift=shift2sq_center * squid_BT_dv_s
            )
            bandages_regs_list.append(bot_bandage_reg)
            self.dc_bandage_reg += bot_bandage_reg

        self._draw_recess(squid=squid)

        return bandages_regs_list

    def _get_bandage_reg(self, center, shift: DVector = DVector(0, 0)):
        center += shift
        rect_lb = center + \
                  DVector(
                      -self.bandage_width / 2,
                      -self.bandage_height / 2
                  )
        bandage_reg = Rectangle(
            origin=rect_lb,
            width=self.bandage_width,
            height=self.bandage_height
        ).metal_region
        bandage_reg.round_corners(
            self.bandage_r_inner,
            self.bandage_r_outer,
            self.bandage_curve_pts_n
        )

        return bandage_reg

    def _draw_recess(self, squid):
        recess_reg = squid.TC.metal_region.dup().size(-1e3)
        self.region_ph -= recess_reg

        # bottom recess(es)
        for i, _ in enumerate(squid.squid_params.bot_wire_x):
                BC = squid.BC_list[i]
                recess_reg = BC.metal_region.dup().size(-1e3)
                self.region_ph -= recess_reg

    def draw_litography_alignment_marks(self):
        marks_centers = [
            DPoint(0.5e6, 13.5e6), DPoint(8.3e6, 13.5e6), DPoint(13.5e6, 13.5e6),
            DPoint(0.5e6, 8.3e6), DPoint(13.5e6, 8.3e6),
            DPoint(0.5e6, 0.5e6), DPoint(8.3e6, 0.5e6), DPoint(13.5e6, 0.5e6)
        ]
        for mark_center in marks_centers:
            self.marks.append(
                MarkBolgar(mark_center)
            )
            self.marks[-1].place(self.region_ph)
            self.marks[-1].place(self.region_bridges1)

    def draw_bridges(self):
        bridges_step = 130e3
        fl_bridges_step = 130e3

        # for readout resonators
        for resonator in self.resonators:
            for name, res_primitive in resonator.primitives.items():
                if "coil" in name:
                    subprimitives = res_primitive.primitives
                    for primitive_name, primitive in subprimitives.items():
                        # place bridges only at arcs of coils
                        # but not on linear segments
                        if "arc" in primitive_name:
                            Bridge1.bridgify_CPW(
                                primitive, bridges_step,
                                gnd2gnd_dy=70e3,
                                dest=self.region_bridges1,
                                dest2=self.region_bridges2
                            )
                    continue
                elif "fork" in name:  # skip fork primitives
                    continue
                else:
                    # bridgify everything else except "arc1"
                    # resonator.primitives["arc1"] is arc that connects
                    # L_coupling with long vertical line for
                    # `EMResonatorTL3QbitWormRLTailXmonFork`
                    if name == "arc1":
                        continue
                    Bridge1.bridgify_CPW(
                        res_primitive, bridges_step,
                        gnd2gnd_dy=70e3,
                        dest=self.region_bridges1,
                        dest2=self.region_bridges2
                    )

        # for contact wires
        for key, val in self.__dict__.items():
            if "cpwrl_md" in key:
                cpwrl_md = val
                Bridge1.bridgify_CPW(
                    cpwrl_md, bridges_step,
                    gnd2gnd_dy=100e3,
                    dest=self.region_bridges1, dest2=self.region_bridges2,
                    avoid_points=[squid.origin for squid in self.squids],
                    avoid_distances=900e3
                )
            elif "cpwrl_fl" in key:
                cpwrl_fl = val
                Bridge1.bridgify_CPW(
                    cpwrl_fl, fl_bridges_step,
                    gnd2gnd_dy=100e3,
                    dest=self.region_bridges1, dest2=self.region_bridges2,
                    avoid_points=[squid.origin for squid in self.squids],
                    avoid_distances=900e3
                )

        # close bridges for cpw_fl line
        for i, cpw_fl in enumerate(self.cpw_fl_lines):
            dy_list = [30e3, 100e3, 235e3, 365e3, 495e3, 625e3, 755e3]
            for dy in dy_list:
                if i < 4:
                    pass
                elif i >= 4:
                    dy = -dy
                bridge_center1 = cpw_fl.end + DVector(0, -dy)
                br = Bridge1(center=bridge_center1, gnd2gnd_dy=70e3,
                             trans_in=Trans.R90)
                br.place(dest=self.region_bridges1,
                         region_id="bridges_1")
                br.place(dest=self.region_bridges2,
                         region_id="bridges_2")

        for i, cpw_md in enumerate(self.cpw_md_lines):
            dy_list = [110e3, 240e3, 370e3, 500e3, 630e3]
            for dy in dy_list:
                if i < 4:
                    pass
                elif i >= 4:
                    dy = -dy
                bridge_center1 = cpw_md.end + DVector(0, -dy)
                br = Bridge1(center=bridge_center1, gnd2gnd_dy=70e3,
                             trans_in=Trans.R90)
                br.place(dest=self.region_bridges1,
                         region_id="bridges_1")
                br.place(dest=self.region_bridges2,
                         region_id="bridges_2")

        # for readout waveguides
        avoid_points = []
        avoid_distances = []
        for res in self.resonators:
            av_pt = res.primitives["coil0"].primitives["cop1"].center()
            avoid_points.append(av_pt)
            av_dist = res.L_coupling / 2 + res.r + res.Z0.b / 2
            avoid_distances.append(av_dist)

        Bridge1.bridgify_CPW(
            self.ro_lines[0], gnd2gnd_dy=100e3,
            bridges_step=bridges_step,
            dest=self.region_bridges1, dest2=self.region_bridges2,
            avoid_points=avoid_points, avoid_distances=avoid_distances
        )
        Bridge1.bridgify_CPW(
            self.ro_lines[1], gnd2gnd_dy=100e3,
            bridges_step=bridges_step,
            dest=self.region_bridges1, dest2=self.region_bridges2,
            avoid_points=avoid_points, avoid_distances=avoid_distances
        )

    def draw_pinning_holes(self):
        # points that select polygons of interest if they were clicked at)
        selection_pts = [
            Point(0.1e6, 0.1e6),
            (self.cpwrl_ro_line1.start + self.cpwrl_ro_line1.end) / 2,
            (self.cpwrl_ro_line2.start + self.cpwrl_ro_line2.end) / 2
        ]

        # creating region of small boxes (polygon selection requires
        # regions)
        dv = DVector(1e3, 1e3)
        selection_regions = [
            Region(pya.Box(pt, pt + dv)) for pt in selection_pts
        ]
        selection_region = Region()
        for reg in selection_regions:
            selection_region += reg

        other_polys_reg = self.region_ph.dup().select_not_interacting(
            selection_region
        )

        reg_to_fill = self.region_ph.dup().select_interacting(
            selection_region
        )
        filled_reg = fill_holes(reg_to_fill, d=40e3, width=15e3,
                                height=15e3)

        self.region_ph = filled_reg + other_polys_reg

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
        design = Design12QStair("testScript")
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


if __name__ == "__main__":
    ''' draw and show design for manual design evaluation '''
    FABRICATION.OVERETCHING = 0.0e3
    design = Design12QStair("testScript")
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
    # design = Design12QStair("testScript")
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
