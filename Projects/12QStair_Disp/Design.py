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
from classLib.coplanars import CPW, CPW2CPW, CPWParameters, DPathCPW
from classLib.chipDesign import ChipDesign
from classLib.marks import MarkBolgar
from classLib.contactPads import ContactPad
from classLib.helpers import fill_holes, split_polygons, extended_region
from classLib.helpers import simulate_cij, save_sim_results, rotate_around
from classLib.shapes import Donut
from classLib.resonators import EMResonatorTL3QbitWormRLTail

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
        self.draw_readout_lines()
        self.draw_microwave_drvie_lines()
        self.draw_flux_control_lines()

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
        p1_end = self.contact_pads[12].end
        p1 = p1_start + DVector(0, -0.5e6)
        p2 = DPoint(6.05e6, 10.2e6)
        p3 = DPoint(7.50e6, 10.2e6)
        p4 = DPoint(8.146e6, 9.6e6)
        p5 = p4 + 6*DVector(ro_line_extension, -ro_line_extension)
        p6 = p1_end + DVector(-0.5e6, 0)
        pts = [p1_start, p1, p2, p3, p4, p5, p6, p1_end]
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
        p2 = DPoint(7.475e6, 9.777e6)
        p3 = DPoint(7.27e6, 8.704e6)
        p_end = self.qubits[8].origin + md_origin_disp
        cpwrl_md = DPathCPW(
            points=[p_start, p1, p2, p3, p_end],
            cpw_parameters=[self.z_md1],
            turn_radii=[r_turn]
        )
        self.cpw_md_lines[8] = cpwrl_md

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
            self.qubits[11].cap_shunt.pars.disk_r + self.qubits[10].cap_shunt.pars.disk_gap
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
            self.qubits[8].cap_shunt.pars.disk_r + self.qubits[10].cap_shunt.pars.disk_gap
        ) + DVector(8.0169e3, 0)
        p2 = DPoint(7.2e6, 9.8e6)
        p3 = DPoint(7.1e6, 8.7e6)
        p_tr_start = DPoint(p_end.x, 7.90e6)
        fl_dpath = DPathCPW(
            points=[p_start, p1, p2, p3, p_tr_start, p_end],
            cpw_parameters=[self.z_fl1],
            turn_radii=[r_turn]
        )
        self.cpw_fl_lines[8] = fl_dpath

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
