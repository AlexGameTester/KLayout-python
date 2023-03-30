__version__ = "12QStair_0.0.0.0"

'''
NOTE:
Every array that regards resonators or qubits structures has first axis ordered
in ascending order to qubit idxs.
E.g. self.resonators[1] - resonator that belongs to qubit №1 (starting from 0)

Changes log
    
'''

# import built-ins
from typing import List
import os
import itertools
import copy

PROJECT_DIR = os.path.dirname(__file__)
import sys
import shutil

sys.path.append(PROJECT_DIR)
from importlib import reload

# import good 3rd party
import numpy as np

# import project specific 3rd party
import pya
from pya import Point, Vector, DPoint, DVector, DEdge, Edges, Region
from pya import DSimplePolygon, SimplePolygon, DPolygon, DBox, Polygon

from pya import DCplxTrans, Trans, DTrans, ICplxTrans

# import project lib
import classLib

reload(classLib)
from classLib.coplanars import CPW, CPW2CPW, CPWParameters, DPathCPW, Bridge1, Intersection
from sonnetSim import SonnetLab, SonnetPort, SimulationBox
from classLib.chipDesign import ChipDesign, GlobalDesignParameters
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
from designElementsGeometry import QubitParams, Qubit, QubitsGrid
from designElementsGeometry import DiskConn8Pars
from designElementsGeometry import ConnectivityMap
from designElementsGeometry import ROResonator, ROResonatorParams
from designElementsGeometry import CqrCouplingParamsType1

import Cqq_couplings

reload(Cqq_couplings)
from Cqq_couplings import CqqCouplingType2, CqqCouplingParamsType2
from Cqq_couplings import CqqCouplingType1, CqqCouplingParamsType1

class Design12QStair(ChipDesign):
    def __init__(
        self, cell_name,
        global_design_params:GlobalDesignParameters=GlobalDesignParameters()
    ):
        super().__init__(cell_name, global_design_params=GlobalDesignParameters())

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

        self.regions: List[Region] = [
            self.region_ph, self.region_bridges1, self.region_bridges2,
            self.region_el, self.dc_bandage_reg,
            self.region_el_protection, self.region_cut
        ]

        # has to call it once more to add new layers
        # defined in this child of `ChipDesign`
        self.lv.add_missing_layers()

        ''' Connectivity map '''
        self.connectivity_map = ConnectivityMap()

        ''' QUBITS GRID '''
        self.qubits_grid: QubitsGrid = QubitsGrid()
        self.NQUBITS = len(self.qubits_grid.pts_grid)  # 12
        self.qubits: List[Qubit] = [None] * self.NQUBITS
        self.squids: List[AsymSquid] = [None] * self.NQUBITS
        ''' QUBIT COUPLINGS '''
        self.q_couplings: np.array = np.empty(
            (self.NQUBITS, self.NQUBITS),
            dtype=object
        )
        self.circle_hull_d = 5e3
        self.qq_coupling_connectors_map: np.ndarray = ConnectivityMap().qq_coupling_connectors_map

        ''' READOUT RESONATORS '''
        self.resonators_params = ROResonatorParams(qubits_grid=self.qubits_grid)
        self.resonators: List[ROResonator] = [None] * self.NQUBITS
        # TODO: hide into designElementsGeomtry
        self.q_res_connector_idxs: np.ndarray = np.array([4, 4, 0, 4, 4, 0, 0, 4, 0, 0, 4, 0])

        self.q_res_connector_roline_map = self.connectivity_map.q_res_connector_roline_map
        self.q_res_coupling_params: List[
            CqrCouplingParamsType1] = self.resonators_params.q_res_coupling_params

        ''' READOUT LINES '''
        self.ro_lines: List[DPathCPW] = [None] * 2
        self.qCenter_roLine_distance = None

        ''' Microwave control lines '''
        # Z = 49.0538 E_eff = 6.25103 (E = 11.45)
        self.z_md1: CPWParameters = CPWParameters(30e3, 15e3)
        self.z_md2: CPWParameters = CPWParameters(4e3, 4e3)
        self.cpw_md_lines: List[DPathCPW] = [None] * self.NQUBITS
        # length of the smoothing part between normal thick and end-thin cpw for md line
        self.md_line_cpw12_smoothhing = 100e3

        # length ofa thin part of the microwave drive line end
        self.md_line_cpw2_len = 300e3

        ''' Flux control lines '''
        self.cpw_fl_lines: List[DPathCPW] = [None] * self.NQUBITS
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

        ''' CHIP PARAMETERS '''
        self.chip = CHIP
        self.chip_box: pya.DBox = self.chip.box
        # Z = 49.5656 E_eff = 6.30782 (E = 11.45)
        self.z_md_fl: CPWParameters = CPWParameters(10e3, 5.7e3)
        self.ro_line_Z = CPWParameters(width=18e3, gap=10e3)
        self.contact_pads: List[ContactPad] = self.chip.get_contact_pads(
            chip_Z_list=[
                self.ro_line_Z, self.z_md1, self.z_fl1, self.z_md1, self.z_fl1,  # left side
                self.z_fl1, self.z_fl1, self.ro_line_Z, self.z_fl1, self.z_fl1,  # bottom
                self.ro_line_Z, self.z_fl1, self.z_fl1, self.z_fl1, self.z_fl1,  # right
                self.ro_line_Z, self.z_fl1, self.z_md1, self.z_fl1, self.z_md1  # top
            ],
            back_metal_gap=200e3,
            back_metal_width=0e3,
            pad_length=700e3,
            transition_len=250e3
        )

        ''' Bandages avoiding points '''
        self.intersection_points: List[DPoint] = []
        self.resonator_avoid_points: List[DPoint] = []
        self.control_lines_avoid_points: List[DPoint] = []

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
        #
        self.draw_readout_resonators()
        self.draw_microwave_drvie_lines()
        self.draw_flux_control_lines()
        self.draw_readout_lines()

        self.resolve_intersections()

        self.draw_test_structures()
        self.draw_express_test_structures_pads()
        self.draw_bandages()
        self.draw_el_protection()

        self.add_chip_marking(text_bl=DPoint(2.6e6, 1.2e6), chip_name="8Q_0.0.0.0", text_scale=200)

        self.draw_litography_alignment_marks()
        self.draw_bridges()
        self.draw_pinning_holes()
        # 4Q_Disp_Xmon v.0.3.0.8 p.12 - ensure that contact pads has no holes
        for contact_pad in self.contact_pads:
            contact_pad.place(self.region_ph)

        self.extend_photo_overetching()
        self.inverse_destination(self.region_ph)
        # convert to gds acceptable polygons (without inner holes)
        self.region_ph.merge()
        self.resolve_holes()
        # convert to litograph readable format. Litograph can't handle
        # polygons with more than 200 vertices.
        # self.split_polygons_in_layers(max_pts=180)

        # for processes after litographies
        self.draw_cutting_marks()

        # requested by fabrication team
        self.draw_additional_boxes()

    def split_polygons_in_layers(self, max_pts=200):
        # TODO: add to parent class
        self.region_ph = split_polygons(self.region_ph, max_pts)
        self.region_bridges2 = split_polygons(
            self.region_bridges2,
            max_pts
        )
        for poly in self.region_ph:
            if poly.num_points() > max_pts:
                print("exists photo")
        for poly in self.region_ph:
            if poly.num_points() > max_pts:
                print("exists bridge2")

    def resolve_holes(self):
        # TODO: add to parent class
        for reg in self.regions:
            tmp_reg = Region()
            for poly in reg:
                tmp_reg.insert(poly.resolved_holes())
            reg.clear()
            reg |= tmp_reg

        # TODO: the following code is not working (region_bridges's polygons remain the same)
        # for poly in chain(self.region_bridges2):
        #     poly.resolve_holes()

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

    def draw_qubits_array(self, new_disk_r=DiskConn8Pars().disk_r, idx_list=[]):

        if len(idx_list) > 0:
            for qubit_idx in idx_list:
                squid_connector_idx = self.connectivity_map.get_squid_connector_idx(qubit_idx)
                pt = self.qubits_grid.get_pt(qubit_idx)
                qubit_pars = QubitParams(
                    squid_params=SQUID_PARS,
                    qubit_cap_params=DiskConn8Pars(disk_r=new_disk_r),
                    squid_connector_idx=squid_connector_idx
                )
                qubit = Qubit(
                    origin=pt,
                    qubit_params=qubit_pars,
                    postpone_drawing=False
                )
                self.qubits[qubit_idx] = qubit
                self.squids[qubit_idx] = qubit.squid

                qubit.place(self.region_ph, region_id="ph")
                qubit.place(self.region_el, region_id="el")
        else:
            for qubit_idx in range(self.NQUBITS):
                squid_connector_idx = self.connectivity_map.get_squid_connector_idx(qubit_idx)
                pt = self.qubits_grid.get_pt(qubit_idx)
                qubit_pars = QubitParams(
                    squid_params=SQUID_PARS,
                    qubit_cap_params=DiskConn8Pars(disk_r=new_disk_r),
                    squid_connector_idx=squid_connector_idx
                )
                qubit = Qubit(
                    origin=pt,
                    qubit_params=qubit_pars,
                    postpone_drawing=False
                )
                self.qubits[qubit_idx] = qubit
                self.squids[qubit_idx] = qubit.squid
                # TODO: not working, qubit.squid.origin is wrong | partially
                #  dealt with defining origin as a connection in `self.init_regions(
                #  `_refresh_named_connections` has to include `self.origin` by default
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

    def draw_qq_couplings(
        self, donut_metal_width=CqqCouplingParamsType1().donut_metal_width, direct_list=[]
    ):
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

            qq_coupling_connectors_idxs = self.qq_coupling_connectors_map[pt1_1d_idx, pt2_1d_idx]

            if all(
                    [
                        (pt1 - pt2).abs() == 1,  # Manhattan qubits-neighbours
                        (qq_coupling_connectors_idxs >= 0).all(),  # if coupled
                        (len(direct_list) == 0 or (
                                (pt1_1d_idx in direct_list) and (pt2_1d_idx in direct_list)))
                    ]
            ):
                q1_connector_idx = qq_coupling_connectors_idxs[0]
                q2_connector_idx = qq_coupling_connectors_idxs[1]
                qq_coupling = CqqCouplingType1(
                    origin=DPoint(0, 0),
                    params=CqqCouplingParamsType1(
                        donut_metal_width=donut_metal_width,
                        disk1=self.qubits[pt1_1d_idx].disk_cap_shunt,
                        disk2=self.qubits[pt2_1d_idx].disk_cap_shunt,
                        disk1_connector_idx=q1_connector_idx,
                        disk2_connector_idx=q2_connector_idx
                    ),
                    region_id="ph"
                )
                qq_coupling.place(self.region_ph, region_id="ph")
                self.q_couplings[pt1_1d_idx, pt2_1d_idx] = qq_coupling

    def draw_readout_resonators(self, q_idx_direct: int=None):
        """

        Parameters
        ----------
        q_idx_direct : int
            Draw only 1 parcticular resonator
        Returns
        -------

        """
        q_idxs = list(range(12))
        resonator_kw_args_list = list(
            map(
                self.resonators_params.get_resonator_params_by_qubit_idx, q_idxs
            )
        )

        if q_idx_direct is not None:
            q_idx = q_idx_direct
            q_res_connector_idx = self.q_res_connector_idxs[q_idx]

            self.draw_readout_resonator(q_idx, q_res_connector_idx, resonator_kw_args_list)
        else:
            for q_idx, _, q_res_connector_idx, _ in self.q_res_connector_roline_map:
                self.draw_readout_resonator(q_idx, q_res_connector_idx, resonator_kw_args_list)

    def draw_readout_resonator(self, q_idx, q_res_connector_idx, resonator_kw_args_list):
        qubit = self.qubits[q_idx]
        resonator_kw_args = resonator_kw_args_list[q_idx]
        self.q_res_coupling_params[q_idx].disk1_connector_idx = q_res_connector_idx
        self.q_res_coupling_params[q_idx].disk1 = qubit.disk_cap_shunt

        # resonator construction and drawing
        res = ROResonator(
            **resonator_kw_args,
            coupling_pars=self.q_res_coupling_params[q_idx]  # contains qubit coordinates
        )
        res.place(self.region_ph)
        self.resonators[q_idx] = res

    def draw_readout_lines(self):
        # readout line is extended around qubit square in order to
        # fit readout resonators `L_couplings` and left a bit more space
        # for consistent and easy simulation of notch port resonator
        self.qCenter_roLine_distance = abs((self.qubits[10].origin - self.resonators[0].start).x) \
                                       + \
                                       ROResonatorParams.to_line_list[10]
        ro_line_extension = self.qCenter_roLine_distance / 2
        turn_radii = self.ro_line_Z.b * 3
        # readout line 1
        p0_start = self.contact_pads[0].end
        p1_start = p0_start + DPoint(turn_radii, 0)
        p0_end = self.contact_pads[7].end
        # RO line extends by 2 curve radii of the resonator along coupling

        pts_middle = []
        for q_idx in [10, 7, 3, 4, 0, 1]:
            # TODO: fix fucking problem with `origin` not tracking changes after construction had
            #  to use `resonator.start` here for god's sake.
            res_i = self.resonators[q_idx]
            res_i_to_line = self.resonators_params.to_line_list[q_idx]
            res_i_Lcoupling = self.resonators_params.L_coupling_list[q_idx]
            res_i_r = self.resonators_params.res_r_list[q_idx]
            res_i_rotation = DCplxTrans(
                1, self.resonators_params.resonator_rotation_angles[q_idx], False, 0, 0
            )
            p_res_start = res_i.start + \
                          res_i_rotation * DVector(res_i_Lcoupling + 3 * res_i_r, res_i_to_line)
            p_res_end = res_i.start + \
                        res_i_rotation * DVector(-3 * res_i_r, res_i_to_line)
            pts_middle += [p_res_start, p_res_end]

        p1_end = p0_end + DVector(0, 3 * turn_radii)
        p2_end = DVector(pts_middle[-1].x, p1_end.y)
        pts = [p0_start, p1_start] + pts_middle + [p2_end, p1_end, p0_end]
        self.ro_lines[0] = DPathCPW(
            points=pts,
            cpw_parameters=[self.ro_line_Z],
            turn_radii=[turn_radii],
            trans_in=None,
            region_id="ph"
        )
        self.ro_lines[0].place(self.region_ph, region_id="ph")

        # readout line 2
        p0_start = self.contact_pads[15].end
        p1_start = p0_start + DPoint(0, -3 * turn_radii)
        p0_end = self.contact_pads[10].end
        # RO line extends by 2 curve radii of the resonator along coupling

        pts_middle = []
        for q_idx in [11, 8, 9, 5, 6, 2]:
            # TODO: fix fucking problem with `origin` not tracking changes after construction had
            #  to use `resonator.start` here for god's sake.
            res_i = self.resonators[q_idx]
            res_i_to_line = self.resonators_params.to_line_list[q_idx]
            res_i_Lcoupling = self.resonators_params.L_coupling_list[q_idx]
            res_i_r = self.resonators_params.res_r_list[q_idx]
            res_i_rotation = DCplxTrans(
                1, self.resonators_params.resonator_rotation_angles[q_idx], False, 0, 0
            )
            p_res_start = res_i.start + \
                          res_i_rotation * DVector(-3 * res_i_r, res_i_to_line)
            p_res_end = res_i.start + \
                        res_i_rotation * DVector(res_i_Lcoupling + 3 * res_i_r, res_i_to_line)
            pts_middle += [p_res_start, p_res_end]

        p1_end = p0_end + DVector(-turn_radii, 0)

        p2_start = DPoint(pts_middle[0].x, pts_middle[0].y + 2 * turn_radii)
        pts = [p0_start, p1_start, p2_start] + pts_middle + [p1_end, p0_end]
        self.ro_lines[1] = DPathCPW(
            points=pts,
            cpw_parameters=[self.ro_line_Z],
            turn_radii=[turn_radii],
            trans_in=None,
            region_id="ph"
        )
        self.ro_lines[1].place(self.region_ph, region_id="ph")

    def draw_microwave_drvie_lines(self):
        r_turn = 100e3

        q_origin_md_end_d = 220e3

        q_idx = 10
        qubit = self.qubits[q_idx]
        p_start = self.contact_pads[19].end
        p1 = p_start + DVector(0, -0.25e6)
        p2 = DPoint(3.1e6, 12.1e6)
        qubit_center_md_dv_n = p2 - qubit.origin
        qubit_center_md_dv_n /= qubit_center_md_dv_n.abs()
        p_end = qubit.origin + q_origin_md_end_d * qubit_center_md_dv_n
        cpwrl_md = DPathCPW(
            points=[p_start, p1, p2, p_end],
            cpw_parameters=[self.z_md1],
            turn_radii=[r_turn]
        )
        self.cpw_md_lines[q_idx] = cpwrl_md

        q_idx = 11
        qubit = self.qubits[q_idx]
        p_start = self.contact_pads[17].end
        p1 = p_start + DVector(0, -0.25e6)
        p2 = DPoint(5.3e6, 12e6)
        p3 = Point(5.85e6, 9.7e6)
        qubit_center_md_dv_n = p2 - qubit.origin
        qubit_center_md_dv_n /= qubit_center_md_dv_n.abs()
        p_end = qubit.origin + q_origin_md_end_d * qubit_center_md_dv_n
        cpwrl_md = DPathCPW(
            points=[p_start, p1, p2, p3, p_end],
            cpw_parameters=[self.z_md1],
            turn_radii=[r_turn]
        )
        self.cpw_md_lines[q_idx] = cpwrl_md

        q_idx = 7
        qubit = self.qubits[q_idx]
        p_start = self.contact_pads[1].end
        p1 = p_start + DVector(0.25e6, 0)
        p2 = DPoint(3.426e6, 7.077e6)
        p3 = DPoint(4.910e6, 7.077e6)
        qubit_center_md_dv_n = p3 - qubit.origin
        qubit_center_md_dv_n /= qubit_center_md_dv_n.abs()
        p_end = qubit.origin + q_origin_md_end_d * qubit_center_md_dv_n
        cpwrl_md = DPathCPW(
            points=[p_start, p1, p2, p3, p_end],
            cpw_parameters=[self.z_md1],
            turn_radii=[r_turn]
        )
        self.cpw_md_lines[q_idx] = cpwrl_md

        q_idx = 3
        qubit = self.qubits[q_idx]
        p_start = self.contact_pads[3].end
        p1 = p_start + DVector(0.25e6, 0)
        p2 = DPoint(3.72e6, 5.68e6)
        p3 = DPoint(4.56e6, 5.92e6)
        qubit_center_md_dv_n = p3 - qubit.origin
        qubit_center_md_dv_n /= qubit_center_md_dv_n.abs()
        p_end = qubit.origin + q_origin_md_end_d * qubit_center_md_dv_n
        cpwrl_md = DPathCPW(
            points=[p_start, p1, p2, p3, p_end],
            cpw_parameters=[self.z_md1],
            turn_radii=[r_turn]
        )
        self.cpw_md_lines[q_idx] = cpwrl_md

        for i, cpw_md_line in enumerate(self.cpw_md_lines):
            if cpw_md_line is not None:
                self.modify_md_line_end_and_place(
                    cpw_md_line, mod_length=self.md_line_cpw2_len,
                    smoothing=self.md_line_cpw12_smoothhing
                )

    def modify_md_line_end_and_place(self, md_line: DPathCPW, mod_length=100e3, smoothing=20e3):
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
             list(last_lines.values())[:-1]]
        )
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

        q_idx = 10
        qubit = self.qubits[q_idx]
        p_start = self.contact_pads[18].end
        p1 = p_start + DVector(0, -0.25e6)
        p_end = qubit.origin + DVector(
            0,
            qubit.disk_cap_shunt.pars.disk_r + \
            qubit.disk_cap_shunt.pars.disk_gap
        ) + DVector(8.0169e3, 0)
        p_tr_start = DPoint(p_end.x, 9.7e6)
        fl_dpath = DPathCPW(
            points=[p_start, p1, p_tr_start, p_end],
            cpw_parameters=[self.z_fl1],
            turn_radii=[r_turn]
        )
        self.cpw_fl_lines[q_idx] = fl_dpath

        q_idx = 11
        qubit = self.qubits[q_idx]
        p_start = self.contact_pads[16].end
        p1 = p_start + DVector(0, -0.25e6)
        p_end = qubit.origin + DVector(
            0,
            qubit.disk_cap_shunt.pars.disk_r + \
            qubit.disk_cap_shunt.pars.disk_gap
        ) + DVector(8.0169e3, 0)
        p2 = DPoint(6e6, 11e6)
        p_tr_start = DPoint(p_end.x, 9e6)
        fl_dpath = DPathCPW(
            points=[p_start, p1, p2, p_tr_start, p_end],
            cpw_parameters=[self.z_fl1],
            turn_radii=[r_turn]
        )
        self.cpw_fl_lines[q_idx] = fl_dpath

        q_idx = 8
        qubit = self.qubits[q_idx]
        p_start = self.contact_pads[14].end
        p1 = p_start + DVector(-0.25e6, 0)
        p_end = qubit.origin + DVector(
            0,
            qubit.disk_cap_shunt.pars.disk_r + qubit.disk_cap_shunt.pars.disk_gap
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
        self.cpw_fl_lines[q_idx] = fl_dpath

        q_idx = 9
        qubit = self.qubits[q_idx]
        p_start = self.contact_pads[13].end
        p1 = p_start + DVector(-0.35e6, 0)
        p_end = qubit.origin + DVector(
            0,
            qubit.disk_cap_shunt.pars.disk_r + qubit.disk_cap_shunt.pars.disk_gap
        ) + DVector(8.0169e3, 0)
        p2 = DPoint(11.6e6, 8.92e6)
        p3 = DPoint(8.842e6, 9.081e6)
        p4 = DPoint(7.902e6, 8.177e6)
        p_tr_start = p_end + DVector(
            0,
            CqqCouplingParamsType1().bendings_disk_center_d
        )
        fl_dpath = DPathCPW(
            points=[p_start, p1, p2, p3, p4, p_tr_start, p_end],
            cpw_parameters=[self.z_fl1],
            turn_radii=[r_turn]
        )
        self.cpw_fl_lines[q_idx] = fl_dpath

        q_idx = 5
        qubit = self.qubits[q_idx]
        p_start = self.contact_pads[12].end
        p1 = p_start + DVector(-3 * r_turn, 0)
        p2 = p1 + DVector(-0.25e6, 0)
        p_end = qubit.origin + DVector(
            0,
            qubit.disk_cap_shunt.pars.disk_r + qubit.disk_cap_shunt.pars.disk_gap
        ) + DVector(8.0169e3, 0)
        p3 = DPoint(11.1e6, 8.31e6)
        p4 = DPoint(9.525e6, 8.414e6)
        p5 = DPoint(8.632e6, 7.466e6)
        p_tr_start = p_end + DVector(
            0,
            CqqCouplingParamsType1().bendings_disk_center_d / 2
        )
        fl_dpath = DPathCPW(
            points=[p_start, p1, p2, p3, p4, p5, p_tr_start, p_end],
            cpw_parameters=[self.z_fl1],
            turn_radii=[r_turn]
        )
        self.cpw_fl_lines[q_idx] = fl_dpath

        q_idx = 6
        qubit = self.qubits[q_idx]
        p_start = self.contact_pads[11].end
        p1 = p_start + DVector(-3 * r_turn, 0)
        p_end = qubit.origin + DVector(
            0,
            qubit.disk_cap_shunt.pars.disk_r + qubit.disk_cap_shunt.pars.disk_gap
        ) + DVector(8.0169e3, 0)
        p2 = DPoint(10.3e6, 7.2e6)
        p_tr_start = p_end + DVector(
            0,
            CqqCouplingParamsType1().bendings_disk_center_d / 2
        )
        fl_dpath = DPathCPW(
            points=[p_start, p1, p2, p_tr_start, p_end],
            cpw_parameters=[self.z_fl1],
            turn_radii=[r_turn]
        )
        self.cpw_fl_lines[q_idx] = fl_dpath

        q_idx = 7
        p_start = self.contact_pads[2].end
        p1 = p_start + DVector(0.25e6, 0)
        p2 = DPoint(3.426e6, 6.905e6)
        qubit = self.qubits[q_idx]
        squid_connector_idx = self.connectivity_map.get_squid_connector_idx(qubit_idx=q_idx)
        fl_line_connector_angle = qubit.qubit_params.qubit_cap_params.connector_angles[
            squid_connector_idx
        ]
        _p_tr_trans = DCplxTrans(1, fl_line_connector_angle, False, 0, 0)
        p_tr_start = qubit.origin + _p_tr_trans * DVector(
            CqqCouplingParamsType1().bendings_disk_center_d, -8.0169e3
        )
        p3 = DPoint(4.9e6, 6.905e6)
        p_end = qubit.origin + _p_tr_trans * DVector(
            qubit.disk_cap_shunt.pars.disk_r +
            qubit.disk_cap_shunt.pars.disk_gap,
            -8.0169e3
        )
        fl_dpath = DPathCPW(
            points=[p_start, p1, p2, p3, p_tr_start, p_end],
            cpw_parameters=[self.z_fl1],
            turn_radii=[r_turn]
        )
        self.cpw_fl_lines[q_idx] = fl_dpath

        q_idx = 3
        p_start = self.contact_pads[4].end
        p1 = p_start + DVector(0.25e6, 0)
        p2 = DPoint(3.87e6, 5.21e6)
        p3 = DPoint(4.89e6, 5.79e6)
        qubit = self.qubits[q_idx]
        squid_connector_idx = self.connectivity_map.get_squid_connector_idx(qubit_idx=q_idx)
        fl_line_connector_angle = qubit.qubit_params.qubit_cap_params.connector_angles[
            squid_connector_idx
        ]
        _p_tr_trans = DCplxTrans(1, fl_line_connector_angle, False, 0, 0)
        p_tr_start = qubit.origin + _p_tr_trans * DVector(
            CqqCouplingParamsType1().bendings_disk_center_d, -8.0169e3
        )
        p_end = qubit.origin + _p_tr_trans * DVector(
            qubit.disk_cap_shunt.pars.disk_r +
            qubit.disk_cap_shunt.pars.disk_gap,
            -8.0169e3
        )
        fl_dpath = DPathCPW(
            points=[p_start, p1, p2, p3, p_tr_start, p_end],
            cpw_parameters=[self.z_fl1],
            turn_radii=[r_turn]
        )
        self.cpw_fl_lines[q_idx] = fl_dpath

        q_idx = 4
        qubit = self.qubits[q_idx]
        p_start = self.contact_pads[5].end
        p1 = p_start + DVector(0, r_turn)
        p_end = qubit.origin + (-1) * DVector(
            8.0169e3,
            qubit.disk_cap_shunt.pars.disk_r + \
            qubit.disk_cap_shunt.pars.disk_gap
        )
        p2 = DPoint(4.70e6, 4.43e6)
        p3 = DPoint(5.57e6, 5.36e6)
        p_tr_start = p_end - DVector(
            0,
            CqqCouplingParamsType1().bendings_disk_center_d
        )
        fl_dpath = DPathCPW(
            points=[p_start, p1, p2, p3, p_tr_start, p_end],
            cpw_parameters=[self.z_fl1],
            turn_radii=[r_turn]
        )
        self.cpw_fl_lines[q_idx] = fl_dpath

        q_idx = 0
        qubit = self.qubits[q_idx]
        p_start = self.contact_pads[6].end
        p1 = p_start + DVector(0, r_turn)
        p_end = qubit.origin + (-1) * DVector(
            8.0169e3,
            qubit.disk_cap_shunt.pars.disk_r + \
            qubit.disk_cap_shunt.pars.disk_gap
        )
        p2 = DPoint(6.7e6, 3.5e6)
        p3 = DPoint(p2.x, 5.0e6)
        p_tr_start = p_end - DVector(
            0,
            CqqCouplingParamsType1().bendings_disk_center_d
        )
        fl_dpath = DPathCPW(
            points=[p_start, p1, p2, p3, p_tr_start, p_end],
            cpw_parameters=[self.z_fl1],
            turn_radii=[r_turn]
        )
        self.cpw_fl_lines[q_idx] = fl_dpath

        q_idx = 1
        qubit = self.qubits[q_idx]
        p_start = self.contact_pads[8].end
        p1 = p_start + DVector(0, r_turn)
        p_end = qubit.origin + (-1) * DVector(
            8.0169e3,
            qubit.disk_cap_shunt.pars.disk_r + \
            qubit.disk_cap_shunt.pars.disk_gap
        )
        p_tr_start = p_end - DVector(
            0,
            CqqCouplingParamsType1().bendings_disk_center_d
        )
        fl_dpath = DPathCPW(
            points=[p_start, p1, p_tr_start, p_end],
            cpw_parameters=[self.z_fl1],
            turn_radii=[r_turn]
        )
        self.cpw_fl_lines[q_idx] = fl_dpath

        q_idx = 2
        qubit = self.qubits[q_idx]
        p_start = self.contact_pads[9].end
        p1 = p_start + DVector(0, r_turn)
        p_end = qubit.origin + (-1) * DVector(
            8.0169e3,
            qubit.disk_cap_shunt.pars.disk_r + \
            qubit.disk_cap_shunt.pars.disk_gap
        )
        p_tr_start = p_end - DVector(
            0,
            CqqCouplingParamsType1().bendings_disk_center_d
        )
        fl_dpath = DPathCPW(
            points=[p_start, p1, p_tr_start, p_end],
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

    def resolve_intersections(self):
        for ctr_lines, ro_line in itertools.product(
                zip(self.cpw_md_lines, self.cpw_fl_lines),
                self.ro_lines
        ):
            for dpathcpw_line in ctr_lines:
                # some qubit do not have control line of the certain type
                if dpathcpw_line is None:
                    continue

                # ro line is placed later and is to be preserved continuous
                result = Intersection.get_intersected_cpws(dpathcpw_line, ro_line)
                if result is not None:
                    cpw1, cpw2 = result
                    intersection_pt = Intersection.resolve_cpw_cpw_intersection(
                        cpw1=cpw1, cpw2=cpw2, cpw_reg=self.region_ph,
                        bridge_reg1=self.region_bridges1, bridge_reg2=self.region_bridges2
                    )
                    self.intersection_points.append(intersection_pt)

    def draw_test_structures(self):
        # TODO: fix test structure SQUIDs
        ''' Triplet of test structures for separate SQUID's JJ and bridges '''
        struct_centers = [
            DPoint(1.8e6, 6.0e6),
            DPoint(11e6, 10.7e6),
            DPoint(4.5e6, 3e6)
        ]
        test_struct_gnd_gap = 20e3  # TODO: move to designElementsGeometry
        test_struct_pads_gap = 40e3
        for struct_center in struct_centers:
            ## JJ test structures ##
            dx = SQUID_PARS.SQB_dx / 2 - SQUID_PARS.SQLBT_dx / 2

            # test structure for left JJ's critical currents
            test_struct1 = TestStructurePadsSquare(
                struct_center,
                # gnd gap in test structure is now equal to
                # the same of first xmon cross, where polygon is placed
                gnd_gap=20e3,
                pads_gap=(
                        self.qubits[0].disk_cap_shunt.pars.disk_gap -
                        (
                                self.qubits[0].squid_pimp.length() -
                                self.qubits[0].disk_cap_shunt.pars.disk_r
                        )
                )
            )
            self.test_squids_pads.append(test_struct1)
            test_struct1.place(self.region_ph)

            text_reg = pya.TextGenerator.default_generator().text(
                "56 nA", 0.001, 25, False, 0, 0
            )
            text_bl = test_struct1.empty_rectangle.p1 - DVector(0, 20e3)
            text_reg.transform(
                ICplxTrans(1.0, 0, False, text_bl.x, text_bl.y)
            )
            self.region_ph -= text_reg

            pars_local = copy.deepcopy(SQUID_PARS)
            pars_local.SQRBT_dx = 0
            pars_local.SQRBJJ_dy = 0
            pars_local.bot_wire_x = [-dx]

            squid_center = test_struct1.center
            test_jj = AsymSquid(
                squid_center + DVector(0, -8.0001234e3) + DVector(0, 5.126e3),
                pars_local
            )
            self.test_squids.append(test_jj)
            test_jj.place(self.region_el)

            # test structure with low critical current (#2)
            test_struct2 = TestStructurePadsSquare(
                struct_center + DPoint(0.3e6, 0),
                gnd_gap=20e3,
                pads_gap=(
                        self.qubits[0].disk_cap_shunt.pars.disk_gap -
                        (
                                self.qubits[0].squid_pimp.length() -
                                self.qubits[0].disk_cap_shunt.pars.disk_r
                        )
                )
            )
            self.test_squids_pads.append(test_struct2)
            test_struct2.place(self.region_ph)

            text_reg = pya.TextGenerator.default_generator().text(
                "11 nA", 0.001, 25, False, 0, 0
            )
            text_bl = test_struct2.empty_rectangle.p1 - DVector(0, 20e3)
            text_reg.transform(
                ICplxTrans(1.0, 0, False, text_bl.x, text_bl.y)
            )
            self.region_ph -= text_reg

            pars_local = copy.deepcopy(SQUID_PARS)
            pars_local.SQLBT_dx = 0
            pars_local.SQLBJJ_dy = 0
            pars_local.bot_wire_x = [dx]

            squid_center = test_struct2.center
            test_jj = AsymSquid(
                squid_center + DVector(0, -8.0001234e3) + DVector(0, 5.126e3),
                pars_local
            )
            self.test_squids.append(test_jj)
            test_jj.place(self.region_el)

            # test structure for bridge DC contact (#3)
            test_struct3 = TestStructurePadsSquare(
                struct_center + DPoint(0.6e6, 0)
            )
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
            DPoint(4.5e6, 2e6)
        ]
        for struct_center in test_dc_el2_centers:
            test_struct1 = TestStructurePadsSquare(struct_center)
            test_struct1.place(self.region_ph)
            text_reg = pya.TextGenerator.default_generator().text(
                "Bandage", 0.001, 40, False, 0, 0
            )
            text_bl = test_struct1.empty_rectangle.origin + DPoint(
                test_struct1.gnd_gap, -4 * test_struct1.gnd_gap
            )
            text_reg.transform(
                ICplxTrans(1.0, 0, False, text_bl.x, text_bl.y)
            )
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

    def draw_el_protection(self):
        protection_a = 300e3
        for squid in (self.squids + self.test_squids):
            self.region_el_protection.insert(
                pya.Box().from_dbox(
                    pya.DBox(
                        squid.center - 0.5 * DVector(
                            protection_a,
                            protection_a
                            ),
                        squid.center + 0.5 * DVector(
                            protection_a,
                            protection_a
                            )
                    )
                )
            )

    def _draw_squid_bandage(self, squid: AsymSquid = None, shift2sq_center=0):
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

        # for readout resonators
        for resonator in self.resonators:
            for name, res_primitive in resonator.primitives.items():
                if "coil" in name:
                    subprimitives = res_primitive.primitives
                    for primitive_name, primitive in subprimitives.items():
                        # place bridges only at arcs of coils
                        # but not on linear segments
                        # Note: "coil0" contains resonator-waveguide coupling as "CPW", hence it
                        # is excluded
                        if "arc" in primitive_name:
                            Bridge1.bridgify_CPW(
                                cpw=primitive,
                                bridges_step=bridges_step, gnd2gnd_dy=70e3,
                                dest=self.region_bridges1, dest2=self.region_bridges2
                            )
                elif any([(skip_name in name) for skip_name in ["fork", "arc"]]):
                    # skip primitives with names above
                    continue
                elif any(
                        [
                            (primitive_type in res_primitive.__class__.__name__)
                            for primitive_type in ["DPathCPW", "CPW", "CPWArc"]
                        ]
                ):
                    # bridgify the rest
                    Bridge1.bridgify_CPW(
                        cpw=res_primitive,
                        bridges_step=bridges_step, gnd2gnd_dy=70e3,
                        dest=self.region_bridges1, dest2=self.region_bridges2
                    )

        ''' contact wires '''
        # TODO: fix this
        self.control_lines_avoid_points += [squid.center for squid in self.squids]
        self.control_lines_avoid_points += self.intersection_points
        for ctr_line in itertools.chain(self.cpw_md_lines, self.cpw_fl_lines):
            if ctr_line is None:
                continue
            else:
                Bridge1.bridgify_CPW(
                    cpw=ctr_line,
                    bridges_step=bridges_step, gnd2gnd_dy=100e3,
                    dest=self.region_bridges1, dest2=self.region_bridges2,
                    avoid_points=self.control_lines_avoid_points,
                    avoid_distances=[200e3]
                )

        # close bridges for flux contact wires
        for q_idx, cpw_fl in enumerate(self.cpw_fl_lines):
            if cpw_fl is None:
                continue
            dy_list = [30e3, 100e3]
            for dy in dy_list:
                squid_connector_idx = self.connectivity_map.get_squid_connector_idx(qubit_idx=q_idx)
                if squid_connector_idx == 2:
                    pass
                else:
                    dy = -dy
                bridge_center1 = cpw_fl.end + DVector(0, dy)
                br = Bridge1(
                    center=bridge_center1, gnd2gnd_dy=70e3,
                    trans_in=Trans.R90
                )
                br.place(
                    dest=self.region_bridges1,
                    region_id="bridges_1"
                )
                br.place(
                    dest=self.region_bridges2,
                    region_id="bridges_2"
                )

        ''' readout lines  '''
        # collecting avoid points
        # TODO: refactor such that resonator bridgifying process utilizes
        #  `self.resonator_avoid_points`
        self.resonator_avoid_points = []
        resonator_avoid_distances = []
        for res in self.resonators:
            # avoid placing bridges close to thin resonator coupling
            av_pt = res.primitives["coil0"].primitives["cop1"].center()
            self.resonator_avoid_points.append(av_pt)
            avoid_distance = res.L_coupling / 2 + res.r + res.Z0.b / 2
            resonator_avoid_distances.append(avoid_distance)

        ro_lines_avoid_points = self.resonator_avoid_points + self.intersection_points
        # 150 um from intersectoins
        resonator_avoid_distances += [150e3] * len(self.intersection_points)
        for ro_line in self.ro_lines:
            Bridge1.bridgify_CPW(
                cpw=ro_line,
                bridges_step=bridges_step, gnd2gnd_dy=100e3,
                dest=self.region_bridges1, dest2=self.region_bridges2,
                avoid_points=ro_lines_avoid_points,
                avoid_distances=resonator_avoid_distances
            )

        ''' Cqq couplings '''
        for cqq_coupling in self.q_couplings.flat:
            if cqq_coupling is not None:
                Bridge1.bridgify_CPW(
                    cpw=cqq_coupling.primitives["cpw_central"],
                    bridges_step=bridges_step, gnd2gnd_dy=100e3,
                    dest=self.region_bridges1, dest2=self.region_bridges2,
                )

    def draw_pinning_holes(self):
        # points that select polygons of interest if they were clicked at)
        selection_pts = [
            Point(0.1e6, 0.1e6),
            Point(self.chip.dx - 0.1e6, self.chip.dy - 0.1e6),
            DPoint(9.7e6, 2.5e6),
            DPoint(0.3e6, 13.7e6)
        ]

        # append grounded polygons inside qubits grid
        for q_idx in [0, 1, 3, 4, 7]:
            selection_pts.append(
                self.qubits[q_idx].origin +
                DVector(self.qubits_grid.dx / 2, self.qubits_grid.dy / 2)
            )

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
        filled_reg = fill_holes(
            reg_to_fill, d=40e3, width=15e3,
            height=15e3
        )

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
        super()._transfer_regs2cell()
        self.cell.shapes(self.dc_bandage_layer).insert(self.dc_bandage_reg)
        self.cell.shapes(self.layer_bridges1).insert(self.region_bridges1)
        self.cell.shapes(self.layer_bridges2).insert(self.region_bridges2)
        self.cell.shapes(self.layer_el_protection).insert(
            self.region_el_protection
        )
        self.lv.zoom_fit()

    def draw_additional_boxes(self):
        abox_top_ph = pya.Box(
            Point(self.chip.dx / 2, self.chip.dy / 2) + Point(
                -self.chip.dx * 0.3, self.chip.dx * 0.52
            ),
            Point(self.chip.dx / 2, self.chip.dy / 2) + Point(
                self.chip.dx * 0.3, self.chip.dx * 0.62
            )
        )
        abox_bot_ph = pya.Box(
            Point(self.chip.dx / 2, self.chip.dy / 2) - Point(
                -self.chip.dx * 0.3, self.chip.dx * 0.52
            ),
            Point(self.chip.dx / 2, self.chip.dy / 2) - Point(
                self.chip.dx * 0.3, self.chip.dx * 0.62
            )
        )
        self.region_ph.insert(abox_top_ph)
        self.region_ph.insert(abox_bot_ph)

        abox_top_el = pya.Box(
            Point(self.chip.dx / 2, self.chip.dy / 2) + Point(
                -self.chip.dx * 0.35, self.chip.dx * 0.54
            ),
            Point(self.chip.dx / 2, self.chip.dy / 2) + Point(
                self.chip.dx * 0.35, self.chip.dx * 0.6
            )
        )
        abox_bot_el = pya.Box(
            Point(self.chip.dx / 2, self.chip.dy / 2) - Point(
                -self.chip.dx * 0.35, self.chip.dx * 0.54
            ),
            Point(self.chip.dx / 2, self.chip.dy / 2) - Point(
                self.chip.dx * 0.35, self.chip.dx * 0.6
            )
        )
        self.region_bridges1.insert(abox_top_el)
        self.region_bridges1.insert(abox_bot_el)

        ext_chip_box = self.chip_box.dup()
        ext_chip_box.left -= 2e6
        ext_chip_box.bottom -= 2e6
        ext_chip_box.top += 2e6
        ext_chip_box.right += 2e6
        ext_chip_box = Region(ext_chip_box)
        ext_chip_box -= Region(self.chip_box)
        self.region_bridges2 += ext_chip_box

    def draw_for_res_Q_sim(self, q_idx, border_down=360e3, border_up=200e3):
        self.draw_chip()
        self.draw_qubits_array()
        self.draw_qq_couplings()

        self.draw_readout_resonators(q_idx_direct=q_idx)

        rotated_angle_deg = ROResonatorParams.resonator_rotation_angles[q_idx]
        rotation_center = self.qubits[q_idx].origin.dup()

        self.region_ph.transform(ICplxTrans(1, 0, False, DVector(-rotation_center)))
        self.region_ph.transform(ICplxTrans(1, -rotated_angle_deg, False, 0, 0))
        self.region_ph.transform(ICplxTrans(1, 0, False, DVector(rotation_center)))

        to_line = ROResonatorParams.to_line_list[q_idx]

        p1 = self.qubits[q_idx].origin
        p2 = self.resonators[q_idx].start
        p2 = ICplxTrans(1, -rotated_angle_deg, False, 0, 0) * (
                p2 - rotation_center) + rotation_center
        p2 += DVector(self.resonators[q_idx].L_coupling, to_line)
        p2 = DPoint(p2)
        print(p1)
        print(p2)
        cent = (p1 + p2) / 2
        bwidth = abs(p1.x - p2.x)
        bheight = abs(p1.y - p2.y)

        readline = CPW(
            start=cent + DVector(bwidth / 2 + border_up, bheight / 2),
            end=cent + DVector(-bwidth / 2 - border_down, bheight / 2),
            cpw_params=CPWParameters(width=20e3, gap=10e3)
        )
        readline.place(self.region_ph)

        dv = DVector(bwidth, bheight)
        crop_box = DBox(
            cent + dv / 2 + DVector(border_up, border_up),
            cent - dv / 2 - DVector(border_down, border_down)
        )
        self.crop(crop_box)

        self.sonnet_ports.append(readline.start)
        self.sonnet_ports.append(readline.end)

        self.transform_region(
            self.region_ph, DTrans(-(cent - dv / 2 - DVector(border_down, border_down))),
            trans_ports=True
        )

        return crop_box


def simulate_res_f_and_Q(q_idx, resolution=(2e3, 2e3), type='freq'):
    ### DRAWING SECTION START ###
    design = Design12QStair("testScript")
    crop_box = design.draw_for_res_Q_sim(q_idx)
    design.show()
    ### DRAWING SECTION END ###

    cur_freqs = [7.395, 7.465, 7.65, 7.485]

    if type == 'freq':
        simulate_S_pars(
            design, crop_box,
            f'res_{q_idx}_{design.resonators_params.L1_list[q_idx] / 1e3:.01f}_S_pars.csv', 7.0, 8.0
        )
    elif type == 'Q':
        simulate_S_pars(
            design, crop_box,
            f'res_{q_idx}_Q_S_pars.csv',
            ROResonatorParams.target_freqs[q_idx] - 0.01,
            ROResonatorParams.target_freqs[q_idx] + 0.01,
            resolution=resolution
        )


def simulate_S_pars(design, crop_box, filename, min_freq=6.0, max_freq=7.0, resolution=(2e3, 2e3)):
    ### SIMULATION SECTION START ###
    ml_terminal = SonnetLab()
    from sonnetSim.cMD import CMD

    resolution_dx = resolution[0]
    resolution_dy = resolution[1]

    ml_terminal._send(CMD.SAY_HELLO)
    ml_terminal.clear()
    simBox = SimulationBox(
        crop_box.width(), crop_box.height(),
        crop_box.width() / resolution_dx,
        crop_box.height() / resolution_dy
    )

    ml_terminal.set_boxProps(simBox)
    from sonnetSim.pORT_TYPES import PORT_TYPES

    ports = [
        SonnetPort(prt, PORT_TYPES.AUTOGROUNDED) for prt in design.sonnet_ports
    ]
    ml_terminal.set_ports(ports)
    ml_terminal.send_polygons(design.cell, design.layer_ph)
    ml_terminal.set_ABS_sweep(min_freq, max_freq)
    # print(f"simulating...{resonator_idx}")
    result_path = ml_terminal.start_simulation(wait=True)
    ml_terminal.release()

    # creating directory with simulation results
    output_projpath = os.path.join(
        PROJECT_DIR,
        filename
    )

    shutil.copy(
        result_path.decode("ascii"),
        output_projpath
    )

    design.layout.write(output_projpath[:-4] + '.gds')

    ### RESULT SAVING SECTION END ###


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


def simulate_Cqr(q_idxs: List[int], resolution=(4e3, 4e3)):
    # TODO: 1. make 2d geometry parameters mesh, for simultaneous finding of C_qr and C_q
    #  2. make 3d geometry optimization inside kLayout for simultaneous finding of C_qr,
    #  C_q and C_qq
    # dl_list = np.linspace(-10e3, 10e3, 21)
    donut_disk_d_list = np.array([20e3], dtype=float)
    donut_metal_width_list = np.linspace(10e3, 80e3, 8)
    # dl_list = [0e3]

    for simulation_i, (donut_disk_d, donut_metal_width, q_idx) in list(
            enumerate(
                list(
                    itertools.product(
                        donut_disk_d_list, donut_metal_width_list, q_idxs
                    )
                )
            )
    ):
        ### DRAWING SECTION START ###
        design = Design12QStair("testScript")
        print("simulation #", simulation_i)
        # exclude coils from simulation (sometimes port is placed onto coil (TODO: fix)
        design.q_res_coupling_params[q_idx].donut_disk_d = donut_disk_d
        design.q_res_coupling_params[q_idx].donut_metal_width = donut_metal_width
        design.resonators_params.N_coils_list = [1] * design.NQUBITS
        design.draw_chip()
        design.draw_qubits_array()
        design.draw_qq_couplings()
        design.draw_readout_resonators()
        # design.draw_readout_lines()

        resonator = design.resonators[q_idx]
        res_reg = resonator.metal_region
        qubit = design.qubits[q_idx]
        q_reg = qubit.disk_cap_shunt.metal_region

        # TODO: make simulation such that all polygons (except those with ports are connected to
        #  ground). Now it is tolerable to have some of them (with large capacity to ground) to
        #  stay with floating potential (it will be close to ground plane potential due to their
        #  respectively large capacitance to ground i.e. low impedance to ground).

        box_region = q_reg.sized(1 * q_reg.bbox().width(), 1 * q_reg.bbox().height())

        C1, C12, C2 = simulate_cij(
            design=design, layer=design.layer_ph,
            subregs=[q_reg, res_reg], env_reg=box_region,
            resolution=resolution, print_values=True
        )
        print()
        '''SAVING REUSLTS SECTION START'''
        output_filepath = os.path.join(
            PROJECT_DIR,
            f"Xmon_Cqr_results.csv"
        )
        design.save_as_gds2(
            filename=os.path.join(
                PROJECT_DIR,
                f"ddd_{int(donut_disk_d / 1e3)}_dmw_{int(donut_metal_width / 1e3)}.gds"
            )
        )

        save_sim_results(
            output_filepath=output_filepath,
            design=design,
            additional_pars={
                "q_idx"                : q_idx,
                "resolution"           : resolution,
                "donut_disk_d, um"     : design.q_res_coupling_params[q_idx].donut_disk_d,
                "donut_metal_width, um": design.q_res_coupling_params[q_idx].donut_metal_width,
                "C1, fF"               : C1,
                "C2, fF"               : C2,
                "C12, fF"              : C12
            }
        )


def simulate_md_Cg(q_idx: int, resolution=(5e3, 5e3)):
    # dl_list = np.linspace(-20e3, 20e3, 3)
    dl_list = [0]
    for dl in dl_list:
        design = Design12QStair("testScript")

        design.draw_chip()
        design.draw_qubits_array()
        design.draw_qq_couplings()
        #
        design.draw_readout_resonators(q_idx_direct=q_idx)
        design.draw_microwave_drvie_lines()

        ''' CROP BOX CALCULATION SECTION START '''
        qubit = design.qubits[q_idx]
        q_reg = qubit.disk_cap_shunt.metal_region
        md_line = design.cpw_md_lines[q_idx]
        md_reg = md_line.metal_region

        # TODO: make simulation such that all polygons (except those with ports are connected to
        #  ground). Now it is tolerable to have some of them (with large capacity to ground) to
        #  stay with floating potential (it will be close to ground plane potential due to their
        #  respectively large capacitance to ground i.e. low impedance to ground).

        md_dv_n = -list(md_line.primitives.values())[-1].dr  # goes out from qubit towards md line
        md_dv_n /= md_dv_n.abs()
        md_importance_length = design.md_line_cpw12_smoothhing + 3/2*design.md_line_cpw2_len

        box_region = q_reg.sized(1 * q_reg.bbox().width(), 1 * q_reg.bbox().height()) + \
                     Region(  # add a box-point to region
                         pya.Box(
                             qubit.origin + 2*md_importance_length*md_dv_n,
                             qubit.origin + 2*md_importance_length*md_dv_n + DVector(1e3, 1e3)
                         )
                     )
        # and take their smallest enclosing Region
        box_region = Region(box_region.bbox())

        C1, C12, C2 = simulate_cij(
            design=design, layer=design.layer_ph,
            subregs=[q_reg, md_reg], env_reg=box_region,
            resolution=resolution, print_values=True
        )
        print()
        '''SAVING REUSLTS SECTION START'''
        design.save_as_gds2(
            filename=os.path.join(
                PROJECT_DIR,
                f"Cg_md_{q_idx}.gds"
            )
        )

        output_filepath = os.path.join(
            PROJECT_DIR,
            f"Cmd_q_results.csv"
        )
        save_sim_results(
            output_filepath=output_filepath,
            design=design,
            additional_pars={
                "q_idx"                       : q_idx,
                "resolution"                  : resolution,
                "md_line_cpw2_len, um"        : design.md_line_cpw2_len,
                "md_line_cpw12_smoothhing, um": design.md_line_cpw12_smoothhing,
                "C1, fF"                      : C1,
                "C2, fF"                      : C2,
                "C12, fF"                     : C12
            }
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
    #         __version__ + "_overetching_0um.gds"
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
    # simulate_Cqr(q_idxs=[0], resolution=(2e3, 2e3))

    ''' Simulation of C_{q1,q2} in fF '''
    # simulate_Cqq(q1_idx=5, q2_idx=6, resolution=(2e3, 2e3))

    ''' MD line C_qd for md1,..., md6 '''
    # for q_idx in [10, 11]:
    #     simulate_md_Cg(q_idx=q_idx, resolution=(4e3, 4e3))

    ''' Resonators Q and f sim'''
    # simulate_resonators_f_and_Q(resolution=(2e3, 2e3))

    ''' Resonators Q and f when placed together'''
    # simulate_resonators_f_and_Q_together()
    # for i in [0, 1, 2, 3]:
    #    simulate_res_f_and_Q(i)
