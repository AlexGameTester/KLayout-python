# import built-ins
from typing import List
from importlib import reload
from dataclasses import dataclass

# import good 3rd party
import numpy as np

# import project specific 3rd party
import pya
from pya import Point, Vector, DPoint, DVector, DEdge, \
    DSimplePolygon, \
    SimplePolygon, DPolygon, DBox, Polygon, Region

from pya import DCplxTrans, DTrans, Trans

# import self-made API
import classLib

reload(classLib)

from classLib.baseClasses import ComplexBase
from classLib.coplanars import CPW, CPWParameters
from classLib.josJ import AsymSquidParams, AsymSquid

# project specific files
from globalDefinitions import CHIP


@dataclass  # dataclass is used for simplicity of class declaration and readability only
class QubitsGrid:
    # in fractions of chip dimensions
    origin: float = DVector(CHIP.dx / 2, CHIP.dy / 2)
    # step of 2D grid in `x` and `y` directions correspondingly
    dx: float = 1e6
    dy: float = 1e6
    pts_grid: np.ndarray = np.array(
        [
            # grid iterates from left to right (direct),
            # from bottom to top (style of recording is inverted),
            # starting from bl corner.
                    (0, 0), (1, 0), (2, 0),
            (-1, 1), (0, 1), (1, 1), (2, 1),
            (-1, 2), (0, 2), (1, 2),
            (-1, 3), (0, 3)
        ],
        dtype=float
    )

    def __post_init__(self):
        self.__centralize_grid()

    def __centralize_grid(self):
        # grid is centralized such that bbox center of the grid has
        # origin (0,0)
        grid_center = np.array([0.5, 1.5], dtype=float)
        for i, pt_pos in enumerate(self.pts_grid):
            self.pts_grid[i] -= grid_center

    def get_pt(self, idx) -> DPoint:
        pt_pos = self.pts_grid[idx]
        pt_x = pt_pos[0] * self.dx
        pt_y = pt_pos[1] * self.dy
        origin = self.origin
        return origin + DVector(pt_x, pt_y)


@dataclass  # dataclass is used for simplicity of class declaration and readability only
class DiskConn8Pars:
    connector_angles = np.linspace(0, 360, 8, endpoint=False)
    disk_r: float = 120e3
    disk_gap = 20e3
    # extension beyond `disk_r` of `CPW` that imitates connection
    pimp_l = 0e3
    conn_width = 0e3
    # side ground-gap of `CPW` that imitates connection
    conn_side_gap = 0e3
    # front gap of the `CPW` that imitates connection
    conn_front_gap = 20e3


class DiskConn8(ComplexBase):
    """
    Single superconducting Island represents shunting ground capacitor
    for qubit
    Has 8 connections distributed equidistantly (in angle) in counter-clockwise direction.
    Connections are represented by zero-width/zero-gap coplanars.
    As a result, coplanars are aligned into star-like structure.
    See relevant parameters in `DiskConn8Pars` definition.
    """

    def __init__(
        self, origin,
        pars: DiskConn8Pars = DiskConn8Pars(),
        trans_in=None,
        region_id="ph"
        ):
        self.pars: DiskConn8Pars = pars
        self.disk: DiskConn8 = None
        self.conn8_list: List[CPW] = []
        super().__init__(
            origin=origin, trans_in=trans_in, region_id=region_id
        )

    def init_primitives(self):
        # print(self.origin)
        origin = DPoint(0, 0)

        from classLib.shapes import Disk
        self.empty_disk = Disk(
            center=origin, r=self.pars.disk_r + self.pars.disk_gap,
            region_id=self.region_id, inverse=True
        )
        self.primitives["empty_disk"] = self.empty_disk

        # draw star-like connection flanges
        self.angle_connections = self.pars.connector_angles/180*np.pi  # degree
        angles_deg = self.pars.connector_angles
        for i, angle_deg in enumerate(angles_deg):
            trans = DCplxTrans(1, angle_deg, False, 0, 0)
            # finger with side gap
            cpw_l = self.pars.disk_r + self.pars.pimp_l
            cpw = CPW(
                start=origin,
                end=origin + DVector(cpw_l, 0),
                width=self.pars.conn_width,
                gap=self.pars.conn_side_gap,
                open_end_gap=self.pars.conn_front_gap,
                trans_in=trans,
                region_id=self.region_id
            )
            self.conn8_list.append(cpw)  # can be redundant
            # TODO: remove hardcode 3e3 from here
            self.connections.append(cpw.open_end_center + cpw.dr/cpw.dr.abs()*3e3)
            self.primitives["conn" + str(i)] = cpw

        self.disk = Disk(
            center=origin, r=self.pars.disk_r,
            region_id=self.region_id
        )
        self.primitives["circle"] = self.disk

        self.connections.append(origin)

    def _refresh_named_connections(self):
        self.origin = self.connections[-1]


class QubitParams:
    def __init__(
        self,
        squid_params: AsymSquidParams = AsymSquidParams(),
        qubit_cap_params: DiskConn8Pars = DiskConn8Pars(),
        squid_connector_idx=4
    ):
        self.squid_connector_idx = squid_connector_idx
        self.squid_params: AsymSquidParams = squid_params
        self.qubit_cap_params: DiskConn8Pars = qubit_cap_params


class Qubit(ComplexBase):
    _shift_into_substrate = 1.5e3

    def __init__(
        self,
        origin: DPoint = DPoint(0, 0),
        qubit_params: QubitParams = QubitParams(),
        trans_in=None, postpone_drawing=False
    ):
        self.qubit_params = qubit_params
        self.squid: AsymSquid = None
        self.cap_shunt: DiskConn8 = None

        super().__init__(
            origin=origin, trans_in=trans_in, region_id="default",
            postpone_drawing=postpone_drawing, region_ids=["ph", "el"]
        )

    def init_primitives(self):
        """
        Qubit is drawn in its coordinate system with SQUID at the bottom and oriented vertically.
        SQUID is oriented from bottom to top.

        Based on `squid_connector_idx` supplied via `QubitParams` constructor,
        qubit rotates squid relative to its center to position squid properly.

        Returns
        -------

        """
        # TODO `multilayer complex objects`: not called twice due to
        #  emergence of `self.initialized` that prevents several
        #  `init_primitives_trans()` calls
        #  but every layer has to be initialized separately in my opinion
        #  nowadays this makes no difference, so this remark will
        #  remain, until such need arises.
        # print("qubit primitives called", self.origin.x, self.origin.y)
        origin = DPoint(0, 0)

        # draw disk with 8 open-ended CPW connectors
        self.cap_shunt = DiskConn8(
            origin=origin,
            pars=self.qubit_params.qubit_cap_params,
            region_id=self.region_ids[0]
        )

        # SQUID is connected to the right of the disk by default (see schematics for details)
        # SQUID is further rotated according to `QubitParams.squid_connector_idx`
        qubit_finger_end = self.cap_shunt.connections[0]

        # draw squid
        # squid_pars = self.qubit_params.squid_params
        # vertical shift of every squid local origin coordinates
        # this line puts squid center on the
        # "outer capacitor plate's edge" of the shunting capacitance.

        # next step is to make additional shift such that top of the
        # BC polygon is located at the former squid center position with
        # additional `_shift_into_substrate = 1.5e3` um shift in direction of substrate
        # _shift_to_BC_center = squid_pars.shadow_gap / 2 + \
        #                       squid_pars.SQLBT_dy + squid_pars.SQB_dy + \
        #                       squid_pars.BCW_dy + \
        #                       self._shift_into_substrate
        # squid_center += DVector(0, _shift_to_BC_center)

        self.squid = AsymSquid(
            origin=qubit_finger_end,
            params=self.qubit_params.squid_params,
            region_id=self.region_ids[1],
            trans_in=DTrans.R90
        )
        connector_angles = self.qubit_params.qubit_cap_params.connector_angles
        connector_idx = self.qubit_params.squid_connector_idx
        angle = connector_angles[connector_idx]
        self.squid.make_trans(DCplxTrans(1, angle, False, 0, 0))

        self.primitives["squid"] = self.squid
        self.primitives["cap_shunt"] = self.cap_shunt

        self.connections.append(origin)

    def _refresh_named_connections(self):
        self.origin = self.connections[-1]


class ResonatorParams:
    """
    Static class that contains information on readout resonators geometry parameters.
    Geometry parameters has to be verified by simulation.
    """
    # see parameters details in `Design_fast.py`
    L_coupling_list = [
        1e3 * x for x in [310, 320, 320, 310] * 3
    ]
    L0_list = [986e3] * 12
    L1_list = [
        1e3 * x for x in
        [
            114.5219, 95.1897, 99.0318, 83.7159,
             114.5219, 95.1897, 99.0318, 83.7159,
             114.5219, 95.1897, 99.0318, 83.7159
        ]
    ]
    res_r_list = [60e3] * 12
    tail_turn_radiuses_list = [60e3] * 12  # res_r_list
    N_coils_list = [3, 3, 3, 3] * 12
    L2_list = [60e3] * 12  # res_r_list
    L3_list = []  # get numericals from Design_fast
    L4_list = [60e3] * 12  # res_r_list
    Z_res_list = [CPWParameters(10e3, 6e3)]*12
    to_line_list = [45e3] * 12

    # fork at the end of resonator parameters
    fork_metal_width_list = [15e3]*12
    fork_gnd_gap_list = [10e3]*12
    xmon_fork_gnd_gap_list = [14e3]*12
    # self.cross_width_y + 2 * (self.xmon_fork_gnd_gap + self.fork_metal_width)
    fork_x_span_list = [45e3 + 2*(14e3 + 15e3)] * 12
    # resonator-fork parameters
    # from simulation of g_qr
    fork_y_span_list = [
        x * 1e3 for x in
        [
            33.18, 91.43, 39.36, 95.31,
            44.34, 96.58, 49.92, 99.59,
            33.18, 91.43, 39.36, 95.31
        ]
    ]
    tail_segments_list = [[60000.0, 215000.0, 60000.0]]*12
    res_tail_shape = "LRLRL"

    tail_turn_angles_list = [
        [np.pi / 2, -np.pi / 2],
        [np.pi / 2, -np.pi / 2],
        [np.pi / 2, -np.pi / 2],
        [np.pi / 2, -np.pi / 2],
        [-np.pi / 2, np.pi / 2],
        [-np.pi / 2, np.pi / 2],
        [-np.pi / 2, np.pi / 2],
        [-np.pi / 2, np.pi / 2],
        [-np.pi / 2, np.pi / 2],
        [-np.pi / 2, np.pi / 2],
        [-np.pi / 2, np.pi / 2],
        [-np.pi / 2, np.pi / 2]
    ]

    @staticmethod
    def get_resonator_params_by_qubit_idx(q_idx):
        return {
            "Z0": ResonatorParams.Z_res_list[q_idx],
            "L_coupling": ResonatorParams.L_coupling_list[q_idx],
            "L0": ResonatorParams.L0_list[q_idx],
            "L1": ResonatorParams.L1_list[q_idx],
            "r": ResonatorParams.res_r_list[q_idx],
            "N": ResonatorParams.N_coils_list[q_idx],
            "tail_shape": ResonatorParams.res_tail_shape,
            "tail_turn_radiuses": ResonatorParams.tail_turn_radiuses_list[q_idx],
            "tail_segment_lengths": ResonatorParams.tail_segments_list[q_idx],
            "tail_turn_angles": ResonatorParams.tail_turn_angles_list[q_idx],
            "tail_trans_in": Trans.R270,
            "fork_x_span": ResonatorParams.fork_x_span_list[q_idx],
            "fork_y_span": ResonatorParams.fork_y_span_list[q_idx],
            "fork_metal_width": ResonatorParams.fork_metal_width_list[q_idx],
            "fork_gnd_gap": ResonatorParams.fork_gnd_gap_list[q_idx]
        }
