# $epilog: --log=DEBUG
# import built-ins
from typing import List, Dict
from importlib import reload
from dataclasses import dataclass, field
import logging

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
logging.debug("loaded classLib from axillary file")
reload(classLib)
logging.debug("reloaded classLib from axillary file")
from classLib.baseClasses import ComplexBase
from classLib.coplanars import CPW, CPWParameters
from classLib.josJ import AsymSquidParams, AsymSquid
logging.debug("used from classLib.baseClasses import ... | total 3 lines")
# project specific files
from globalDefinitions import CHIP


@dataclass  # dataclass is used for simplicity of class declaration and readability only
class QubitsGrid:
    # in fractions of chip dimensions
    origin: float = DVector(CHIP.dx / 2, CHIP.dy / 2)
    # step of 2D grid in `x` and `y` directions correspondingly
    dx: float = 1e6
    dy: float = 1e6
    pts_grid: np.ndarray = None  # initialized after construction since it is mutable
    NQUBITS: int = None  # number of qubits in `pts_grid`
    # {q_idx: flux} correspondance map. Flux normalized to [0,1] interval
    q_fl_idle: Dict[int, float] = None
    def __post_init__(self):
        # mutable input argument has to be initialized here
        # otherwise this argument will depend on the history of it's modifications
        # look into in PEP 557 or in the docs:
        # https://docs.python.org/3.7/library/dataclasses.html#mutable-default-values
        self.pts_grid = np.array(
            [
                # grid iterates from left to right (direct),
                # from top to bottom
                # starting from the left-most qubit of the top row.
                (-1, 3), (0, 3),
                (-1, 2), (0, 2), (1, 2),
                (-1, 1), (0, 1), (1, 1), (2, 1),
                (0, 0), (1, 0), (2, 0), (3, 0),
                (1, -1), (2, -1), (3, -1)
            ],
            dtype=float
        )
        self.NQUBITS = len(self.pts_grid)

        # used 2 non-neighbour diagonals that qubits placed along. Choosing I and III diagonals,
        # see schematics
        q_uss_list = [
            1, 4, 8,
            12, 3, 7, 11, 15
        ]
        self.q_fl_idle = {}
        for q_idx in range(self.NQUBITS):
            if q_idx in q_uss_list:
                self.q_fl_idle[q_idx] = 0
            else:
                self.q_fl_idle[q_idx] = 0.5

        # shift qubits grid to match with the chip center
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
    disk_r: float = 150e3
    disk_gap = 30e3
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
        self.angle_connections = self.pars.connector_angles / 180 * np.pi  # degree
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
            # center of the disk-ground gap in direction that corresponds to `anlge_deg`
            self.connections.append(
                cpw.dr / cpw.dr.abs() * (self.pars.disk_r +
                                         self.pars.disk_gap / 2)
            )
            self.primitives["conn" + str(i)] = cpw

        self.disk = Disk(
            center=origin, r=self.pars.disk_r,
            region_id=self.region_id
        )
        self.primitives["circle"] = self.disk

        self.connections.append(origin)

    def get_connector_dv(self, connector_idx, normalized=False):
        dr = self.connections[connector_idx] - self.origin
        if normalized:
            dr /= dr.abs()

        return dr

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
        self.disk_cap_shunt: DiskConn8 = None

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
        self.disk_cap_shunt = DiskConn8(
            origin=origin,
            pars=self.qubit_params.qubit_cap_params,
            region_id=self.region_ids[0]
        )

        # SQUID is connected to the right of the disk by default (see schematics for details)
        # SQUID is further rotated according to `QubitParams.squid_connector_idx`
        qubit_finger_end = self.disk_cap_shunt.connections[0]

        # draw squid
        self.squid = AsymSquid(
            origin=qubit_finger_end,
            params=self.qubit_params.squid_params,
            region_id=self.region_ids[1],
            trans_in=DTrans.R90
        )
        connector_angles = self.qubit_params.qubit_cap_params.connector_angles
        connector_idx = self.qubit_params.squid_connector_idx
        angle = connector_angles[connector_idx]
        # shift squid to make its bottom contact pads intersected by ground plane
        # exactly through center
        outer_disk_r = self.disk_cap_shunt.pars.disk_r + self.disk_cap_shunt.pars.disk_gap
        shift_dx = outer_disk_r - self.squid.BC_list[0].center().x
        self.squid.make_trans(DCplxTrans(1, 0, False, shift_dx, 0))
        self.squid.make_trans(DCplxTrans(1, angle, False, 0, 0))

        self.primitives["squid"] = self.squid
        self.primitives["disk_cap_shunt"] = self.disk_cap_shunt

        # draw pimp in order to shorten SQUID vertical dimensions
        self.squid_pimp = CPW(
            start=origin,
            end=self.squid.TC.center(),
            width=2 * self.squid.squid_params.squid_dx,
            gap=0,
            region_id=self.region_ids[0]
        )
        self.primitives["squid_pimp"] = self.squid_pimp

        self.connections.append(origin)

    def get_connector_dr(self, connector_idx, normalized=False):
        return self.disk_cap_shunt.get_connector_dv(
            connector_idx=connector_idx, normalized=normalized
        )

    def _refresh_named_connections(self):
        self.origin = self.connections[-1]


from classLib.resonators import EMResonatorTL3QbitWormRLTail
from classLib.shapes import Donut
from classLib.helpers import rotate_around
from classLib.coplanars import DPathCPW


@dataclass()
class CqrCouplingParamsType1:
    # distance betweeen bendings of the coupling cpw and center of the qubit disks.
    bendings_disk_center_d = 320e3

    # distance between resonator's end and qubit center
    q_res_d = 500e3

    donut_delta_alpha_deg = 360 / 9 * 2 / 3
    donut_metal_width: float = 50e3
    donut_gnd_gap = 15e3
    donut_disk_d: float = 10e3

    disk1: DiskConn8 = None
    disk1_connector_idx: int = 1


@dataclass()
class ConnectivityMap:
    """
    This datastructure holds elements connections (basically,
    """
    # incidence matrix for qubits graph
    # incidence matrix entries consists of 2 numbers - corresponding
    # qubits connectors idxs (see schematics for details)
    # if any of connectors idxs is equal to `-1` then qubit pair considered disconnected
    qq_coupling_connectors_map: np.ndarray = np.zeros((16, 16, 2), dtype=int) - 1

    # squids directed from bottom to top (qubit indexes)
    direct_squids_idxs = [2, 5, 6, 9, 10, 13, 14, 15]
    # upside-down squids (qubit indexes)
    upsidedown_squids_idxs = [0, 1, 3, 4, 7, 8, 11, 12]

    # TODO: fill structure automatically for more qubits
    # horizontal
    # row 0
    qq_coupling_connectors_map[0, 1] = np.array((0, 4))
    # row 1
    qq_coupling_connectors_map[2, 3] = np.array((0, 4))
    qq_coupling_connectors_map[3, 4] = np.array((7, 4))
    # row 2
    qq_coupling_connectors_map[5, 6] = np.array((0, 3))
    qq_coupling_connectors_map[6, 7] = np.array((0, 4))
    qq_coupling_connectors_map[7, 8] = np.array((7, 4))
    # row 3
    qq_coupling_connectors_map[9, 10] = np.array((0, 3))
    qq_coupling_connectors_map[10, 11] = np.array((0, 4))
    qq_coupling_connectors_map[11, 12] = np.array((7, 4))
    # row 4
    qq_coupling_connectors_map[13, 14] = np.array((0, 3))
    qq_coupling_connectors_map[14, 15] = np.array((0, 4))

    # vertical
    # col 0
    qq_coupling_connectors_map[0, 2] = np.array((6, 2))
    qq_coupling_connectors_map[2, 5] = np.array((7, 2))
    # col 2
    qq_coupling_connectors_map[1, 3] = np.array((6, 3))
    qq_coupling_connectors_map[3, 6] = np.array((6, 2))
    qq_coupling_connectors_map[6, 9] = np.array((7, 2))
    # col 3
    qq_coupling_connectors_map[4, 7] = np.array((6, 3))
    qq_coupling_connectors_map[7, 10] = np.array((6, 2))
    qq_coupling_connectors_map[10, 13] = np.array((7, 2))
    # сol 4
    qq_coupling_connectors_map[8, 11] = np.array((6, 3))
    qq_coupling_connectors_map[11, 14] = np.array((6, 2))
    # col 5
    qq_coupling_connectors_map[12, 15] = np.array((6, 2))

    # q_idx, res_idx, q_connector_idx, ro_line_idx
    # for `q_connector_idx` see schematics .drawio
    # `-1` in element entry means that no such element exists
    # grouped in rows by direction diagonals (diagonals go from top-left to bottom-right),
    # starting with topmost diagonal corresponding to the first row.
    # ATTENTION: GROUPING ABOVE FOR CODE READABILITY ONLY. THIS STRUCTURE IS SORTED IN ASCENDING
    # ORDER BY qubit idx (i.e. pre-sorted ascending, first col) (see `self.__post_init`)
    q_res_connector_roline_map: np.ndarray = np.array(
        [
            (1, 4, 0, 0), (4, 5, 0, 0), (8, 7, 0, 0), (12, 2, 0, 0),
            (0, 1, 4, 0), (3, 6, 0, 0), (7, 0, 0, 0), (11, 3, 0, 0), (15, 7, 0, 1),
            (2, 3, 4, 1), (6, 2, 4, 1), (10, 4, 4, 1), (14, 1, 4, 1),
            (5, 0, 4, 1), (9, 5, 4, 1), (13, 6, 4, 1)

        ],
    )
    q_idx_map: np.ndarray = None

    def __post_init__(self):
        # sort by `q_idx` equivalent to sort by first entry in every row of
        # `q_res_connector_roline_map` in ascending order.
        self.q_idx_map = self.q_res_connector_roline_map[:, 0].argsort()
        self.q_res_connector_roline_map = self.q_res_connector_roline_map[self.q_idx_map]

    def get_squid_connector_idx(self, qubit_idx: int):
        """
        Returns qubit capacitor disk connector idx for squid
        Parameters
        ----------
        qubit_idx: int
            integer from `range(0,12)`

        Returns
        -------
        int
            qubit disk's connector idx
        """
        if qubit_idx in self.direct_squids_idxs:
            return 6
        elif qubit_idx in self.upsidedown_squids_idxs:
            return 2
        else:
            raise Exception(
                "ConnectivityMap.get_squid_connector_idx:"
                f"squid connector idx for qubit_idx={qubit_idx}"
                f"is not defined."
            )

    def get_md_connector_idx(self, qubit_idx):
        if qubit_idx in self.direct_squids_idxs:
            return 5
        elif qubit_idx in self.upsidedown_squids_idxs:
            return 2
        else:
            raise Exception(
                "ConnectivityMap.get_md_connector_idx:"
                f"squid connector idx for qubit_idx={qubit_idx}"
                f"is not defined."
            )


@dataclass()
class ROResonatorParams:
    """
        Static class that contains information on readout resonators geometry parameters.
        Geometry parameters has to be verified by simulation.
        """
    # see parameters details in `Design_fast.py`
    current_sim_freqs = [7.58338877, 7.34696578, 7.43494798, 7.26783703, 7.44070826,
                         7.50407024, 7.58251541, 7.50828703, 7.3511825, 7.16000063,
                         7.18802701, 7.26911275,

                         7.26911275,  7.26911275,  7.26911275,  7.26911275]
    target_freqs = [7.6, 7.36, 7.44, 7.28, 7.44, 7.52, 7.6, 7.52, 7.36, 7.2, 7.2, 7.28,
                    7.28, 7.28, 7.28, 7.28]
    target_qfactor = [10e3] * 16

    L_coupling_list = [
        1e3 * x for x in [310, 310, 310, 310] * 4
    ]
    L0_list = [986e3] * 16
    L1_list = [
        1e3 * x for x in
        [
            135.68041077160032,
            114.75000254984994,
            174.67609952692922,
            174.73092406491017,

            145.19432055474881,
            151.651946572570495,
            149.0199335546049,
            157.12727225058603,

            142.20935545352937,
            147.35206171474732,
            165.73535930161088,
            125.11764705882364,

            125.11764705882364,
            125.11764705882364,
            125.11764705882364,
            125.11764705882364
        ]
    ]
    res_r_list = [40e3] * 16  # [60e3] * 16
    tail_turn_radiuses_list = [60e3] * 16  # res_r_list
    N_coils_list = [3, 3, 3, 2,
                    3, 3, 2, 2,
                    3, 3, 2, 2,
                    3, 3, 3, 3]
    L2_list = [60e3] * 16  # res_r_list
    L3_list = []  # get numericals from Design_fast
    L4_list = [60e3] * 16  # res_r_list
    Z_res_list = [CPWParameters(width=20e3, gap=10e3)] * 16
    to_line_list = [54e3] * 16

    tail_segments_list = [[60000.0, 215000.0, 60000.0]] * 16
    res_tail_shapes_list = ["LRLRL"] * 16

    tail_turn_angles_list = [
        [np.pi / 2, -np.pi / 2],
        [np.pi / 2, -np.pi / 2],
        [np.pi / 2, -np.pi / 2],
        [np.pi / 2, -np.pi / 2],
        [np.pi / 2, -np.pi / 2],
        [np.pi / 2, -np.pi / 2],
        [np.pi / 2, -np.pi / 2],
        [np.pi / 2, -np.pi / 2],
        [np.pi / 2, -np.pi / 2],
        [np.pi / 2, -np.pi / 2],
        [np.pi / 2, -np.pi / 2],
        [np.pi / 2, -np.pi / 2],


        [np.pi / 2, -np.pi / 2],
        [np.pi / 2, -np.pi / 2],
        [np.pi / 2, -np.pi / 2],
        [np.pi / 2, -np.pi / 2]
    ]
    resonator_rotation_angles: np.ndarray = np.array(
        [
            0, 0,
            90, -45, -45,
            90, 135, -45, -45,
                135, 135, -45, 270,
                    180, 180, 180
        ],
        dtype=float
    )  # addressed by qubit index

    # corresponding to 8 RO frequencies from 7.2 to 7.76 both inclusive
    _donut_metal_width_USS = 1e3 * np.array([9.86, 11.24, 12.69, 14.20, 15.76, 17.36,
                                             17.36, 17.36])
    _donut_disk_d_USS = [20e3] * 8
    _donut_metal_width_LSS = 1e3 * np.array([52.52, 53.86, 55.08, 56.2, 57.21, 58.12,
                                             58.12, 58.12])
    _donut_disk_d_LSS = [10e3] * 8
    # see correspondence of q_idx and freq_idx in the working journal
    q_idx_ro_freq_idx = None
    q_res_coupling_params: List[CqrCouplingParamsType1] = None

    qubits_grid: QubitsGrid = None

    NQUBITS = None

    def __post_init__(self):
        self.NQUBITS = len(self.qubits_grid.pts_grid)

        # see correspondence of q_idx and freq_idx in the working journal
        # `self.q_idx_ro_freq_idx[q_idx]` contains freq idx of qubit referenced by `q_idx`
        # freq_idx instead of frequency is used, due to discrete amount of resonators with
        # different frequencies for each readout line.
        # Rules are reflected in work journal
        self.q_idx_ro_freq_idx = [
            1, 4,
            3, 6, 5,
            0, 2, 0, 7,
               5, 4, 3, 2,
                  6, 1, 7
        ]

        self.q_res_coupling_params = [None] * self.NQUBITS

        for q_idx in range(self.NQUBITS):
            freq_idx = self.q_idx_ro_freq_idx[q_idx]
            q_fl_idle = self.qubits_grid.q_fl_idle[q_idx]
            if q_fl_idle == 0:
                self.q_res_coupling_params[q_idx] = CqrCouplingParamsType1(
                    donut_metal_width=self._donut_metal_width_USS[freq_idx],
                    donut_disk_d=self._donut_disk_d_USS[freq_idx]
                )
            elif q_fl_idle == 0.5:
                self.q_res_coupling_params[q_idx] = CqrCouplingParamsType1(
                    donut_metal_width=self._donut_metal_width_LSS[freq_idx],
                    donut_disk_d=self._donut_disk_d_LSS[freq_idx]
                )

    def get_resonator_params_by_qubit_idx(self, q_idx: int):
        return {
            "Z0"                           : self.Z_res_list[q_idx],
            "L_coupling"                   : self.L_coupling_list[q_idx],
            "L0"                           : self.L0_list[q_idx],
            "L1"                           : self.L1_list[q_idx],
            "r"                            : self.res_r_list[q_idx],
            "N"                            : self.N_coils_list[q_idx],
            "tail_shape"                   : self.res_tail_shapes_list[q_idx],
            "tail_turn_radiuses"           : self.tail_turn_radiuses_list[q_idx],
            "tail_segment_lengths"         : self.tail_segments_list[q_idx],
            "tail_turn_angles"             : self.tail_turn_angles_list[q_idx],
            "tail_trans_in"                : Trans.R270,
            "to_line"                      : self.to_line_list[q_idx],
            "resonator_rotation_trans"     : DCplxTrans(
                1, self.resonator_rotation_angles[q_idx], False, 0, 0
            ),
            "additional_displacement_trans": self.get_adjusting_trans_by_qubit_idx(
                q_idx=q_idx
            )
        }

    def get_adjusting_trans_by_qubit_idx(self, q_idx: int):
        additional_displacement_trans = DCplxTrans(1, 0, False, 0, 0)
        if q_idx in [13, 14, 15]:  # bottom row
            if q_idx in [13, 14]:
                additional_displacement_trans = DCplxTrans(
                    1, 0, False,
                    DVector(-self.L_coupling_list[q_idx], 0)
                ) * additional_displacement_trans
            elif q_idx in [15]:
                additional_displacement_trans = DCplxTrans(
                    1, 0, False,
                    DVector(self.L_coupling_list[q_idx], 0)
                ) * additional_displacement_trans
        elif q_idx in [3, 4, 7, 8, 11, 12, 14, 15]:  # upper diagonal (both inner and outer)
            additional_displacement_trans = DCplxTrans(
                1, 0, False,
                DVector(self.L_coupling_list[q_idx] / 2, -self.L_coupling_list[q_idx] / 2)
            ) * additional_displacement_trans
            # align inner qubit resonators bases with outer ones
            if q_idx in [3, 7, 11]:
                additional_displacement_trans = DCplxTrans(
                    1, 0, False,
                    DVector(self.qubits_grid.dx/2, self.qubits_grid.dy/2)
                ) * additional_displacement_trans
            elif q_idx in [12]:  # rightmost column
                additional_displacement_trans = DCplxTrans(
                    1, 0, False,
                    DVector(0, self.L_coupling_list[q_idx] / 2)
                ) * additional_displacement_trans
        elif q_idx in [6, 9, 10]:  # lower diagonal (both inner and outer)
            additional_displacement_trans = DCplxTrans(
                1, 0, False,
                DVector(-self.L_coupling_list[q_idx] / 2, self.L_coupling_list[q_idx] / 2)
            ) * additional_displacement_trans
            # align inner qubit resonators bases with outer ones
            if q_idx in [6, 10]:
                additional_displacement_trans = DCplxTrans(
                    1, 0, False,
                    DVector(-self.qubits_grid.dx/2, -self.qubits_grid.dy/2)
                ) * additional_displacement_trans
        elif q_idx in [1]:
            additional_displacement_trans = DCplxTrans(
                1, 0, False,
                DVector(self.L_coupling_list[q_idx], 0)
            ) * additional_displacement_trans
        elif q_idx in [0]:
            additional_displacement_trans = DCplxTrans(
                1, 0, False,
                DVector(-self.L_coupling_list[q_idx], 0)
            ) * additional_displacement_trans

        return additional_displacement_trans


class ROResonator(EMResonatorTL3QbitWormRLTail):
    def __init__(
        self, Z0, L_coupling, L0, L1, r, N,
        tail_shape, tail_turn_radiuses,
        tail_segment_lengths, tail_turn_angles, tail_trans_in=None,
        coupling_pars: CqrCouplingParamsType1 = CqrCouplingParamsType1(),
        to_line: float = None,
        resonator_rotation_trans=DCplxTrans(1, 0, False, 0, 0),
        additional_displacement_trans=DCplxTrans(1, 0, False, 0, 0),
        trans_in=None
    ):
        """
        Resonator is placed at origin and then moved to their corresponding qubit, with coupling
        parameters drawn.

        1. Draw resonator utilizing its base class
        2. Rotate around `resonator_rotation_trans`
        3. Translate resonator such that `self.end` is at the coupling direction of a qubit
        at a distance `self.coupling_pars.q_res_d` from the qubit's disk center.
        4. Make additional displacement descripted in `additional_displacement_trans` of a
        resonator to provide enough space for control lines to reach a qubit.
        5. Draw donut sector coupling at the qubit's connector idx
        from `coupling_pars.disk1_connector_idx`
        6. Draw coplanar waveguide that connects end of a resonator to center of a donut sector.
        7. Patch a bandage such that donut sector and cpw connecting it to resonator are
        connected robustly.


        Parameters
        ----------
        Z0
        L_coupling
        L0
        L1
        r
        N
        tail_shape
        tail_turn_radiuses
        tail_segment_lengths
        tail_turn_angles
        tail_trans_in
        coupling_pars
        to_line
        resonator_rotation_trans: DCplxTrans
            rotation of the resonator itself
        additional_displacement_trans: DCplxTrans
            additional displacement, after resonators body is put near it's corresponding qubit
        trans_in
        """
        self.coupling_pars = coupling_pars
        self.resonator_rotation_trans = resonator_rotation_trans
        self.additional_displacement_trans = additional_displacement_trans
        self.to_line = to_line
        super().__init__(
            Z0, start=DPoint(0, 0),
            L_coupling=L_coupling,
            L0=L0, L1=L1, r=r, N=N,
            tail_shape=tail_shape, tail_turn_radiuses=tail_turn_radiuses,
            tail_segment_lengths=tail_segment_lengths, tail_turn_angles=tail_turn_angles,
            tail_trans_in=tail_trans_in,
            trans_in=trans_in
        )
        self._geometry_parameters["donut_metal_width, um"] = coupling_pars.donut_metal_width / 1e3
        self._geometry_parameters["donut_disk_d, um"] = coupling_pars.donut_disk_d / 1e3
        self._geometry_parameters["donut_gnd_gap, um"] = coupling_pars.donut_gnd_gap / 1e3
        self._geometry_parameters["donut_delta_alpha_deg, deg"] = \
            coupling_pars.donut_delta_alpha_deg

    def init_primitives(self):
        super().init_primitives()

        ''' move resonator to it's corresponding qubit. Resonator is properly rotated on the way '''
        qubit_res_d = self.coupling_pars.q_res_d
        self.make_trans(self.resonator_rotation_trans)
        dv = self.start - self.end + self.coupling_pars.disk1.origin + \
             self.resonator_rotation_trans * DVector(0, qubit_res_d)
        self.make_trans(DCplxTrans(1, 0, False, dv))
        self.make_trans(self.additional_displacement_trans)

        """ add fork to the end of the base class resonator """
        # adding fork horizontal part
        self.draw_cpw_path_arc_coupling_end()

        # remove open end from the base class resonator
        del self.primitives["cpw_end_open_gap"]
        del self.cpw_end_open_gap

    def draw_cpw_path_arc_coupling_end(self):
        angle1 = self.coupling_pars.disk1.angle_connections[
                     self.coupling_pars.disk1_connector_idx
                 ] / (2 * np.pi) * 360

        self.arc_coupler = Donut(
            origin=self.coupling_pars.disk1.origin,
            inner_r=self.coupling_pars.donut_disk_d + self.coupling_pars.disk1.pars.disk_r,
            ring_width=self.coupling_pars.donut_metal_width,
            alpha_start=-self.coupling_pars.donut_delta_alpha_deg / 2,
            alpha_end=self.coupling_pars.donut_delta_alpha_deg / 2,
            region_id=self.region_id
        )
        rotate_around(self.arc_coupler, self.arc_coupler.origin, angle1)
        d_alpha_deg = 360 / 2 / np.pi * self.coupling_pars.donut_gnd_gap / (
                (self.arc_coupler.inner_r +
                 self.arc_coupler.outer_r)
                / 2)
        self.arc_coupler_empty = Donut(
            origin=self.arc_coupler.origin,
            inner_r=self.arc_coupler.inner_r - self.coupling_pars.donut_disk_d,
            outer_r=self.arc_coupler.outer_r + self.coupling_pars.donut_gnd_gap,
            alpha_start=self.arc_coupler.alpha_start - d_alpha_deg,
            alpha_end=self.arc_coupler.alpha_end + d_alpha_deg,
            inverse=True,
            region_id=self.region_id
        )
        rotate_around(self.arc_coupler_empty, self.arc_coupler_empty.origin, angle1)
        # empty first, filling metal later
        self.primitives["arc_coupler_empty"] = self.arc_coupler_empty
        self.primitives["arc_coupler"] = self.arc_coupler

        resonator_end = self.end
        resonator_end_dv_n = DVector(np.cos(self.alpha_end), np.sin(self.alpha_end))
        # 1*self.r restricts first bending angle to be < 45 deg (see DPathCPW class, bendings
        # drawing section)
        resonator_end_bending_pt = resonator_end + 1 * self.r * resonator_end_dv_n

        connector_dv_n = DVector(np.cos(np.pi * angle1 / 180), np.sin(np.pi * angle1 / 180))
        disk_far_bending_point = self.coupling_pars.disk1.origin + \
                                 self.coupling_pars.bendings_disk_center_d * connector_dv_n
        # print("qubit center:", self.coupling_pars.disk1.origin)
        # print("coupling line arc angle:", angle1, "deg")
        # print((disk_far_bending_point - resonator_end).abs())
        # print((self.arc_coupler.outer_arc_center - disk_far_bending_point).abs())
        # print(self.r)
        # print()
        self.res_coulingArc_cpw_path = DPathCPW(
            points=[resonator_end, resonator_end_bending_pt,
                    disk_far_bending_point, self.arc_coupler.outer_arc_center],
            cpw_parameters=[self.Z0],
            turn_radii=[self.r],
            region_id=self.region_id
        )
        self.primitives["res_coulingArc_cpw_path"] = self.res_coulingArc_cpw_path

        self.coupling_bandage = CPW(
            start=self.res_coulingArc_cpw_path.end,
            end=(self.arc_coupler.inner_arc_center + self.arc_coupler.outer_arc_center) / 2,
            width=self.Z0.width,
            gap=0,
            region_id=self.region_id
        )
        self.primitives["coupling_bandage"] = self.coupling_bandage

    def get_resonator_ro_connections(self):
        p_res_start = self.start + \
                      self.resonator_rotation_trans * DVector(-3 * self.r, self.to_line)
        p_res_end = (self.start +
                     self.resonator_rotation_trans * DVector(self.L_coupling + 3 * self.r, self.to_line))
        return [p_res_start, p_res_end]
