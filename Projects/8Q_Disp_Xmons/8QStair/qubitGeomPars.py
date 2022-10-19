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

from pya import DCplxTrans

# import self-made API
import classLib

reload(classLib)

from classLib.baseClasses import ComplexBase
from classLib.coplanars import CPW
from classLib.josJ import AsymSquidParams, AsymSquid


# project specific files
from globalDefinitions import CHIP

@dataclass()
class QubitsGrid:
    # in fractions of chip dimensions
    origin: float = DVector(CHIP.dx / 2, CHIP.dy / 2)
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
    disk_r = 200e3
    disk_gap = 20e3
    pimp_l = 0e3
    conn_width = 0e3
    conn_side_gap = 0e3
    conn_front_gap = 0e3


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
        # print(self.origin)
        origin = DPoint(0, 0)

        from classLib.shapes import Disk
        self.empty_disk = Disk(
            center=origin, r=self.pars.disk_r + self.pars.disk_gap,
            region_id=self.region_id, inverse=True
        )
        self.primitives["empty_disk"] = self.empty_disk

        # draw star-like connection flanges
        angles = np.linspace(0, 360, 8, endpoint=False)  # degree
        for i, angle in enumerate(angles):
            trans = DCplxTrans(1, angle, False, 0, 0)
            # finger with side gap
            cpw_l = self.pars.disk_r + self.pars.pimp_l
            cpw = CPW(
                start=origin,
                end=origin + DVector(cpw_l, 0),
                width=self.pars.conn_width,
                gap=self.pars.conn_side_gap,
                open_end_l=self.pars.conn_front_gap,
                trans_in=trans,
                region_id=self.region_id
            )
            self.conn8_list.append(cpw)
            self.primitives["conn" + str(i)] = cpw

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

        qubit_finger = self.cap_shunt.conn8_list[6]

        # draw squid
        squid_pars = self.qubit_params.squid_params
        # vertical shift of every squid local origin coordinates
        # this line puts squid center on the
        # "outer capacitor plate's edge" of the shunting capacitance.
        squid_center = qubit_finger.open_end_end
        # next step is to make additional shift such that top of the
        # BC polygon is located at the former squid center position with
        # additional `_shift_into_substrate = 1.5e3` um shift in direction of substrate
        _shift_to_BC_center = squid_pars.shadow_gap / 2 + \
                              squid_pars.SQLBT_dy + squid_pars.SQB_dy + \
                              squid_pars.BCW_dy + \
                              self._shift_into_substrate
        squid_center += DVector(0, _shift_to_BC_center)

        self.squid = AsymSquid(
            origin=squid_center,
            params=self.qubit_params.squid_params,
            region_id=self.region_ids[1]
        )

        self.primitives["squid"] = self.squid
        self.primitives["cap_shunt"] = self.cap_shunt