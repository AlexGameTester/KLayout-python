# import built-ins
from typing import List
import os
import itertools
import sys
from dataclasses import dataclass

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
from classLib.baseClasses import ComplexBase
from classLib.coplanars import CPW, CPWParameters, DPathCPW
from classLib.chipDesign import ChipDesign
from classLib.contactPads import ContactPad

from designElementsGeometry import DiskConn8, DiskConn8Pars


# TODO: use `Fork` and `Meander` to implement new, simpler, curved resonator drawing class.
@dataclass()
class TwoTeethForkParams:
    # TODO: class annotations for all this file
    # fork's base
    metal_width: float = 20e3
    foundation_length = 100e3
    gnd_gap = 25e3

    # fork's teeth
    teeth_metal_width = None
    teeth_length = 35e3
    teeth_gnd_gap = None

    def __post_init__(self):
        if self.teeth_metal_width is None:
            self.teeth_metal_width = self.metal_width
        if self.teeth_gnd_gap is None:
            self.teeth_gnd_gap = self.gnd_gap


class TwoTeethFork(ComplexBase):
    # TODO: figure how to handle `region_id` and `region_ids` with complex base class
    # TODO: document this shape - assigned to @Shamil777
    def __init__(
        self, origin, params: TwoTeethForkParams = TwoTeethForkParams(), trans_in=None,
        region_id="ph"
    ):
        self.params: TwoTeethForkParams = params
        super().__init__(
            origin=origin, trans_in=trans_in, region_id=region_id,
            region_ids=["ph"]
        )

    def init_primitives(self):
        origin = DPoint(0, 0)

        # both teeth are looking towards Ox direction by default
        # and foundation of the fork's teeth is adjacent to Oy axis from right-hand side.
        foundation_cpw = CPW(
            width=self.params.metal_width,
            gap=self.params.gnd_gap,
            start=origin + DVector(
                self.params.metal_width / 2,
                -self.params.foundation_length / 2
            ),
            end=origin + DVector(
                self.params.metal_width / 2,
                self.params.foundation_length / 2
            ),
            open_end_gap=self.params.gnd_gap,
            open_start_gap=self.params.gnd_gap,
            region_id=self.region_id

        )
        self.primitives["foundation_cpw"] = foundation_cpw

        # upper fork tooth
        # foundation goes from bottom to top, hence start is at the end of
        # the foundation of the fork
        start = foundation_cpw.end + DVector(
            self.params.metal_width/2,
            -self.params.teeth_metal_width / 2
        )
        upper_fork_cpw = CPW(
            width=self.params.teeth_metal_width,
            gap=self.params.teeth_gnd_gap,
            start=start,
            end=start + DVector(self.params.teeth_length, 0),
            open_end_gap=self.params.teeth_gnd_gap,
            region_id=self.region_id
        )
        self.primitives["upper_fork_cpw"] = upper_fork_cpw
        #
        # lower fork cpw
        start = foundation_cpw.start + DVector(
            self.params.metal_width/2,
            self.params.teeth_metal_width / 2
        )
        lower_fork_cpw = CPW(
            width=self.params.teeth_metal_width,
            gap=self.params.teeth_gnd_gap,
            start=start,
            end=start + DVector(self.params.teeth_length, 0),
            open_end_gap=self.params.teeth_gnd_gap,
            region_id=self.region_id
        )
        self.primitives["lower_fork_cpw"] = lower_fork_cpw


@dataclass()
class CqqCouplingParamsType2:
    # central 2 sided fork with 2 teeth
    gnd_gap = 20e3
    length_without_forks = 600e3
    forks_length = 200e3
    metal_width = 20e3

    # finger growing from qubit's disk
    finger_extension_l = 50e3
    finger_metal_width = 20e3
    finger_gnd_gap = 20e3

    disk1: DiskConn8 = None
    disk1_connector_idx = 0
    disk2: DiskConn8 = None
    disk2_connector_idx = 1


class CqqCouplingType2(ComplexBase):
    # TODO: document this class properly if in production even once - assigned to @Shamil777
    def __init__(
        self, origin, params:
        CqqCouplingParamsType2 = CqqCouplingParamsType2(),
        trans_in=None, region_id="ph", postpone_drawing=False, region_ids=["ph"]
    ):
        self.params: CqqCouplingParamsType2 = params
        super().__init__(
            origin.dup(), trans_in=trans_in, region_id=region_id, postpone_drawing=postpone_drawing,
            region_ids=region_ids
        )

    def init_primitives(self):
        origin = (self.params.disk1.origin + self.params.disk2.origin)/2

        # declare helping structures
        q12_dv = self.params.disk2.origin - self.params.disk1.origin
        q12_dv /= q12_dv.abs()
        angle = np.arctan2(q12_dv.y, q12_dv.x)
        rotation = DCplxTrans(1, angle*180/np.pi, False, 0, 0)


        # ### DRAW Cqq coupler central part ###
        # central horizontal cpw
        p1 = (self.params.disk1.origin + self.params.disk2.origin) / 2 - \
             self.params.length_without_forks / 2 * q12_dv
        p2 = p1 + self.params.length_without_forks * q12_dv

        # right fork #2
        right_fork = TwoTeethFork(
            origin=p2,
            region_id=self.region_id,
            trans_in=rotation
        )
        self.primitives["right_fork"] = right_fork

        left_fork = TwoTeethFork(
            origin=p1,
            region_id=self.region_id,
            trans_in=DCplxTrans(1, 180, False, 0, 0) * rotation
        )
        self.primitives["left_fork"] = left_fork

        cpw_central = CPW(
            start=p1, end=p2, width=self.params.metal_width, gap=self.params.gnd_gap,
            region_id=self.region_id
        )
        self.primitives["cpw_central"] = cpw_central

        finger1 = CPW(
            start=self.params.disk1.origin + self.params.disk1.pars.disk_r * q12_dv,
            end=self.params.disk1.origin + \
                (self.params.finger_extension_l + self.params.disk1.pars.disk_r) * q12_dv,
            width=self.params.finger_metal_width,
            gap=self.params.finger_gnd_gap,
            region_id=self.region_id
        )

        ### DRAW disk fingers ###
        self.primitives["finger1"] = finger1
        # eliminate gap between finger and circle
        finger1_bandage = CPW(
            start=finger1.start,
            end=self.params.disk1.origin,
            width=self.params.finger_metal_width,
            gap=0,
            region_id=self.region_id
        )
        self.primitives["finger1_bandage"] = finger1_bandage

        finger2 = CPW(
            start=self.params.disk2.origin - self.params.disk2.pars.disk_r * q12_dv,
            end=self.params.disk2.origin - \
                (self.params.finger_extension_l + self.params.disk2.pars.disk_r) * q12_dv,
            width=self.params.finger_metal_width,
            gap=self.params.finger_gnd_gap,
            region_id=self.region_id
        )
        self.primitives["finger2"] = finger2
        finger2_bandage = CPW(
            start=finger2.start,
            end=self.params.disk2.origin,
            width=self.params.finger_metal_width,
            gap=0,
            region_id=self.region_id
        )
        self.primitives["finger2_bandage"] = finger2_bandage

        self.connections.append(origin)

    def _refresh_named_connections(self):
        self.origin = self.connections[-1]
