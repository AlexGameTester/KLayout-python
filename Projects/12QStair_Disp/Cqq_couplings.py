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
from classLib.helpers import rotate_around
from classLib.chipDesign import ChipDesign
from classLib.contactPads import ContactPad
from classLib.shapes import Donut

from designElementsGeometry import DiskConn8, DiskConn8Pars


# TODO: use `Fork` and `Meander` to implement new, simpler, curved resonator drawing class.
@dataclass()
class TwoTeethForkParams:
    # TODO: class annotations for all this file
    # fork's base
    metal_width: float = 10e3
    foundation_length = 46e3
    gnd_gap = 15e3

    # fork's teeth
    teeth_metal_width = None
    teeth_length = 70e3
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
        self.foundation_cpw: CPW = None
        self.upper_tooth_cpw: CPW = None
        self.lower_tooth_cpw: CPW = None
        self.params: TwoTeethForkParams = params
        super().__init__(
            origin=origin, trans_in=trans_in, region_id=region_id,
            region_ids=["ph"]
        )

    def init_primitives(self):
        origin = DPoint(0, 0)

        # both teeth are looking towards Ox direction by default
        # and foundation of the fork's teeth is adjacent to Oy axis from right-hand side.
        self.foundation_cpw = CPW(
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
        self.primitives["foundation_cpw"] = self.foundation_cpw

        # upper fork tooth
        # foundation goes from bottom to top, hence start is at the end of
        # the foundation of the fork
        start = self.foundation_cpw.end + DVector(
            self.params.metal_width / 2,
            -self.params.teeth_metal_width / 2
        )
        self.upper_tooth_cpw = CPW(
            width=self.params.teeth_metal_width,
            gap=self.params.teeth_gnd_gap,
            start=start,
            end=start + DVector(self.params.teeth_length, 0),
            open_end_gap=self.params.teeth_gnd_gap,
            region_id=self.region_id
        )
        self.primitives["upper_tooth_cpw"] = self.upper_tooth_cpw
        #
        # lower fork cpw
        start = self.foundation_cpw.start + DVector(
            self.params.metal_width / 2,
            self.params.teeth_metal_width / 2
        )
        self.lower_tooth_cpw = CPW(
            width=self.params.teeth_metal_width,
            gap=self.params.teeth_gnd_gap,
            start=start,
            end=start + DVector(self.params.teeth_length, 0),
            open_end_gap=self.params.teeth_gnd_gap,
            region_id=self.region_id
        )
        self.primitives["lower_tooth_cpw"] = self.lower_tooth_cpw


@dataclass()
class CqqCouplingParamsType2:
    # central 2 sided fork with 2 teeth
    length_without_forks = 500e3
    bendings_disk_d = 280e3
    fork_params = TwoTeethForkParams()

    # finger growing from qubit's disk
    finger_extension_l = 110e3
    finger_metal_width = 10e3
    finger_gnd_gap = 15e3

    finger_fork_d = 10e3

    disk1: DiskConn8 = None
    disk1_connector_idx: int = 1
    disk2: DiskConn8 = None
    disk2_connector_idx: int = 1


class CqqCouplingType2(ComplexBase):
    # TODO: document this class properly if in production more than once - assigned to @Shamil777
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
        origin = (self.params.disk1.origin + self.params.disk2.origin) / 2
        # declare helping structures
        q12_dv = self.params.disk2.origin - self.params.disk1.origin
        q12_dv /= q12_dv.abs()
        angle = np.arctan2(q12_dv.y, q12_dv.x)

        ### DRAW disk fingers ###
        angle1 = self.params.disk1.angle_connections[self.params.disk1_connector_idx]
        self.finger1 = CPW(
            start=self.params.disk1.origin + DVector(self.params.disk1.pars.disk_r, 0),
            end=self.params.disk1.origin + \
                DVector((self.params.finger_extension_l + self.params.disk1.pars.disk_r), 0),
            width=self.params.finger_metal_width,
            gap=self.params.finger_gnd_gap,
            open_end_gap=self.params.finger_gnd_gap,
            region_id=self.region_id
        )
        rotate_around(
            primitive=self.finger1,
            rotation_center=self.params.disk1.origin, angle_deg=angle1 * 180 / np.pi
        )
        self.primitives["finger1"] = self.finger1

        angle2 = self.params.disk2.angle_connections[self.params.disk2_connector_idx]
        self.finger2 = CPW(
            start=self.params.disk2.origin + DVector(self.params.disk2.pars.disk_r, 0),
            end=self.params.disk2.origin + \
                DVector((self.params.finger_extension_l + self.params.disk2.pars.disk_r), 0),
            width=self.params.finger_metal_width,
            gap=self.params.finger_gnd_gap,
            open_end_gap=self.params.finger_gnd_gap,
            region_id=self.region_id
        )
        rotate_around(
            primitive=self.finger2,
            rotation_center=self.params.disk2.origin, angle_deg=angle2 * 180 / np.pi
        )
        self.primitives["finger2"] = self.finger2

        # fork 1
        finger1_dv_n = self.finger1.dr/self.finger1.dr.abs()
        # origin and transformation
        self.fork1 = TwoTeethFork(
            origin=self.finger1.end + (self.params.finger_fork_d +
                                       self.params.fork_params.metal_width)*finger1_dv_n,
            params=self.params.fork_params,
            region_id=self.region_id,
            trans_in=DCplxTrans(1, angle1*180/np.pi + 180, False, 0, 0)
        )
        self.primitives["fork1"] = self.fork1

        # fork 2
        finger2_dv_n = self.finger2.dr / self.finger2.dr.abs()
        self.fork2 = TwoTeethFork(
            origin=self.finger2.end + (self.params.finger_fork_d +
                                       self.params.fork_params.metal_width)*finger2_dv_n,
            params=self.params.fork_params,
            region_id=self.region_id,
            trans_in=DCplxTrans(1, angle2*180/np.pi + 180, False, 0, 0)
        )
        self.primitives["fork2"] = self.fork2

        ### DRAW Cqq coupler central part ###
        p_start = self.fork1.origin
        p1 = self.finger1.start + self.params.bendings_disk_d*finger1_dv_n
        p2 = self.finger2.start + self.params.bendings_disk_d*finger2_dv_n
        p_end = self.fork2.origin
        cpw_central = DPathCPW(
            points=[p_start, p1, p2, p_end],
            cpw_parameters=[
                CPWParameters(
                    width=self.params.fork_params.metal_width,
                    gap=self.params.fork_params.gnd_gap
                )
            ],
            turn_radii=[self.params.bendings_disk_d/5],
            region_id=self.region_id
        )
        self.primitives["cpw_central"] = cpw_central

        '''
        eliminate gap between finger and circle
        and draw finger, if fork tooth `self.params.fork_params.gnd_gap` is large
        such that tooth's ground gap erased fingers
        '''
        self.finger1_bandage = CPW(
            start=self.params.disk1.origin,
            end=self.finger1.end,
            width=self.params.finger_metal_width,
            gap=0,
            region_id=self.region_id
        )
        self.primitives["coupling1_bandage"] = self.finger1_bandage

        self.finger2_bandage = CPW(
            start=self.params.disk2.origin,
            end=self.finger2.end,
            width=self.params.finger_metal_width,
            gap=0,
            region_id=self.region_id
        )
        self.primitives["coupling2_bandage"] = self.finger2_bandage

        self.connections.append(origin)

    def _refresh_named_connections(self):
        # TODO: make `self.origin` the default point which coordinates follows all the
        #  transformations of the object it attributes to
        self.origin = self.connections[-1]


@dataclass()
class CqqCouplingParamsType1:
    # central 2 sided fork with 2 teeth
    length_without_forks = 500e3

    # distance betweeen bendings of the coupling cpw and center of the qubit disks.
    bendings_disk_center_d = 280e3

    central_metal_width = 40e3
    central_gnd_gap = 10e3

    donut_delta_alpha_deg = 360/9 * 2/3
    donut_metal_width = 40e3
    donut_gnd_gap = 15e3
    donut_disk_d = 10e3

    disk1: DiskConn8 = None
    disk1_connector_idx: int = 1
    disk2: DiskConn8 = None
    disk2_connector_idx: int = 1


class CqqCouplingType1(ComplexBase):
    # TODO: document this class properly if in production more than once - assigned to @Shamil777
    def __init__(
        self, origin, params:
        CqqCouplingParamsType1 = CqqCouplingParamsType1(),
        trans_in=None, region_id="ph", postpone_drawing=False, region_ids=["ph"]
    ):
        self.params: CqqCouplingParamsType1 = params
        super().__init__(
            origin.dup(), trans_in=trans_in, region_id=region_id, postpone_drawing=postpone_drawing,
            region_ids=region_ids
        )

    def init_primitives(self):
        origin = (self.params.disk1.origin + self.params.disk2.origin) / 2
        # declare helping structures
        q12_dv = self.params.disk2.origin - self.params.disk1.origin
        q12_dv /= q12_dv.abs()
        angle = np.arctan2(q12_dv.y, q12_dv.x)

        donut_delta_alpha_deg = self.params.donut_delta_alpha_deg

        ### DRAW disk fingers ###
        angle1 = self.params.disk1.angle_connections[self.params.disk1_connector_idx]/(2*np.pi)*360
        self.arc_coupler1 = Donut(
            origin=self.params.disk1.origin,
            inner_r=self.params.donut_disk_d + self.params.disk1.pars.disk_r,
            ring_width=self.params.donut_metal_width,
            alpha_start=-donut_delta_alpha_deg/2,
            alpha_end=donut_delta_alpha_deg/2,
            region_id=self.region_id
        )
        rotate_around(self.arc_coupler1, self.arc_coupler1.origin, angle1)
        # calculate angle required for gnd gap at the faces to be proper
        d_alpha_deg = 360/2/np.pi * self.params.donut_gnd_gap/((self.arc_coupler1.inner_r +
                                                               self.arc_coupler1.outer_r)/2)
        self.arc_coupler1_empty = Donut(
            origin = self.arc_coupler1.origin,
            inner_r=self.arc_coupler1.inner_r - self.params.donut_disk_d,
            outer_r=self.arc_coupler1.outer_r + self.params.donut_gnd_gap,
            alpha_start=self.arc_coupler1.alpha_start - d_alpha_deg,
            alpha_end=self.arc_coupler1.alpha_end + d_alpha_deg,
            inverse=True,
            region_id=self.region_id
        )
        rotate_around(self.arc_coupler1_empty, self.arc_coupler1_empty.origin, angle1)
        # empty first, filling metal later
        self.primitives["arc_coupler1_empty"] = self.arc_coupler1_empty
        self.primitives["arc_coupler1"] = self.arc_coupler1

        angle2 = self.params.disk2.angle_connections[self.params.disk2_connector_idx]/(2*np.pi)*360
        self.arc_coupler2 = Donut(
            origin=self.params.disk2.origin,
            inner_r=self.params.donut_disk_d + self.params.disk2.pars.disk_r,
            ring_width=self.params.donut_metal_width,
            alpha_start=-donut_delta_alpha_deg / 2,
            alpha_end=donut_delta_alpha_deg / 2,
            region_id=self.region_id
        )
        rotate_around(self.arc_coupler2, self.arc_coupler2.origin, angle2)
        # calculate angle required for gnd gap at the faces to be proper
        d_alpha_deg = 360 / 2 / np.pi * self.params.donut_gnd_gap / ((self.arc_coupler2.inner_r +
                                                                      self.arc_coupler2.outer_r)
                                                                     / 2)
        self.arc_coupler2_empty = Donut(
            origin=self.arc_coupler2.origin,
            inner_r=self.arc_coupler2.inner_r - self.params.donut_disk_d,
            outer_r=self.arc_coupler2.outer_r + self.params.donut_gnd_gap,
            alpha_start=self.arc_coupler2.alpha_start - d_alpha_deg,
            alpha_end=self.arc_coupler2.alpha_end + d_alpha_deg,
            inverse=True,
            region_id=self.region_id
        )
        rotate_around(self.arc_coupler2_empty, self.arc_coupler2_empty.origin, angle2)
        # empty first, filling metal later
        self.primitives["arc_coupler2_empty"] = self.arc_coupler2_empty
        self.primitives["arc_coupler2"] = self.arc_coupler2


        ### DRAW Cqq coupler central part ###
        dv1_n = self.arc_coupler1.outer_arc_center - self.arc_coupler1.origin
        dv1_n = dv1_n/dv1_n.abs()
        dv2_n = self.arc_coupler2.outer_arc_center - self.arc_coupler2.origin
        dv2_n = dv2_n/dv2_n.abs()
        p_start = self.arc_coupler1.outer_arc_center
        p1 = self.arc_coupler1.origin + self.params.bendings_disk_center_d * dv1_n
        p2 = self.arc_coupler2.origin + self.params.bendings_disk_center_d * dv2_n
        p_end = self.arc_coupler2.outer_arc_center
        self.cpw_central = DPathCPW(
            points=[p_start, p1, p2, p_end],
            cpw_parameters=[
                CPWParameters(
                    width=self.params.central_metal_width,
                    gap=self.params.central_gnd_gap
                )
            ],
            turn_radii=[self.params.bendings_disk_center_d / 5],
            region_id=self.region_id
        )
        self.primitives["cpw_central"] = self.cpw_central

        '''
        eliminate gap between donut sector and coupling stick
        and draw finger, if fork tooth `self.params.fork_params.gnd_gap` is large
        '''
        self.coupling1_bandage = CPW(
            start=self.cpw_central.start,
            end=(self.arc_coupler1.inner_arc_center + self.arc_coupler1.outer_arc_center)/2,
            width=self.params.central_metal_width,
            gap=0,
            region_id=self.region_id
        )
        self.primitives["coupling1_bandage"] = self.coupling1_bandage

        self.coupling2_bandage = CPW(
            start=self.cpw_central.end,
            end=(self.arc_coupler2.inner_arc_center + self.arc_coupler2.outer_arc_center) / 2,
            width=self.params.central_metal_width,
            gap=0,
            region_id=self.region_id
        )
        self.primitives["coupling2_bandage"] = self.coupling2_bandage

        self.connections.append(origin)

    def _refresh_named_connections(self):
        self.origin = self.connections[-1]