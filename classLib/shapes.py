# import built-ins
from typing import List, Union

# import good 3rd party
import numpy as np
# import project specific 3rd party
from pya import DPoint, DSimplePolygon, SimplePolygon, \
    DPolygon, Polygon, Region, DPath, DVector, DBox

# import project lib
from classLib.baseClasses import ElementBase


class Rectangle(ElementBase):
    def __init__(self, origin, width, height, trans_in=None,
                 inverse=False, region_id="default", ):
        """
        Represents rectanlge object

        Parameters
        ----------
        origin : DPoint
            left-bottom point of rectangle
        width : float
            rectangle's width
        height : b
            rectangle's b
        trans_in : DcplxTrans
            Initial transformation
        inverse : bool
            If the rectangle should be added or substracted from destination
            during the call of `self.place`.
            True - rectangle is substracted
            False - rectangle is added (default)
        """
        self.width = width
        self.height = height
        self.p1 = origin
        super().__init__(
            origin=origin, trans_in=trans_in,
            inverse=inverse, region_id=region_id
        )
        self.p1 = self.connections[0]
        self.p2 = self.connections[1]

    def init_regions(self):
        origin = DPoint(0, 0)
        p1 = origin + DPoint(self.width, 0)
        self.p2 = p1 + DPoint(0, self.height)
        p3 = self.p2 + DPoint(-self.width, 0)
        pts_arr = [origin, p1, self.p2, p3]
        if self.inverse:
            self.empty_region.insert(
                SimplePolygon(DSimplePolygon(pts_arr)))
        else:
            self.metal_region.insert(
                SimplePolygon(DSimplePolygon(pts_arr))
            )
        self.connections = [origin, self.p2]

    def _refresh_named_connections(self):
        self.p1 = self.connections[0]
        self.p2 = self.connections[1]

    def center(self):
        return (self.p1 + self.p2) / 2


class Cross(ElementBase):
    def __init__(self,
                 origin,
                 inner_square_a,
                 outer_square_a,
                 trans_in=None,
                 inverse=False,
                 region_id="default"):
        """

        Parameters
        ----------
        origin : DPoint
        inner_square_a : float
            width of cross
        outer_square_a : float
            length of square that encloses the cross
            2*cross_width + cross_branch_len
        trans_in : Union[DCplxTrans, CplxTrans, DTrans, Trans]
            Represent object's local coordinate system transformation
            see more at
            https://www.klayout.de/doc-qt4/code/class_DCplxTrans.html
        inverse : bool
            Whether to interchange `empty_polygon` and
            `metal_polygon` before placing
        region_id : str
            region identifier to object if object is part of a compound
            object and/or resides in objects tree.
        """
        self.in_a = inner_square_a
        self.out_a = outer_square_a
        self.center = origin + DPoint(self.out_a / 2, self.out_a / 2)
        super().__init__(
            origin=origin, trans_in=trans_in,
            inverse=inverse, region_id=region_id
        )
        self.center = self.connections[0]

    def init_regions(self):
        origin = DPoint(0, 0)
        w = self.out_a / 2 - self.in_a / 2

        rec1 = Rectangle(origin, w, w)
        p2 = origin + DPoint(self.in_a + w, 0)
        rec2 = Rectangle(p2, w, w)
        p3 = origin + DPoint(self.in_a + w, self.in_a + w)
        rec3 = Rectangle(p3, w, w)
        p4 = origin + DPoint(0, self.in_a + w)
        rec4 = Rectangle(p4, w, w)

        tmp_reg = Region()

        rec1.place(tmp_reg)
        rec2.place(tmp_reg)
        rec3.place(tmp_reg)
        rec4.place(tmp_reg)

        rec = Rectangle(origin, self.out_a, self.out_a)
        rec.place(self.metal_region)

        self.empty_region = tmp_reg
        self.connections = [self.center]


class Disk(ElementBase):
    def __init__(self, center, r, n_pts=50, inverse=False,
                 offset_angle=0, trans_in=None, region_id="default"):
        """
        Represents disk of radius `r` and `n_pts` points constituting
        Parameters
        ----------
        center : DPoint
            center of the circle
        r : float
            circle radius
        trans_in : DcplxTrans
            initial transformation, None by default
        n_pts : int
            number of points comprising the circumference of the ring (50 by default)
        inverse : bool
            if True then the ring is subtracted from width layer (False by default)
        offset_angle : float
            Angle in radians where the first point of the circle will be placed.
            Makes sense for small `n_pts` values.
        trans_in : Union[DCplxTrans, CplxTrans, DTrans, Trans]
            Represent object's local coordinate system transformation
            see more at
            https://www.klayout.de/doc-qt4/code/class_DCplxTrans.html
        inverse : bool
            Whether to interchange `empty_polygon` and
            `metal_polygon` before placing
        region_id : str
            region identifier to object if object is part of a compound
            object and/or resides in objects tree.
        """
        self.center = center
        self._offset_angle = offset_angle
        self.r = r
        self.n_pts = n_pts
        super().__init__(origin=center, trans_in=trans_in, inverse=inverse,
                         region_id=region_id)

    def init_regions(self):
        origin = DPoint(0, 0)
        alphas = np.linspace(0, 2 * np.pi, self.n_pts + 1)
        dpts_arr = self.r * np.transpose([np.cos(alphas), np.sin(alphas)])
        dpts_arr = [DPoint(x, y) for x, y in dpts_arr]
        circle_polygon = SimplePolygon(DSimplePolygon(dpts_arr))
        if self.inverse:
            self.empty_region.insert(circle_polygon)
        else:
            self.metal_region.insert(circle_polygon)

        self.connections.extend([origin, origin + DVector(0, -self.r)])
        self.angle_connections.extend([0, 0])

    def _refresh_named_connections(self):
        self.center = self.connections[0]


class DiskSector(ElementBase):
    def __init__(self, center, r, alpha_start=0, alpha_end=np.pi,
                 n_pts=50, trans_in=None, inverse=False,
                 region_id="default"):
        """

        Parameters
        ----------
        center : DPoint
            center of the disk
        r : float
            disk radius
        alpha_start : float
            starting angle
        alpha_end : float
            ending angle
        n_pts : int
            number of consecutive fixed step-size points array
            constituting disk perimeter
        trans_in : Union[DCplxTrans, CplxTrans, DTrans, Trans]
            Represent object's local coordinate system transformation
            see more at
            https://www.klayout.de/doc-qt4/code/class_DCplxTrans.html
        inverse : bool
            Whether to interchange `empty_polygon` and
            `metal_polygon` before placing
        region_id : str
            region identifier to object if object is part of a compound
            object and/or resides in objects tree.
        """
        self.center = center
        self.r = r
        self.alpha_start = alpha_start
        self.alpha_end = alpha_end
        self.n_pts = n_pts
        super().__init__(
            origin=center, trans_in=trans_in, inverse=inverse,
            region_id=region_id
        )

    def init_regions(self):
        alphas = np.linspace(self.alpha_start, self.alpha_end, self.n_pts)
        dpts_arr = self.r * np.transpose([np.cos(alphas), np.sin(alphas)])
        dpts_arr = [DPoint(x, y) for x, y in dpts_arr]
        dpts_arr.append(DPoint(0, 0))

        disk_arc_polygon = SimplePolygon(DSimplePolygon(dpts_arr))
        self.metal_regions[self.region_id].insert(disk_arc_polygon)


class Ring(ElementBase):
    def __init__(self, origin, outer_r, thickness, n_pts=50, trans_in=None,
                 inverse=False, region_id="default"):
        """

        Parameters
        ----------
        origin : DPoint
            the center of the ring
        outer_r : float
            outer radius of the ring
        thickness : float
            thickness of the ring
        n_pts : int
            number of points comprising the circumference of the ring (50 by default)
        trans_in : Union[DCplxTrans, CplxTrans, DTrans, Trans]
            Represent object's local coordinate system transformation
            see more at
            https://www.klayout.de/doc-qt4/code/class_DCplxTrans.html
        inverse : bool
            Whether to interchange `empty_polygon` and
            `metal_polygon` before placing
        region_id : str
            region identifier to object if object is part of a compound
            object and/or resides in objects tree.
        """
        self.r = outer_r
        self.t = thickness
        self.n_pts = n_pts
        super().__init__(
            origin=origin, trans_in=trans_in, inverse=inverse,
            region_id=region_id
        )

    def init_regions(self):
        Rout = self.r
        Rin = self.r - self.t

        alphas = np.linspace(0, 2 * np.pi, self.n_pts + 1)
        dpts_arr_Rout = Rout * np.transpose([np.cos(alphas), np.sin(
            alphas)])
        dpts_arr_Rout = [DPoint(x, y) for x, y in dpts_arr_Rout]

        dpts_arr_Rin = Rout * np.transpose([np.cos(alphas), np.sin(
            alphas)])
        dpts_arr_Rin = [DPoint(x, y) for x, y in dpts_arr_Rin]

        ring_dpoly = DPolygon(dpts_arr_Rout)
        ring_dpoly.insert_hole(dpts_arr_Rin)
        ring_poly = Polygon(ring_dpoly)

        self.metal_region.insert(ring_poly)


class IsoTrapezoid(ElementBase):
    """@brief: class represents an isosceles trapezoid
        @params:  DPoint origin - position of the left bottom node
                        float b - b of the trapezoid
                        float bottom - length of the bottom side
                        float top - length of the top side
                        Trans trans_in - initial transformation (None by default)
                        bool inverse - if True then the ring is subtracted from width layer (False by default)
    """

    def __init__(self, origin, height, bottom, top, trans_in=None,
                 inverse=False):
        self.h = height
        self.b = bottom
        self.t = top
        super().__init__(origin, trans_in, inverse)

    def init_regions(self):
        origin = DPoint(0, 0)
        dx = (self.b - self.t) / 2
        p1 = origin + DPoint(dx, self.h)
        p2 = p1 + DPoint(self.t, 0)
        p3 = p2 + DPoint(dx, -self.h)
        pts_arr = [origin, p1, p2, p3]
        if self.inverse:
            self.empty_region.insert(
                SimplePolygon(DSimplePolygon(pts_arr)))
        else:
            self.metal_region.insert(
                SimplePolygon(DSimplePolygon(pts_arr)))


class Cross2(ElementBase):
    """@brief: class represents width cross
        @params:  DPoint origin - center of the cross
                        float thickness - thickness of the line
                        float length - size of the cross
                        Trans trans_in - initial transformation (None by default)
                        bool inverse - if True then the ring is subtracted from width layer (False by default)
    """

    def __init__(self, origin, thickness, length, trans_in=None,
                 inverse=False):
        self.l = length
        self.t = thickness
        super().__init__(origin, trans_in, inverse)

    def init_regions(self):
        origin = DPoint(0, 0)
        hor_box = DBox(origin - DPoint(self.l / 2, self.t / 2),
                       DPoint(self.l / 2, self.t / 2))
        vert_box = DBox(origin - DPoint(self.t / 2, self.l / 2),
                        DPoint(self.t / 2, self.l / 2))
        cross = (Region(hor_box) + Region(vert_box)).merge()
        if self.inverse:
            self.empty_region.insert(cross)
        else:
            self.metal_region.insert(cross)


''' NOT SO FREQUENTLY USED '''


class Kolbaska(ElementBase):
    """
    Warnings
    ---------
    Deprecated (Shrinks already existing DPath functionality).
    """

    def __init__(self, origin, stop, width, r, trans_in=None):
        self._width = width
        self._vec = stop - origin
        self._r = r
        super().__init__(origin, trans_in)
        self.start = self.connections[0]
        self.end = self.connections[1]

    def init_regions(self):
        ext_start = self._r
        ext_end = self._r
        shps = [self._width, ext_start, ext_end]
        kolb = DPath([DPoint(0, 0), DPoint(0, 0) + self._vec], *shps, True)
        self.metal_region.insert(kolb)
        self.connections.extend([DPoint(0, 0), DPoint(0, 0) + self._vec])
        self.angle_connections.extend([0, 0])

    def _refresh_named_connections(self):
        self.start = self.connections[0]
        self.end = self.connections[1]


class DPathCL(DPath, ElementBase):
    def __init__(self, pts: List[DPoint], width: float, bgn_ext: float = 0,
                 end_ext: float = 0, is_round: bool = False, bendings_r=0,
                 trans_in=None, inverse=None):
        """
        Constructor given the points of the path's spine, the width,
        the extensions and the round end flag.

        Parameters
        ----------
        pts : List[DPoint]
            The points forming the spine of the path
        width : float
            The width of the path
        bgn_ext : float
            The start extension of the path
        end_ext : float
            The end extension of the path
        is_round : bool
            If this flag is true, the path will get rounded ends
        bendings_r : float
            Creates bent DPath by calling `DPath.round(r=bendings_r)`
        """
        if bendings_r < width:
            Warning("inner roundings will be rendered incorrectly due to "
                    "`bendings_r` is lesser than width `width` of width `DPath`")
        self.pts: List[DPoint] = pts
        self.width = width
        self.bgn_ext = bgn_ext
        self.end_ext = end_ext
        self.round = is_round
        self.start: DPoint = DPoint(0, 0)
        self.end: DPoint = DPoint(0, 0)
        self.connections = self.pts
        # TODO: implement connection angles
        DPath.__init__(self, pts, width, bgn_ext, end_ext, is_round)
        self.polygon = self.round_corners(
            bendings_r, PROGRAM.ARC_PTS_N, 1
        ).polygon()
        ElementBase.__init__(self, DPoint(0, 0), trans_in=trans_in)

    def init_regions(self):
        self.connections = self.pts
        self.start = self.connections[0]
        self.end = self.connections[-1]
        self.metal_region.insert(self.polygon)

    def _refresh_named_connections(self):
        self.pts = self.connections
        self.start = self.connections[0]
        self.end = self.connections[-1]
