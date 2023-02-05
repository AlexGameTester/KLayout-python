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
    def __init__(
        self, origin, width, height, trans_in=None,
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
                SimplePolygon(DSimplePolygon(pts_arr))
            )
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
    def __init__(
        self,
        origin,
        inner_square_a,
        outer_square_a,
        trans_in=None,
        inverse=False,
        region_id="default"
    ):
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


# DEPRECATED: USE `classlib.shapes.Donut` instead
class Disk(ElementBase):
    def __init__(
        self, center, r, n_pts=50, inverse=False,
        offset_angle=0, trans_in=None, region_id="default"
    ):
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
        super().__init__(
            origin=center, trans_in=trans_in, inverse=inverse,
            region_id=region_id
        )

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


# DEPRECATED: USE `classlib.shapes.Donut` instead
class DiskSector(ElementBase):
    def __init__(
        self, center, r, alpha_start=0, alpha_end=np.pi,
        n_pts=50, trans_in=None, inverse=False,
        region_id="default"
    ):
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


# DEPRECATED: USE `classlib.shapes.Donut` instead
class RingSector(ElementBase):
    def __init__(
        self, origin, outer_r, thickness, n_pts=50, trans_in=None,
        inverse=False, region_id="default"
    ):
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
        dpts_arr_Rout = Rout * np.transpose(
            [np.cos(alphas), np.sin(
                alphas
            )]
        )
        dpts_arr_Rout = [DPoint(x, y) for x, y in dpts_arr_Rout]

        dpts_arr_Rin = Rout * np.transpose(
            [np.cos(alphas), np.sin(
                alphas
            )]
        )
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
                        bool inverse - if True then the ring is subtracted from width layer (
                        False by default)
    """

    def __init__(
        self, origin, height, bottom, top, trans_in=None,
        inverse=False
    ):
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
                SimplePolygon(DSimplePolygon(pts_arr))
            )
        else:
            self.metal_region.insert(
                SimplePolygon(DSimplePolygon(pts_arr))
            )


class Cross2(ElementBase):
    """@brief: class represents width cross
        @params:  DPoint origin - center of the cross
                        float width - width of the line
                        float length - size of the cross
                        Trans trans_in - initial transformation (None by default)
                        bool inverse - if True then the ring is subtracted from width layer (
                        False by default)
    """

    def __init__(
        self, origin, thickness, length, trans_in=None,
        inverse=False
    ):
        self.l = length
        self.t = thickness
        super().__init__(origin, trans_in, inverse)

    def init_regions(self):
        origin = DPoint(0, 0)
        hor_box = DBox(
            origin - DPoint(self.l / 2, self.t / 2),
            DPoint(self.l / 2, self.t / 2)
        )
        vert_box = DBox(
            origin - DPoint(self.t / 2, self.l / 2),
            DPoint(self.t / 2, self.l / 2)
        )
        cross = (Region(hor_box) + Region(vert_box)).merge()
        if self.inverse:
            self.empty_region.insert(cross)
        else:
            self.metal_region.insert(cross)


class Donut(ElementBase):
    def __init__(
        self, origin, inner_r=None, outer_r=None, ring_width=None,
        alpha_start=0, alpha_end=360, n_pts=50,
        trans_in=None,
        inverse=False, region_id="default"
    ):
        """
        Class that draw circles/rings and their angled sectors.
        Angles supplied is counted counter-clockwise from Ox axis in radii.

        Inspired by the structure of the same name from the sonnet 14.52 geometry design toolbox.

        Parameters
        ----------
        origin : DPoint
            the center of the ring
        inner_r: float
            inner radius of the ring
        outer_r : float
            outer radius of the ring
        ring_width : float
            thickness of the ring
            If all 3 parameters from the set {`inner_r`, `outer_r`, `ring_width`} are supplied
            than `ring_width` will be ignored.
        alpha_start : float
            start angle in degree. For ring sector. Count starts from x axis in counter-clockwise
            direction.
            Default = 0
        alpha_end : float
            end angle in degree. For ring sector
            Default = 360
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

        def parseRinRoutRwidth(inner_r, outer_r, ring_width):
            # parsing {`inner_r`, `outer_r`, `ring_width`} set to be self-consistent
            # If 2 out of 3 are supplied - third is calculated automatically
            # If <2 or all 3 are supplied - ValueError will be thrown
            counter = 0
            for entry in [outer_r, inner_r, ring_width]:
                if entry is not None:
                    counter += 1
            if (counter == 3) or (counter < 2):
                raise ValueError(
                    "Only 2 out of 3 keyword arguments with following names has to be \n"
                    "suuplied: {`inner_r`, `outer_r`, `ring_width`}\n"
                    f"You supplied {counter} arguments"
                )

            if inner_r is None:
                inner_r = outer_r - ring_width
            if outer_r is None:
                outer_r = inner_r + ring_width
            if ring_width is None:
                ring_width = outer_r - inner_r

            return inner_r, outer_r, ring_width

        inner_r, outer_r, ring_width = parseRinRoutRwidth(
            inner_r=inner_r, outer_r=outer_r, ring_width=ring_width
        )
        self.inner_r = inner_r
        self.outer_r = outer_r
        self.ring_width = ring_width
        self.n_pts = n_pts
        self.alpha_start = alpha_start
        self.alpha_end = alpha_end
        self.outer_arc_center: DPoint = None
        self.inner_arc_center: DPoint = None
        super().__init__(
            origin=origin, trans_in=trans_in, inverse=inverse,
            region_id=region_id
        )

    def init_regions(self):
        # name aliases for short adressing
        Rout = self.outer_r
        Rin = self.inner_r

        ''' Check if values are self-consistent '''

        # unpredicted behaviour in supplied parameters point
        if Rin >= Rout:
            print(
                f"{self.__class__.__name__} error:\n"
                f"{Rin} = inner_r >= outer_r = {Rout}\n"
                f"unpredicted behaviour for this supplied parameters region"
            )
        if self.alpha_start > self.alpha_end:
            print(
                f"{self.__class__.__name__} error:\n"
                f"{self.alpha_start} = alpha_start >= alpha_end = {self.alpha_end}\n"
                f"unpredicted behaviour for this supplied parameters region"
            )
        if abs(self.alpha_end - self.alpha_start) >= 360:
            print(
                f"{self.__class__.__name__} error:\n"
                f"angle difference surpass 2pi: alpha difference\n"
                f" = {abs(self.alpha_end-self.alpha_start)}\n"
                f"underfined behaviour for this parameter point"
            )

        # empty region returned
        if Rout == 0:
            print(
                f"{self.__class__.__name__} error:\n"
                f"outer radius = {Rout} <= 0\n"
                f"empty region will be constructed"
            )
        elif self.alpha_start == self.alpha_end:
            print(
                f"{self.__class__.__name__} error:\n"
                f"alpha_start =  alpha_end = {self.alpha_end}\n"
                f"empty region will be constructed"
            )
        else:  # proper region returned
            ''' Drawing multiple scenarious '''
            alphas = 2*np.pi/360*np.linspace(self.alpha_start, self.alpha_end, self.n_pts)
            # always draw hull (exterior part)
            dpts_arr_Rout = Rout * np.transpose(
                [np.cos(alphas), np.sin(alphas)]
            )
            dpts_arr_Rout = [DPoint(x, y) for x, y in dpts_arr_Rout]
            donut_dpolygon = DPolygon()
            if abs((self.alpha_start - self.alpha_end) % (360)) <= 1e-8:  # complete donut
                # exclude endpoint for closed hull contour since first and last points coincide
                donut_dpolygon.hull = dpts_arr_Rout[:-1]  # DPolygon allow holes inside
                if Rin != 0:  # insert disk hole
                    dpts_arr_Rin = Rin * np.transpose([np.cos(alphas), np.sin(alphas)])
                    # exclude endpoint for closed inner contour since first and last points simply
                    # coincide
                    dpts_arr_Rin = [DPoint(x, y) for x, y in dpts_arr_Rin[:-1]]
                    donut_dpolygon.insert_hole(dpts_arr_Rin)
                else:  # complete disk
                    # no holes needed
                    pass
            else:  # angle sector of a donut
                # no endpoints excluded
                angle_sector_pts = dpts_arr_Rout
                if Rin == 0:  # disk sector
                    # add center to the hull to complete contour
                    angle_sector_pts += [DPoint(0, 0)]
                else:
                    # donut angle sector
                    dpts_arr_Rin = Rin * np.transpose([np.cos(alphas), np.sin(alphas)])
                    # clockwise direction due to `dpts_arr_Rin[::-1]`
                    dpts_arr_Rin = [DPoint(x, y) for x, y in dpts_arr_Rin[::-1]]
                    angle_sector_pts += dpts_arr_Rin
                donut_dpolygon = DPolygon(angle_sector_pts)

            if self.inverse:  # TODO: implement inverse in `ElementBase` and `ComplexBase` probably.
                self.empty_region.insert(donut_dpolygon)
            else:
                self.metal_region.insert(donut_dpolygon)
            alpha_center_rad = (self.alpha_start + self.alpha_end)/2*(2*np.pi/360)
            outer_arc_center = Rout*DVector(np.cos(alpha_center_rad), np.sin(alpha_center_rad))
            inner_arc_center = Rin*DVector(np.cos(alpha_center_rad), np.sin(alpha_center_rad))
            self.connections = [DPoint(0, 0), outer_arc_center, inner_arc_center]

    def _refresh_named_connections(self):
        self.center = self.connections[0]
        self.outer_arc_center = self.connections[1]
        self.inner_arc_center = self.connections[2]


''' NOT SO FREQUENTLY USED '''


# DEPRECETED: Use classlib.coplanars.DPathCPW instead
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


# DEPRECATED: Use `classlib.shapes.DPathCPW` instead
class DPathCL(DPath, ElementBase):
    def __init__(
        self, pts: List[DPoint], width: float, bgn_ext: float = 0,
        end_ext: float = 0, is_round: bool = False, bendings_r=0,
        trans_in=None, inverse=None
    ):
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
            Warning(
                "inner roundings will be rendered incorrectly due to "
                "`bendings_r` is lesser than width `width` of width `DPath`"
            )
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
