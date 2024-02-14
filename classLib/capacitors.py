# import built-ins
# ---
# import good 3rd party
# ---
# import project specific 3rd party
import pya
from math import sqrt, cos, sin, tan, atan2, pi, copysign
from pya import Point, DPoint, Vector, DVector, DSimplePolygon, SimplePolygon, DPolygon, Polygon, Region
from pya import Trans, DTrans, CplxTrans, DCplxTrans, ICplxTrans

# import project lib
from classLib.baseClasses import ComplexBase
from classLib.coplanars import *  # TODO: get rid of * statements
from classLib.shapes import *


class CWave2CPW(ElementBase):
    '''
    Draws width semi-circle coupler from coplanar waveguide to jelly capacitance plates.
    '''

    def __init__(self, c_wave_cap, params, n_pts=50, trans_in=None):
        self.c_wave_ref = c_wave_cap
        if isinstance(params, dict):
            self.Z1 = params['Z1']
            self.d_alpha1 = params['d_alpha1']
            self.width1 = params['width1']
            self.gap1 = params['gap1']
            self.Z2 = params['Z2']
            self.d_alpha2 = params['d_alpha2']
            self.width2 = params['width2']
            self.gap2 = params['gap2']
        else:
            # not recommended
            self.Z1 = params[0]
            self.d_alpha1 = params[1]
            self.width1 = params[2]
            self.gap1 = params[3]
            self.Z2 = params[4]
            self.d_alpha2 = params[5]
            self.width2 = params[6]
            self.gap2 = params[7]
        self.n_pts = n_pts
        super().__init__(self.c_wave_ref.origin.dup(), trans_in)

    def _get_solid_arc(self, center, R, width, alpha_start, alpha_end, n_inner, n_outer):
        pts = []

        d_alpha_inner = (alpha_end - alpha_start) / (n_inner - 1)
        d_alpha_outer = -(alpha_end - alpha_start) / (n_outer - 1)

        for i in range(0, n_inner):
            alpha = alpha_start + d_alpha_inner * i
            pts.append(center + DPoint(cos(alpha), sin(alpha)) * (R - width / 2))
        for i in range(0, n_outer):
            alpha = alpha_end + d_alpha_outer * i
            pts.append(center + DPoint(cos(alpha), sin(alpha)) * (R + width / 2))

        return DSimplePolygon(pts)

    def init_regions(self):
        origin = DPoint(0, 0)
        arc_1_solid = self._get_solid_arc(origin, self.c_wave_ref.in_circle.r + self.gap1 + self.width1 / 2,
                                          self.width1,
                                          pi / 2 - self.d_alpha1 / 2, pi / 2 + self.d_alpha1 / 2,
                                          self.n_pts, self.n_pts)
        arc_1_empty = self._get_solid_arc(origin, self.c_wave_ref.in_circle.r + self.gap1 / 2, self.gap1,
                                          pi / 2 - self.d_alpha1 / 2, pi / 2 + self.d_alpha1 / 2,
                                          self.n_pts, self.n_pts)

        arc_2_solid = self._get_solid_arc(origin, self.c_wave_ref.in_circle.r + self.gap2 + self.width2 / 2,
                                          self.width2,
                                          3 / 2 * pi - self.d_alpha2 / 2, 3 / 2 * pi + self.d_alpha2 / 2,
                                          self.n_pts, self.n_pts)
        arc_2_empty = self._get_solid_arc(origin, self.c_wave_ref.in_circle.r + self.gap2 / 2, self.gap2,
                                          3 / 2 * pi - self.d_alpha2 / 2, 3 / 2 * pi + self.d_alpha2 / 2,
                                          self.n_pts, self.n_pts)
        self.metal_region.insert(arc_1_solid)
        self.metal_region.insert(arc_2_solid)
        self.empty_region.insert(arc_1_empty)
        self.empty_region.insert(arc_2_empty)


class CWave(ComplexBase):
    '''
    Draws width condensator from width circle cutted into 2 pieces.
    '''

    def __init__(self, center, r_out, dr, n_segments, s, alpha, r_curve, delta=40e3, n_pts=50, solid=True,
                 trans_in=None):
        '''
        Parameters:
        center: DPoint
            A center of circle.
        r_out: float
            The outer outer_r of width circle (along to the edge of the ground).
        dr: float
            The gap in between the common ground and the circle perimeter.
        n_segments: int
            The number of segments (without turns) composed into width condensator gap.
        s: float
            The half-width of the gap.
        alpha: rad
            The angle of single arc of slice.
        r_curve: float
            The outer_r of single arc.
        delta: float
            length of the horizontal lines on the ends of the cut
        n_pts: int
            The number of points on the perimeter of the circle.
        inverse: ???
        trans_in: Bool
            Initial transformation
        '''
        self.r_out = r_out
        self.dr = dr
        self.r_in = self.r_out - self.dr
        self.n_segments = n_segments
        self.s = s
        self.alpha = alpha
        self.r_curve = r_curve
        self.n_pts = n_pts
        self.delta = delta
        self.L_full = 2 * (self.r_out - self.dr) - 2 * self.delta
        # calculating parameters of the CPWDPath #
        L_full = self.L_full
        alpha = self.alpha
        if abs(self.alpha) != pi:
            self.L1 = L_full / ((self.n_segments + 1) * cos(alpha))
            self.L0 = self.L1 / 2
        else:
            raise ValueError("180 degrees turns in CWave are not supported.")

        if (self.L0 < 0):
            print("CPWDPath: impossible parameters combination")

        super().__init__(center, trans_in)

    def init_primitives(self):
        origin = DPoint(0, 0)
        # erased line params
        Z = CPWParameters(0, self.s / 2)
        shapes = ''
        angles = []
        lengths = []
        # placing circle r_out with dr clearance from ground polygon
        self.empt_circle = Disk(origin, self.r_out, n_pts=self.n_pts, inverse=True)
        self.in_circle = Disk(origin, self.r_out - self.dr, n_pts=self.n_pts, inverse=False)
        self.empt_circle.empty_region -= self.in_circle.metal_region
        self.primitives["empt_circle"] = self.empt_circle
        self.primitives["in_circle"] = self.in_circle

        self.RL_start = origin + DPoint(-self.in_circle.r, 0)
        shapes += 'LRL'
        angles.extend([self.alpha])
        lengths.extend([self.delta, self.L0])
        # rl_path_start = CPWDPath(self.RL_start, "LRLR", Z, self.r_curve, [self.delta, self.L0], [self.alpha,-self.alpha] )
        # self.primitives["rl_start"] = rl_path_start

        # intermidiate RLRs
        for i in range(self.n_segments):
            if (i % 2 == 1):
                m_x = -1
            else:
                m_x = 1
            shapes += 'RL'
            angles.extend([-2 * m_x * self.alpha])
            lengths.extend([self.L1])
            # prev_path = list(self.primitives.values())[-1]
            # rl_path_p = CPWDPath( prev_path.end, "RLR", Z, self.r_curve, [self.L1], [-m_x*self.alpha,m_x*self.alpha])
            # self.primitives["rl_path_" + str(k)] = rl_path_p
            # ending RLR
        if (self.n_segments % 2 == 1):
            m_x = 1
        else:
            m_x = -1
        shapes += 'RLRL'
        angles.extend([2 * m_x * self.alpha, -m_x * self.alpha])
        lengths.extend([self.L0, self.delta])
        cut = CPWRLPath(self.RL_start, shapes, Z, self.r_curve, lengths, angles)
        # prev_path = list(self.primitives.values())[-1]
        # rl_path_end = CPWDPath( prev_path.end, "RLRL", Z, self.r_curve, [self.L0, self.delta], [m_x*self.alpha,-m_x*self.alpha])
        self.primitives["cut"] = cut


class XmonCross(ComplexBase):
    def __init__(self, origin,
                 sideX_length, sideX_width, sideX_gnd_gap,
                 sideX_face_gnd_gap=None,
                 sideY_length=None, sideY_width=None, sideY_gnd_gap=None,
                 sideY_face_gnd_gap=None,
                 trans_in=None):
        """
        Draws cross for xmon qubit with width lot of customization parameters.

        Parameters
        ----------
        origin : DPoint
            center of the cross
        sideX_length : float
            length of the cross extensions along x-axis
        sideX_width : float
            width of the cross extensions along x-axis
        sideX_gnd_gap : float
            ground gap between longs sides of extensions along x-axis and ground (gap along y-axis)
        sideX_face_gnd_gap : float
            ground gap between face end of extensions along x-axis and ground (gap along x-axis)
            default - `sideX_gnd_gap`
        sideY_length : float
            length of the cross extensions along y-axis
            default - `sideX_length`
        sideY_width : float
            width of the cross extensions along y-axis
            default - `sideX_width`
        sideY_gnd_gap : float
            ground gap between face end of extensions along y-axis and ground (gap along x-axis)
            default - `sideX_gnd_gap`
        sideY_face_gnd_gap : float
            ground gap between face end of extensions along y-axis and ground (gap along y-axis)
            default - `sideX_face_gnd_gap`
        trans_in : DCplxTrans
            transformation of the cross

        Notes
        -----------
        if `sideX_face_gnd_gap` and `sideY_face_gnd_gap` both are ommited, then the latter
        will be equal `sideX_face_gnd_gap` default value that is `sideX_gnd_gap`
        """
        self.sideX_length = sideX_length
        self.sideY_length = None
        if sideY_length is None:
            self.sideY_length = self.sideX_length
        else:
            self.sideY_length = sideY_length

        self.sideX_width = sideX_width
        self.sideY_width = None
        if sideY_width is None:
            self.sideY_width = self.sideX_width
        else:
            self.sideY_width = sideY_width

        self.sideX_gnd_gap = sideX_gnd_gap
        self.sideY_gnd_gap = None
        if sideY_gnd_gap is None:
            self.sideY_gnd_gap = self.sideX_gnd_gap
        else:
            self.sideY_gnd_gap = sideY_gnd_gap

        if sideX_face_gnd_gap is None:
            self.sideX_face_gnd_gap = self.sideX_gnd_gap
        else:
            self.sideX_face_gnd_gap = sideX_face_gnd_gap

        if sideY_face_gnd_gap is None:
            self.sideY_face_gnd_gap = self.sideY_gnd_gap
        else:
            self.sideY_face_gnd_gap = sideY_face_gnd_gap

        # for saving
        self.center = origin
        super().__init__(self.center.dup(), trans_in)
        self._geometry_parameters[
            "sideX_length, um"] = self.sideX_length / 1e3
        self._geometry_parameters[
            "sideY_length, um"] = self.sideY_length / 1e3
        self._geometry_parameters[
            "sideX_width, um"] = self.sideX_width / 1e3
        self._geometry_parameters[
            "sideY_width, um"] = self.sideY_width / 1e3
        self._geometry_parameters[
            "sideX_gnd_gap, um"] = self.sideX_gnd_gap / 1e3
        self._geometry_parameters[
            "sideY_gnd_gap, um"] = self.sideY_gnd_gap / 1e3
        self._geometry_parameters[
            "sideX_face_gnd_gap, um"] = self.sideX_face_gnd_gap / 1e3
        self._geometry_parameters[
            "sideY_face_gnd_gap, um"] = self.sideY_face_gnd_gap / 1e3
        self.center = self.connections[0]

    def init_primitives(self):
        origin = DPoint(0, 0)

        # draw central square
        from classLib.shapes import Rectangle
        lb_corner = DPoint(-self.sideY_width / 2, -self.sideX_width / 2)
        center_square = Rectangle(lb_corner, self.sideY_width, self.sideX_width)
        self.primitives["center_square"] = center_square

        """ left part of Xmon cross """
        p1 = origin + DPoint(-self.sideY_width / 2, 0)
        p2 = p1 + DPoint(-self.sideX_length, 0)
        self.cpw_l = CPW(self.sideX_width, self.sideX_gnd_gap, p1, p2)
        self.primitives["cpw_l"] = self.cpw_l
        p3 = p2 + DPoint(-self.sideX_face_gnd_gap, 0)
        self.cpw_lempt = CPW(0, self.cpw_l.b / 2, p2, p3)
        self.primitives["cpw_lempt"] = self.cpw_lempt

        """ right part of Xmon cross """
        p1 = origin + DPoint(self.sideY_width / 2, 0)
        p2 = p1 + DPoint(self.sideX_length, 0)
        self.cpw_r = CPW(self.sideX_width, self.sideX_gnd_gap, p1, p2)
        self.primitives["cpw_r"] = self.cpw_r
        p3 = p2 + DPoint(self.sideX_face_gnd_gap, 0)
        self.cpw_rempt = CPW(0, self.cpw_r.b / 2, p2, p3)
        self.primitives["cpw_rempt"] = self.cpw_rempt

        """ top part of Xmon cross """
        p1 = origin + DPoint(0, self.sideX_width / 2)
        p2 = p1 + DPoint(0, self.sideY_length)
        self.cpw_t = CPW(self.sideY_width, self.sideY_gnd_gap, p1, p2)
        self.primitives["cpw_t"] = self.cpw_t
        p3 = p2 + DPoint(0, self.sideY_face_gnd_gap)
        self.cpw_tempt = CPW(0, self.cpw_t.b / 2, p2, p3)
        self.primitives["cpw_tempt"] = self.cpw_tempt

        """ bottom part of Xmon cross """
        p1 = origin + DPoint(0, -self.sideX_width / 2)
        p2 = p1 + DPoint(0, -self.sideY_length)
        self.cpw_b = CPW(self.sideY_width, self.sideY_gnd_gap, p1, p2)
        self.primitives["cpw_b"] = self.cpw_b
        p3 = p2 + DPoint(0, -self.sideY_face_gnd_gap)
        self.cpw_bempt = CPW(0, self.cpw_b.b / 2, p2, p3)
        self.primitives["cpw_bempt"] = self.cpw_bempt

        self.connections = [origin]

    def _refresh_named_connections(self):
        self.center = self.connections[0]


class TmonT(ComplexBase):
    def __init__(self, origin,
                 sideX_length, sideX_width, sideX_gnd_gap,
                 sideX_face_gnd_gap=None,
                 sideY_length=None, sideY_width=None, sideY_gnd_gap=None,
                 sideY_face_gnd_gap=None,
                 trans_in=None):
        """
        Draws T-shape for Tmon qubit (transmon with T-shaped shunting
        capacitor) with a lot of customization parameters.

        Parameters
        ----------
        origin : DPoint
            center of the cross
        sideX_length : float
            length of the cross extensions along x-axis
        sideX_width : float
            width of the cross extensions along x-axis
        sideX_gnd_gap : float
            ground gap between longs sides of extensions along x-axis and ground (gap along y-axis)
        sideX_face_gnd_gap : float
            ground gap between face end of extensions along x-axis and ground (gap along x-axis)
            default - `sideX_gnd_gap`
        sideY_length : float
            length of the cross extensions along y-axis
            default - `sideX_length`
        sideY_width : float
            width of the cross extensions along y-axis
            default - `sideX_width`
        sideY_gnd_gap : float
            ground gap between face end of extensions along y-axis and ground (gap along x-axis)
            default - `sideX_gnd_gap`
        sideY_face_gnd_gap : float
            ground gap between face end of extensions along y-axis and ground (gap along y-axis)
            default - `sideX_face_gnd_gap`
        trans_in : DCplxTrans
            transformation of the cross

        Notes
        -----------
        if `sideX_face_gnd_gap` and `sideY_face_gnd_gap` both are ommited, then the latter
        will be equal `sideX_face_gnd_gap` default value that is `sideX_gnd_gap`
        """
        self.sideX_length = sideX_length
        self.sideY_length = None
        if sideY_length is None:
            self.sideY_length = self.sideX_length
        else:
            self.sideY_length = sideY_length

        self.sideX_width = sideX_width
        self.sideY_width = None
        if sideY_width is None:
            self.sideY_width = self.sideX_width
        else:
            self.sideY_width = sideY_width

        self.sideX_gnd_gap = sideX_gnd_gap
        self.sideY_gnd_gap = None
        if sideY_gnd_gap is None:
            self.sideY_gnd_gap = self.sideX_gnd_gap
        else:
            self.sideY_gnd_gap = sideY_gnd_gap

        if sideX_face_gnd_gap is None:
            self.sideX_face_gnd_gap = self.sideX_gnd_gap
        else:
            self.sideX_face_gnd_gap = sideX_face_gnd_gap

        if sideY_face_gnd_gap is None:
            self.sideY_face_gnd_gap = self.sideY_gnd_gap
        else:
            self.sideY_face_gnd_gap = sideY_face_gnd_gap

        # for saving
        self.center = origin
        super().__init__(origin.dup(), trans_in)
        self._geometry_parameters[
            "sideX_length, um"] = self.sideX_length / 1e3
        self._geometry_parameters[
            "sideY_length, um"] = self.sideY_length / 1e3
        self._geometry_parameters[
            "sideX_width, um"] = self.sideX_width / 1e3
        self._geometry_parameters[
            "sideY_width, um"] = self.sideY_width / 1e3
        self._geometry_parameters[
            "sideX_gnd_gap, um"] = self.sideX_gnd_gap / 1e3
        self._geometry_parameters[
            "sideY_gnd_gap, um"] = self.sideY_gnd_gap / 1e3
        self._geometry_parameters[
            "sideX_face_gnd_gap, um"] = self.sideX_face_gnd_gap / 1e3
        self._geometry_parameters[
            "sideY_face_gnd_gap, um"] = self.sideY_face_gnd_gap / 1e3
        self.center = self.connections[0]

    def init_primitives(self):
        origin = DPoint(0, 0)

        # draw central square
        from classLib.shapes import Rectangle
        lb_corner = DPoint(-self.sideX_width / 2, -self.sideY_width / 2)
        center_square = Rectangle(lb_corner, self.sideX_width,
                                  self.sideY_width)
        self.primitives["center_square"] = center_square

        """ left part of Xmon cross """
        p1 = origin + DPoint(-self.sideY_width / 2, 0)
        p2 = p1 + DPoint(-self.sideX_length, 0)
        self.cpw_l = CPW(self.sideX_width, self.sideX_gnd_gap, p1, p2)
        self.primitives["cpw_l"] = self.cpw_l
        p3 = p2 + DPoint(-self.sideX_face_gnd_gap, 0)
        self.cpw_lempt = CPW(0, self.cpw_l.b / 2, p2, p3)
        self.primitives["cpw_lempt"] = self.cpw_lempt

        """ right part of Tmon T """
        p1 = origin + DPoint(self.sideY_width / 2, 0)
        p2 = p1 + DPoint(self.sideX_length, 0)
        self.cpw_r = CPW(self.sideX_width, self.sideX_gnd_gap, p1, p2)
        self.primitives["cpw_r"] = self.cpw_r
        p3 = p2 + DPoint(self.sideX_face_gnd_gap, 0)
        self.cpw_rempt = CPW(0, self.cpw_r.b / 2, p2, p3)
        self.primitives["cpw_rempt"] = self.cpw_rempt

        """ top part of Tmon T """
        p1 = origin + DPoint(0, self.sideX_width / 2)
        p2 = p1 + DPoint(0, self.sideY_length)
        self.cpw_t = CPW(self.sideY_width, self.sideY_gnd_gap, p1, p2)
        self.primitives["cpw_t"] = self.cpw_t
        p3 = p2 + DPoint(0, self.sideY_face_gnd_gap)
        self.cpw_tempt = CPW(0, self.cpw_t.b / 2, p2, p3)
        self.primitives["cpw_tempt"] = self.cpw_tempt

        """ bottom part of Tmon T (for backward compatibility) """
        p1 = origin + DPoint(0, -self.sideX_width / 2)
        p2 = p1
        self.cpw_b = CPW(self.sideY_width, self.sideY_gnd_gap, p1, p2)
        self.primitives["cpw_b"] = self.cpw_b
        p3 = p2 + DPoint(0, -self.sideX_gnd_gap)
        self.cpw_bempt = CPW(0, self.cpw_b.b / 2, p2, p3)
        self.primitives["cpw_bempt"] = self.cpw_bempt

        self.connections = [origin]

    def _refresh_named_connections(self):
        self.center = self.connections[0]
