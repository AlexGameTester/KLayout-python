import pya
from math import sqrt, cos, sin, atan2, pi, copysign
from pya import Point, DPoint, DSimplePolygon, SimplePolygon, DPolygon, Polygon, Region
from pya import Trans, DTrans, CplxTrans, DCplxTrans, ICplxTrans

from classLib.baseClasses import *  # TODO: get rid of `*` - assign to @vdrhtc


class Claw(ElementBase):
    # TODO: add usage example, docs and schematic - assign to @vdrhtc
    def __init__(
        self, origin, resonator_cpw_params, claw_aperture,
        w_thin_ground=5e3, w_claw_pad=20e3, l_claw_pad=60e3, w_claw=None,
        trans_in=None
    ):
        self._resonator_cpw_params = resonator_cpw_params
        self._claw_aperture = claw_aperture
        self._w_thin_ground = w_thin_ground
        self._w_claw_pad = w_claw_pad
        self._l_claw_pad = l_claw_pad
        self._w_claw = resonator_cpw_params.width if w_claw is None else w_claw

        super().__init__(origin.dup(), trans_in)

        self.start = self.connections[0]
        self.end = self.connections[1]

    @staticmethod
    def get_phase_shift(frequency):
        return 0.014 + 1.436 * frequency;

    def init_regions(self):
        w_res, g_res = self._resonator_cpw_params.width, self._resonator_cpw_params.gap
        w_protect_pad = self._w_claw_pad + g_res * 2
        protect_aperture = self._claw_aperture + 2 * self._w_thin_ground
        protect_width = protect_aperture + 2 * w_protect_pad
        l_claw_pad = self._l_claw_pad
        w_claw = self._w_claw

        protect_points = [DPoint(0, -g_res),
                          DPoint(-protect_width / 2, -g_res),
                          DPoint(-protect_width / 2, w_claw + l_claw_pad + g_res),
                          DPoint(-protect_width / 2 + w_protect_pad, w_claw + l_claw_pad + g_res),
                          DPoint(-protect_width / 2 + w_protect_pad, w_claw + g_res),
                          DPoint(protect_width / 2 - w_protect_pad, w_claw + g_res),
                          DPoint(protect_width / 2 - w_protect_pad, w_claw + l_claw_pad + g_res),
                          DPoint(protect_width / 2, w_claw + l_claw_pad + g_res),
                          DPoint(protect_width / 2, -g_res)]

        protect_region = Region(DSimplePolygon(protect_points))

        metal_points = [DPoint(0, 0),
                        DPoint(-protect_width / 2 + g_res, 0),
                        DPoint(-protect_width / 2 + g_res, w_claw + l_claw_pad),
                        DPoint(-protect_width / 2 + w_protect_pad - g_res, w_claw + l_claw_pad),
                        DPoint(-protect_width / 2 + w_protect_pad - g_res, w_claw),
                        DPoint(protect_width / 2 - w_protect_pad + g_res, w_claw),
                        DPoint(protect_width / 2 - w_protect_pad + g_res, w_claw + l_claw_pad),
                        DPoint(protect_width / 2 - g_res, w_claw + l_claw_pad),
                        DPoint(protect_width / 2 - g_res, 0)]

        metal_poly = DSimplePolygon(metal_points)
        self.metal_region.insert(metal_poly)

        empty_region = protect_region - self.metal_region

        self.empty_region.insert(empty_region)

        self.connections = [DPoint(0, 0), DPoint(0, w_claw + g_res + self._w_thin_ground)]
