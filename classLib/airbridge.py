import pya
from pya import Point, DPoint, DSimplePolygon, SimplePolygon, DPolygon, Polygon, Region, Box
from pya import Trans, DTrans, CplxTrans, DCplxTrans, ICplxTrans, CellInstArray

from classLib.baseClasses import *


class Airbridge(ElementBase):

    def __init__(self, origin, trans_in):
        super().__init__(origin.dup(), trans_in)

    def init_regions(self):
        bridge_length = 85e3
        bridge_width = 20e3
        patch_width = 40e3
        patch_length = 30e3

        bridge = Box(DPoint(-bridge_width / 2, -bridge_length / 2), DPoint(bridge_width / 2, bridge_length / 2))

        patch1 = Box(DPoint(-patch_width / 2, -bridge_length / 2 - patch_length / 3),
                     DPoint(patch_width / 2, -bridge_length / 2 + 2 * patch_length / 3))
        patch2 = Box(DPoint(-patch_width / 2, bridge_length / 2 + patch_length / 3),
                     DPoint(patch_width / 2, bridge_length / 2 - 2 * patch_length / 3))

        self.metal_regions["bridges"] = Region(bridge)
        self.metal_regions["bridge_patches"] = Region(patch1) + Region(patch2)
        self.empty_regions["bridges"] = Region()
        self.empty_regions["bridge_patches"] = Region()
