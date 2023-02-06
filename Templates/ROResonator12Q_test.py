import pya
from pya import Point, DPoint, Vector, DVector, DSimplePolygon, SimplePolygon, DPolygon, Polygon, Region, Box, DBox
from pya import Trans, DTrans, CplxTrans, DCplxTrans, ICplxTrans

from importlib import reload

import classLib
reload(classLib)
from classLib.chipDesign import ChipDesign
from classLib.shapes import Donut


class MyDesign(ChipDesign):
    def draw(self):
        origin = DPoint(0, 0)
        self.region_ph.insert(DBox(DPoint(-1e6, -1e6), DPoint(1e6, 1e6)))
        self.donut = Donut(
            origin, inner_r=200e3, outer_r=300e3,
            alpha_start=0, alpha_end=90,
            n_pts=100,
            inverse=True
        )
        self.donut.place(self.region_ph)




import numpy as np
### MAIN FUNCTION ###
if __name__ == "__main__":
    my_design = MyDesign("testScript")
    my_design.draw()
    my_design.show()



