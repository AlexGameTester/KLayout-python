import numpy as np

import pya
from pya import Point, DPoint, Vector, DVector, DSimplePolygon, SimplePolygon, DPolygon, Polygon, \
    Region, Box, DBox
from pya import Trans, DTrans, CplxTrans, DCplxTrans, ICplxTrans

from importlib import reload

import classLib

reload(classLib)
from classLib.chipDesign import ChipDesign
from classLib.shapes import Donut


class MyDesign(ChipDesign):
    def draw(self):
        origin = DPoint(0, 0)
        donut = Donut(
            origin, inner_r=200e3, outer_r=300e3,
            alpha_start=0, alpha_end=10*np.pi/2,
            n_pts=100,
            inverse=False
        )
        donut.place(self.region_ph)


### MAIN FUNCTION ###
if __name__ == "__main__":
    my_design = MyDesign("testScript")
    my_design.draw()
    my_design.show()
