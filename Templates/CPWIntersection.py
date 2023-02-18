from importlib import reload
import numpy as np

import pya
from pya import Point, DPoint, Vector, DVector, DSimplePolygon, SimplePolygon, DPolygon, Polygon, \
    Region, Box, DBox
from pya import Trans, DTrans, CplxTrans, DCplxTrans, ICplxTrans
from pya import DEdge

import classLib

reload(classLib)
from classLib.chipDesign import ChipDesign
from classLib.shapes import Donut
from classLib.coplanars import CPW, Bridge1


class Intersection:
    @staticmethod
    def resolve_cpw_cpw_intersection(
        cpw1: CPW, cpw2: CPW, cpw_reg: Region, bridge_reg1: Region, bridge_reg2: Region,
        clearance_mul=1,
        debug_reg: Region = None
    ):
        """
        Places bridges between all grounds and slightly modifies geometry of the cutted coplanar
        edges in order to make a proper bridged CPW-CPW intersection

        Parameters
        ----------
        cpw1 : CPW
            was drawn first, so it's ends near the `cpw2` will be "repaired" a bit
        cpw2 : CPW
            was drawn second
        cpw_reg : Region
            Region where both cpws are going to be placed
        clearance_mul : float
            distance between cpw2 center line and cpw1 `repaired` edges centers in units of cpw2.b.
            (point to line distance)
            default = 1

        Returns
        -------

        Notes
        -------------
        Yes, repaired intersection could have its "repaired" cpw1's edges aligned along cpw2,
        but this is not that I have straight and simple solution to implement this.
        """
        edge1 = DEdge(cpw1.start, cpw1.end)
        edge2 = DEdge(cpw2.start, cpw2.end)
        intersection = edge1.intersection_point(edge2)
        if intersection is None:
            print(
                "Intersection resolve error, coplanars supplied does not intersect\n"
                "returning with None"
            )
            return None
        s1 = cpw1.dr / cpw1.dr.abs()
        s2 = cpw2.dr / cpw2.dr.abs()
        # 90 deg clockwise rotated tangent vectors
        rot90 = DCplxTrans(1, 90, False, 0, 0)
        n1 = rot90 * s1
        n2 = rot90 * s2

        # broken cpw endpoints
        # TODO: probably optimize intersection finding
        #  with `DEgde.cut_point()`, introduced in KLaout 0.27.1
        edge1 = DEdge(cpw1.start, cpw1.end)  # along the center of cpw1
        # edge of the one of the cpw2 grounds
        edge2 = DEdge(
            cpw2.start + n2 * (cpw2.gap + cpw2.width / 2),
            cpw2.end + n2 * (cpw2.gap + cpw2.width / 2)
        )
        # edge on the second of the cpw2 grounds
        edge3 = DEdge(cpw2.start - n2 * cpw2.b / 2, cpw2.end - n2 * cpw2.b / 2)

        broken_cpw_p1 = edge1.intersection_point(edge2)  # to the left from `s2`
        broken_cpw_p2 = edge1.intersection_point(edge3)  # to the right from `s2`
        # c1 = Donut(
        #     origin=broken_cpw_p1,
        #     inner_r=0,
        #     outer_r=cpw1.b/2
        # )
        # c1.place(debug_reg)
        # calculating points where broken cpw1 must be cutted along `s2` direction
        # shorted angle to rotate from `s1` to `s2` that lies in (-np.pi/2, np.pi/2)
        # `intersection_angle` > 0 if rotation from s1 to s2 is clockwise
        # and < 0 if anti-clockwise
        intersection_angle = np.arcsin(np.cross((s1.x, s1.y, 0), (s2.x, s2.y, 0))[-1])

        cut_repair_p1 = intersection - clearance_mul * cpw2.b / np.sin(intersection_angle) * s1
        cut_repair_p2 = intersection + clearance_mul * cpw2.b / np.sin(intersection_angle) * s1

        ''' Repairing broken cpw1 ends '''
        # filling ground for cpw2
        ground_filling_cpw1 = CPW(
            start=cut_repair_p2,
            end=cut_repair_p1,
            width=cpw1.b + 1,  # 1.1 to ensure filling is done well
            gap=0,
            open_end_gap=cpw1.b / 2,
            open_start_gap=cpw1.b / 2
        )
        ground_filling_cpw1.place(cpw_reg)
        # restoring `cpw2` state
        cpw2.place(cpw_reg)
        ''' Placing bridges '''
        # central cpw1 brigde
        bridge1 = Bridge1(
            center=intersection,
            gnd2gnd_dy=ground_filling_cpw1.dr.abs() + cpw1.b,
            trans_in=DCplxTrans(1, -90 + 180*np.arctan2(s1.y, s1.x)/np.pi, False, 0, 0)
        )
        bridge1.place(bridge_reg1, region_id="bridges_1")
        bridge1.place(bridge_reg2, region_id="bridges_2")

        bridge2 = Bridge1(
            center=intersection + 2*cpw1.b * s2,
            gnd2gnd_dy=ground_filling_cpw1.dr.abs() + cpw1.b,
            trans_in=DCplxTrans(1, -90 + 180 * np.arctan2(s1.y, s1.x) / np.pi, False, 0, 0)
        )
        bridge2.place(bridge_reg1, region_id="bridges_1")
        bridge2.place(bridge_reg2, region_id="bridges_2")

        bridge3 = Bridge1(
            center=intersection - 2*cpw1.b * s2,
            gnd2gnd_dy=ground_filling_cpw1.dr.abs() + cpw1.b,
            trans_in=DCplxTrans(1, -90 + 180 * np.arctan2(s1.y, s1.x) / np.pi, False, 0, 0)
        )
        bridge3.place(bridge_reg1, region_id="bridges_1")
        bridge3.place(bridge_reg2, region_id="bridges_2")


class MyDesign(ChipDesign):
    def __init__(self, cell_name: str = "testScript"):
        super().__init__(cell_name=cell_name)
        info_bridges1 = pya.LayerInfo(4, 0)  # bridge photo layer 1
        self.region_bridges1 = Region()
        self.layer_bridges1 = self.layout.layer(info_bridges1)

        info_bridges2 = pya.LayerInfo(5, 0)  # bridge photo layer 2
        self.region_bridges2 = Region()
        self.layer_bridges2 = self.layout.layer(info_bridges2)

    def draw(self):
        origin = DPoint(0, 0)
        self.region_ph.insert(DBox(DPoint(-5e6, -5e6), DPoint(5e6, 5e6)))
        self.region_bridges2.insert(DBox(DPoint(-5e6, -5e6), DPoint(5e6, 5e6)))

        cpw1 = CPW(
            start=DPoint(-2e6, -1e6),
            end=DPoint(3e6, 2e6),
            width=20e3,
            gap=10e3
        )
        cpw1.place(self.region_ph)
        cpw2 = CPW(
            start=DPoint(2e6, -1e6),
            end=DPoint(-3e6, 1e6),
            width=30e3,
            gap=15e3
        )
        cpw2.place(self.region_ph)

        Intersection.resolve_cpw_cpw_intersection(
            cpw1, cpw2, cpw_reg=self.region_ph,
            bridge_reg1=self.region_bridges1,
            bridge_reg2=self.region_bridges2,
            debug_reg=self.region_el
        )

    def _transfer_regs2cell(self):
        super()._transfer_regs2cell()
        self.cell.shapes(self.layer_bridges1).insert(self.region_bridges1)
        self.cell.shapes(self.layer_bridges2).insert(self.region_bridges2)


### MAIN FUNCTION ###
if __name__ == "__main__":
    my_design = MyDesign("testScript")
    my_design.draw()
    my_design.show()
