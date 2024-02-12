from importlib import reload
from Projects.Kinemon.DesignCopy import DesignDmon
from pya import DPoint
import Projects.Kinemon.KmonDesign
from Projects.Kinemon.KmonDesign import KinIndMeander, MeanderParams, Kinemon
from classLib import ChipDesign
from pya import Region, DCplxTrans, DVector
import pya

from classLib.helpers.polygon_splitting import split_polygons

from Projects.Kinemon.DesignCopy import DefaultParams


class TestDesign(DesignDmon):
    def __init__(self, cell_name="testScript"):
        super().__init__(cell_name)


    def draw_kin_ind(self):
        meander_params = MeanderParams(**DefaultParams.meander_params_dict, line_length=3e5)
        meander_params.dr = DVector(0, 20e3)
        trans = DCplxTrans(DVector(3e3, 0))
        self.kin_ind = KinIndMeander(meander_params, trans_in=trans, region_id="kinInd")
        self.kin_ind.place(self.region_kinInd, region_id="kinInd")

    def split_polygons_in_layers(self, max_pts=200):
        super().split_polygons_in_layers(max_pts)
        self.region_kinInd = split_polygons(self.region_kinInd, max_pts)

        for poly in self.region_kinInd:
            if poly.num_points() > max_pts:
                print(f'Splitting kinetic inductance region is not done. Problem with polygon: {poly}')

    def draw_kinemon(self):
        kinemon = Kinemon(DPoint(10e3, 50e3), )
    def draw(self, design_params=None):
        self.draw_chip()

        # self.draw_readout_waveguide()

        self.draw_kin_ind()

        self.split_polygons_in_layers(max_pts=180)


if __name__ == "__main__":
    reload(Projects.Kinemon.KmonDesign)
    test_design = TestDesign("testScript")
    test_design.draw()
