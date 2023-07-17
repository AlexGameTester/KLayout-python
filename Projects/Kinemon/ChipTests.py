from importlib import reload

from Projects.Dmon.Design import DesignDmon
from pya import DPoint
import Projects.Kinemon.KmonDesign
from Projects.Kinemon.KmonDesign import KinIndMeander, MeanderParams
from classLib import ChipDesign
from pya import Region
import pya



class TestDesign(DesignDmon):
    def __init__(self, cell_name="testScript"):
        super().__init__(cell_name)

        dc_bandage_layer_i = pya.LayerInfo(3,
                                           0)  # for DC contact deposition
        self.dc_bandage_reg = Region()
        self.dc_bandage_layer = self.layout.layer(dc_bandage_layer_i)

        info_bridges1 = pya.LayerInfo(4, 0)  # bridge photo layer 1
        self.region_bridges1 = Region()
        self.layer_bridges1 = self.layout.layer(info_bridges1)

        info_bridges2 = pya.LayerInfo(5, 0)  # bridge photo layer 2
        self.region_bridges2 = Region()
        self.layer_bridges2 = self.layout.layer(info_bridges2)

        # layer with polygons that will protect structures located
        # on the `self.region_el` - e-beam litography layer
        info_el_protection = pya.LayerInfo(6, 0)
        self.region_el_protection = Region()
        self.layer_el_protection = self.layout.layer(info_el_protection)

        info_kinInd_layer = pya.LayerInfo(7, 0)
        self.region_kinInd = Region()
        self.layer_kinInd = self.layout.layer(info_kinInd_layer)

        # has to call it once more to add new layers
        self.lv.add_missing_layers()

    def draw_kin_ind(self):
        meander_params = MeanderParams(DPoint(0, 5e3), 10, 100, 5, 5, 20e3, 0)
        self.kin_ind = KinIndMeander(meander_params, region_id="kinInd")
        self.kin_ind.place(self.region_kinInd, region_id="kinInd")

    def draw(self, design_params=None):
        self.draw_chip()

        self.draw_readout_waveguide()

        self.draw_kin_ind()




if __name__ == "__main__":
    reload(Projects.Kinemon.KmonDesign)
    test_design = TestDesign("testScript")
    test_design.draw()
    test_design.show()