from enum import Enum

from pya import Region, DBox

from classLib import ChipDesign
from classLib.chipTemplates import CHIP_14x14_20pads


class StartMode(Enum):
   SHOW = 0
class ProductionParams:
    start_mode = StartMode.SHOW

class DesignKmon(ChipDesign):
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
        # on the `self.region_el` - e-beam lithography layer
        info_el_protection = pya.LayerInfo(6, 0)
        self.region_el_protection = Region()
        self.layer_el_protection = self.layout.layer(info_el_protection)

        info_kinInd_layer = pya.LayerInfo(7, 0)
        self.region_kinInd = Region()
        self.layer_kinInd = self.layout.layer(info_kinInd_layer)

        # has to call it once more to add new layers
        self.lv.add_missing_layers()

        # Unchangeable parameters
        self.chip = CHIP_14x14_20pads
        self.chip_box = self.chip.box

    def draw(self, design_params=None):
        self.draw_chip()

        self.draw_contact_pads()

    def draw_chip(self):
        self.region_bridges2.insert(self.chip_box)
        self.region_ph.insert(self.chip_box)

    def draw_contact_pads(self):
        pass

    def _transfer_regs2cell(self):
        # this too methods assumes that all previous drawing
        # functions are placing their object on regions
        # in order to avoid extensive copying of the polygons
        # to/from cell.shapes during the logic operations on
        # polygons
        self.cell.shapes(self.layer_ph).insert(self.region_ph)
        self.cell.shapes(self.layer_el).insert(self.region_el)
        self.cell.shapes(self.dc_bandage_layer).insert(self.dc_bandage_reg)
        self.cell.shapes(self.layer_bridges1).insert(self.region_bridges1)
        self.cell.shapes(self.layer_bridges2).insert(self.region_bridges2)
        self.cell.shapes(self.layer_el_protection).insert(
            self.region_el_protection)
        self.cell.shapes(self.layer_kinInd).insert(self.region_kinInd)
        self.lv.zoom_fit()


if __name__ == "__main__":
    if ProductionParams.start_mode == StartMode.SHOW:
        design = DesignKmon()
        design.draw()
        design.show()



