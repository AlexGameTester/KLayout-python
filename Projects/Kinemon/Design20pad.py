import os
import sys
import logging
from datetime import datetime
from enum import Enum
from importlib import reload
from typing import List

from pya import Region, DBox

import classLib

logging.debug("imported classLib from main")
reload(classLib)
logging.debug("reloaded classLib from main")
from classLib.chipDesign import ChipDesign
from classLib.chipTemplates import CHIP_14x14_20pads


# Configuring logging
LOGS_FILEPATH = "c:/klayout_dev/logs/Design20pad/"
stdout_handler = logging.StreamHandler(sys.stdout)
file_handler = logging.FileHandler(filename=os.path.join(LOGS_FILEPATH, datetime.now().strftime("%Y-%m-%d %H-%M-%S") + ".log"), mode='w')
logging.basicConfig(level=logging.DEBUG,
                    handlers=[stdout_handler, file_handler],
                    format="%(asctime)s %(levelname)s: %(message)s",
                    force=True)

class StartMode(Enum):
    """
    This class represents different script execution modes
    """
    SHOW = 0


class ProductionParams:
    """
    This class contains all variable parameters of the design and script execution parameters.
    """
    start_mode = StartMode.SHOW


class DesignKmon(ChipDesign):
    """
    This class represents a 14x14 chip with kinemons that are grouped into one, two and three qubit schemes.
    """

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

        self.regions: List[Region] = [
            self.region_ph, self.region_bridges1, self.region_bridges2,
            self.region_el, self.dc_bandage_reg,
            self.region_el_protection, self.region_cut, self.region_kinInd
        ]
        # has to call it once more to add new layers
        self.lv.add_missing_layers()

        # Unchangeable parameters
        self.chip = CHIP_14x14_20pads
        self.chip_box = self.chip.box

        # None parameters initialization
        self.contact_pads = None
        self._init_contact_pads()



    def _init_contact_pads(self):
        logging.debug("Creating contact pads")
        self.contact_pads = CHIP_14x14_20pads.get_contact_pads()

    def draw(self, design_params=None):
        self.draw_chip()

        self.draw_contact_pads()

        logging.debug("Drawing completed")

    def draw_chip(self):
        logging.debug("Drawing chip box")
        self.region_bridges2.insert(self.chip_box)
        self.region_ph.insert(self.chip_box)

    def draw_contact_pads(self):
        logging.debug("Drawing contact pads")
        for i, contact_pad in enumerate(self.contact_pads):
            contact_pad.place(self.region_ph)

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
        logging.info("Executing script in showing mode")
        design = DesignKmon()
        design.draw()
        design.show()
