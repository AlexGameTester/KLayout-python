import os
import sys
import logging
import numpy as np
from datetime import datetime
from enum import Enum
from importlib import reload
from typing import List

from pya import Region, DBox, DPoint

import Projects.Kinemon.KmonDesign

logging.debug("Imported KmonDesign")
reload(Projects.Kinemon.KmonDesign)
logging.debug("Reloaded KmonDesign")
from Projects.Kinemon.KmonDesign import MeanderParams, KinemonParams
from Projects.Kinemon.KmonDesign import Kinemon

import classLib
from classLib.josJ import AsymSquidParams

logging.debug("Imported classLib from main")
reload(classLib)
logging.debug("Reloaded classLib from main")
from classLib.chipDesign import ChipDesign
from classLib.chipTemplates import CHIP_14x14_20pads
from classLib.coplanars import CPWParameters, CPW

# Configuring logging
LOGS_FOLDER = "c:/klayout_dev/logs/Design20pad/"
stdout_handler = logging.StreamHandler(sys.stdout)
file_handler = logging.FileHandler(
    filename=os.path.join(LOGS_FOLDER, datetime.now().strftime("%Y-%m-%d %H-%M-%S") + ".log"), mode='w')
logging.basicConfig(level=logging.DEBUG,
                    handlers=[stdout_handler, file_handler],
                    format="%(asctime)s %(levelname)s: %(message)s",
                    force=True)


class StartMode(Enum):
    """
    This class represents different script execution modes
    """
    SHOW = 0


class DefaultParams:
    """
    This class contains design parameters that are expected not to be changed during simulations, etc.
    """
    ro_line_Z = CPWParameters(width=18e3, gap=10e3)
    z_fl1 = CPWParameters(width=30e3, gap=15e3)

    asquid_dx = 35e3 / 2

    asquid_params = {
        # "band_ph_tol": 1e3,
        "squid_dx": 2 * asquid_dx,
        "squid_dy": 13e3,
        "TC_dx": 2.5e3 * np.sqrt(2) + 1e3,
        "TC_dy": 4e3,
        "TCW_dy": 0,
        "BCW_dy": 0e3,
        "BC_dx": [2.5e3 * np.sqrt(2) + 1e3],
        "BC_dy": 4e3,
    }

    contact_pad_params = {"back_metal_gap": 200e3,
                          "back_metal_width": 0e3,
                          "pad_length": 700e3,
                          "transition_len": 250e3
                          }
    meander_params = {
        'dr': DPoint(0, 0),
        'line_width_dx': 0.260e3,
        'line_width_dy': 0.120e3,
        'add_dx_mid': 8e3,
        'line_gap': 0.4e3
    }

    kinemon_params = {
        'area_ratio': 0.99,
        'MC_dy': None,
        'MC_dx': None,
        'KI_bridge_width': 1e3,
        'KI_bridge_height': 0.9e3,
        'KI_pad_y_offset': 0.2e3,
        'KI_pad_width': 3e3,
        'KI_ledge_y_offset': 0.3e3,
        'KI_JJ_ledge_height': 4e3,
        'KI_JJ_ledge_width': 2e3,
    }


class ProductionParams:
    """
    This class contains all variable parameters of the design and script execution parameters.
    """

    _instance = None

    @classmethod
    def instance(cls):
        if cls._instance is None:
            cls._instance = ProductionParams()
            logging.debug("ProductionParams instance is created")

        return cls._instance

    @property
    def chip_Z_list(self):
        return np.copy(self._chip_Z_list)

    @property
    def meander_length_list(self):
        return np.copy(self._meander_length_list)

    @property
    def big_jj_dx_list(self):
        return np.copy(self._big_jj_dx_list)

    @property
    def big_jj_dy_list(self):
        return np.copy(self._big_jj_dy_list)

    @property
    def small_jj_dx_list(self):
        return np.copy(self._small_jj_dx_list)

    @property
    def small_jj_dy_list(self):
        return np.copy(self._small_jj_dy_list)

    @property
    def ro_lines_ends_indexes(self):
        return np.copy(self._ro_lines_ends_indexes)

    def __init__(self):
        self.start_mode = StartMode.SHOW

        '''General geometry params'''
        self.NQUBITS = 8

        self._chip_Z_list = [
            DefaultParams.ro_line_Z, DefaultParams.z_fl1, DefaultParams.z_fl1, DefaultParams.z_fl1, DefaultParams.z_fl1,
            # left side
            DefaultParams.z_fl1, DefaultParams.z_fl1, DefaultParams.z_fl1, DefaultParams.z_fl1, DefaultParams.ro_line_Z,
            # bottom
            DefaultParams.ro_line_Z, DefaultParams.z_fl1, DefaultParams.z_fl1, DefaultParams.z_fl1, DefaultParams.z_fl1,
            # right
            DefaultParams.z_fl1, DefaultParams.z_fl1, DefaultParams.z_fl1, DefaultParams.z_fl1, DefaultParams.ro_line_Z
            # top
        ]

        # Indexes of contact pads that are connected to readout waveguides, grouped in pairs - start and end of each waveguide
        self._ro_lines_ends_indexes = np.array([(1, 13), (3, 11)])

        '''KinInd Meander params'''
        meander_length = 262949.689
        self._meander_length_list = np.array([meander_length] * self.NQUBITS)

        '''Josephson junction params'''
        self._big_jj_dx_list = np.array([120] * 8)
        self._big_jj_dy_list = np.array([
            287.65,
            168.00,
            236.45,
            311.81,
            264.64,
            377.67,
            170.30,
            292.83,
        ])
        self._small_jj_dx_list = np.array([
            114.85,
            114.85,
            114.85,
            114.85,
            119.66,
            141.90,
            119.66,
            119.66,
        ])
        self._small_jj_dy_list = np.array([
            104.85,
            104.85,
            104.85,
            104.85,
            109.66,
            131.90,
            109.66,
            109.66,
        ])



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

        logging.debug("Finished region initialization")
        self.lv.add_missing_layers()

        # Unchangeable parameters
        self.chip = CHIP_14x14_20pads
        self.chip_box = self.chip.box

        # None parameters initialization
        self.contact_pads = None
        self._init_contact_pads()

        self.readout_lines = None
        self._init_readout_lines()

        self.qubit_params_list = None
        self.qubit_origins = None
        self.qubit_trans = None
        self._init_qubit_params()

        self.qubits = None
        self._init_qubits()

    def _init_contact_pads(self):
        logging.debug("Initializing contact pads")
        self.contact_pads = CHIP_14x14_20pads.get_contact_pads(ProductionParams.instance().chip_Z_list,
                                                               **DefaultParams.contact_pad_params)

    def _init_readout_lines(self):
        logging.debug("Initializing readout lines")
        self.readout_lines = []
        for idx1, idx2 in ProductionParams.instance().ro_lines_ends_indexes:
            self.readout_lines.append(CPW(start=self.contact_pads[idx1].end,
                                          end=self.contact_pads[idx2].end,
                                          cpw_params=self.chip.chip_Z))

        if len(self.readout_lines) != len(ProductionParams.instance().ro_lines_ends_indexes):
            raise ValueError('')

    def _init_qubit_params(self):
        logging.debug("Initializing qubit parameters")
        meander_params_list = [MeanderParams(**DefaultParams.meander_params,
                                             line_length=length) for length in
                               ProductionParams.instance().meander_length_list]

        jj_params_iter = zip(ProductionParams.instance().big_jj_dx_list,
                             ProductionParams.instance().big_jj_dy_list,
                             ProductionParams.instance().small_jj_dx_list,
                             ProductionParams.instance().small_jj_dy_list)

        asquid_params_list = [AsymSquidParams(SQLTJJ_dx=big_jj_dx,
                                              SQLBJJ_dy=big_jj_dy,
                                              SQRTJJ_dx=small_jj_dx,
                                              SQRBJJ_dy=small_jj_dy,
                                              **DefaultParams.asquid_params)
                              for (big_jj_dx, big_jj_dy, small_jj_dx, small_jj_dy) in jj_params_iter]

        for asquid_params in asquid_params_list:
            asquid_params.bot_wire_x = [-DefaultParams.asquid_dx, DefaultParams.asquid_dx]
            asquid_params.BC_dx = [asquid_params.BC_dx[0]] * 2
            asquid_params.BCW_dx = [asquid_params.BCW_dx[0]] * 2

        if len(asquid_params_list) != ProductionParams.instance().NQUBITS:
            raise ValueError(
                f'Incorrect asquid_params_list length {len(asquid_params_list)}. NQUBITS = {ProductionParams.instance().NQUBITS}')

        if len(meander_params_list) != ProductionParams.instance().NQUBITS:
            raise ValueError(
                f'Incorrect meander_params_list length {len(meander_params_list)}. NQUBITS = {ProductionParams.instance().NQUBITS}')

        logging.debug("Preliminary list initialization finished. Creating qubit_params list")

        self.qubit_params_list = [
            KinemonParams(asym_squid_params=asquid_params,
                          meander_params=meander_params_list,
                          **DefaultParams.kinemon_params)
            for (asquid_params, meander_params_list) in zip(asquid_params_list, meander_params_list)]

        # TODO: Calculate qubit origins and trans here
        self.qubit_origins = [DPoint(0,0)] * ProductionParams.instance().NQUBITS

    def _init_qubits(self):
        logging.debug("Initializing qubits")
        self.qubits = [Kinemon(origin, kinemon_params)
                       for (origin, kinemon_params) in
                       zip(self.qubit_origins, self.qubit_params_list)]

    def draw(self, design_params=None):
        logging.debug("Drawing started")
        self.draw_chip()

        self.draw_contact_pads()

        self.draw_readout_lines()

        self.draw_qubits()

        logging.debug("Drawing completed")

    def draw_chip(self):
        logging.debug("Drawing chip box")
        self.region_bridges2.insert(self.chip_box)
        self.region_ph.insert(self.chip_box)

    def draw_contact_pads(self):
        logging.debug("Drawing contact pads")
        for i, contact_pad in enumerate(self.contact_pads):
            contact_pad.place(self.region_ph)

    def draw_readout_lines(self):
        for ro_line in self.readout_lines:
            ro_line.place(self.region_ph)

    def draw_qubits(self):
        logging.debug("Drawing qubits")
        for qubit in self.qubits:
            qubit.place(self.region_el, region_id="default")
            qubit.place(self.region_el, region_id="default_empty")
            qubit.place(self.region_kinInd, region_id="kinInd")

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
    if ProductionParams.instance().start_mode == StartMode.SHOW:
        logging.info("Executing script in showing mode")
        design = DesignKmon()
        design.draw()
        design.show()
