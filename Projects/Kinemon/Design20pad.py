import os
import sys
import logging
import numpy as np
from datetime import datetime
from enum import Enum
from importlib import reload
from typing import List

from pya import Region, DBox, DPoint, DVector, Trans

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
from classLib.capacitors import XmonCross
from classLib.resonators import EMResonatorTL3QbitWormRLTailXmonFork

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

    Z_res = CPWParameters(10e3, 6e3)
    res_r = 60e3
    n_coils = 3
    res_tail_shape = "LRLRL"

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

    @property
    def resonator_above_line_list(self):
        return np.copy(self._resonator_above_line_list)

    @property
    def cross_len_x_list(self):
        return np.copy(self._cross_len_x_list)

    @property
    def cross_width_x_list(self):
        return np.copy(self._cross_width_x_list)

    @property
    def cross_gnd_gap_x_list(self):
        return np.copy(self._cross_gnd_gap_x_list)

    @property
    def cross_len_y_list(self):
        return np.copy(self._cross_len_y_list)

    @property
    def cross_width_y_list(self):
        return np.copy(self._cross_width_y_list)

    @property
    def cross_gnd_gap_y_list(self):
        return np.copy(self._cross_gnd_gap_y_list)

    @property
    def cross_gnd_gap_face_y_list(self):
        return np.copy(self._cross_gnd_gap_face_y_list)

    @property
    def resonator_L_coupling_list(self):
        return np.copy(self._resonator_L_coupling_list)

    @property
    def resonator_origin_on_line_list(self):
        return np.copy(self._resonator_origin_on_line_list)

    @property
    def resonator_trans_list(self):
        return np.copy(self._resonator_trans_list)

    @property
    def resonator_above_line_list(self):
        return np.copy(self._resonator_above_line_list)

    @property
    def resonator_L0_list(self):
        return np.copy(self._resonator_L0_list)

    @property
    def resonator_L1_list(self):
        return np.copy(self._resonator_L1_list)

    @property
    def fork_metal_width_list(self):
        return np.copy(self._fork_metal_width_list)

    @property
    def fork_gnd_gap_list(self):
        return np.copy(self._fork_gnd_gap_list)

    @property
    def xmon_fork_gnd_gap_list(self):
        return np.copy(self._xmon_fork_gnd_gap_list)

    @property
    def fork_y_span_list(self):
        return np.copy(self._fork_y_span_list)

    resonator_to_line_list = property(fget=lambda self: np.copy(self._resonator_to_line_list))
    xmon_res_d_list = property(fget=lambda self: np.copy(self._xmon_res_d_list))

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

        '''Resonator params'''
        # It's not real origin of resonator, it'll be shifted according to to_line parameter.
        self._resonator_origin_on_line_list = np.array([DPoint(1.5e6 * i, 0) for i in range(self.NQUBITS)])
        self._resonator_trans_list = np.array([Trans.R0] * self.NQUBITS)
        # Additional variable that shows if resonator is above or below corresponding readout line. This variable is defined by corresponding trans
        self._resonator_above_line_list = np.array([True] * self.NQUBITS)

        self._resonator_L_coupling_list = np.array([
            1e3 * x for x in [310, 320, 320, 310] * 2
        ])

        L0 = 986e3
        # long vertical line length
        self._resonator_L0_list = [L0] * self.NQUBITS

        self._resonator_L1_list = np.array(
            [141080.2779, 132432.5366, 125948.9511, 106176.2057,
             115540.3816, 89988.1556, 77104.8224, 46724.314]
        )

        self._fork_metal_width_list = np.array(
            [1e3 * x for x in ([10] * 3 + [6] * 2 + [10] * 3)])
        self._fork_gnd_gap_list = np.array([10e3] * self.NQUBITS)
        self._xmon_fork_gnd_gap_list = np.array([10e3] * self.NQUBITS)
        self._fork_y_span_list = np.array(
            [
                x * 1e3 for x in
                [42.0, 47, 52.0, 64.0, 70.0, 62.0, 68.0, 81.0]
            ]
        )

        self._resonator_to_line_list = np.array([45e3] * 8)

        '''X-mon params'''
        self._cross_len_x_list = np.array(
            [1e3 * x for x in [129.905, 65.602, 35.098, 0, 0, 0,
                               0, 0]]
        )
        self._cross_width_x_list = np.array(
            [1e3 * x for x in [16, 16, 16, 16, 16, 32, 56, 56]]
        )
        self._cross_gnd_gap_x_list = np.array([60e3] * self.NQUBITS)
        self._cross_len_y_list = np.array(
            [1e3 * x for x in
             [170.0, 232.0, 281.0, 256.0, 267.0, 186.0, 195.0, 223.0]]
        )
        self._cross_width_y_list = np.array(
            [1e3 * x for x in [16, 16, 16, 16, 16, 32, 32, 32]]
        )
        self._cross_gnd_gap_y_list = np.array([60e3] * self.NQUBITS)
        self._cross_gnd_gap_face_y_list = np.array([20e3] * self.NQUBITS)

        self._xmon_res_d_list = [40e3] * self.NQUBITS

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

    def _check_assertions_list(self, lst):
        assert isinstance(lst, np.ndarray)
        assert len(lst) == self.NQUBITS
        assert np.all(lst >= 0)

    def check_assertions(self):
        # TODO: Maybe add assert(all > 0) for the lists here
        logging.info("Checking assertions in ProductionParams")
        assert isinstance(self.start_mode, StartMode)

        assert isinstance(self.NQUBITS, int)
        assert self.NQUBITS > 0

        assert len(self._chip_Z_list) == 20

        assert all([len(ind_pair) == 2 and all([0 <= ind < 20 for ind in ind_pair]) for ind_pair in
                    self.ro_lines_ends_indexes])

        assert len(self._meander_length_list) == self.NQUBITS

        assert len(self._resonator_trans_list) == self.NQUBITS
        assert len(self._resonator_above_line_list) == self.NQUBITS

        assert len(self._cross_len_x_list) == self.NQUBITS
        assert len(self._cross_width_x_list) == self.NQUBITS
        assert len(self._cross_gnd_gap_x_list) == self.NQUBITS
        assert len(self._cross_len_y_list) == self.NQUBITS
        assert len(self._cross_width_y_list) == self.NQUBITS
        assert len(self._cross_gnd_gap_y_list) == self.NQUBITS
        assert len(self._cross_gnd_gap_face_y_list) == self.NQUBITS

        assert len(self._big_jj_dx_list) == self.NQUBITS
        assert len(self._big_jj_dy_list) == self.NQUBITS
        assert len(self._small_jj_dx_list) == self.NQUBITS
        assert len(self._small_jj_dx_list) == self.NQUBITS

        logging.info("Assertion check completed")


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

        self.resonators = None
        self._init_resonators()

        self.xmons = None
        self._init_xmons()

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

    def _init_resonators(self):
        logging.debug("Initializing resonators")
        self.resonators = []
        self._L2_list = [DefaultParams.res_r] * ProductionParams.instance().NQUBITS
        # horizontal line connected to L2
        self._L3_list = [2 * DefaultParams.res_r] * ProductionParams.instance().NQUBITS
        # vertical line connected to L3
        self._L4_list = [DefaultParams.res_r] * ProductionParams.instance().NQUBITS

        self._tail_segment_lengths_list = [[L2, L3, L4]
                                           for L2, L3, L4 in
                                           zip(self._L2_list, self._L3_list,
                                               self._L4_list)]
        self._tail_turn_angles_list = [
                                          [np.pi / 2, -np.pi / 2]
                                      ] * ProductionParams.instance().NQUBITS
        self._tail_trans_in_list = [Trans.R270] * ProductionParams.instance().NQUBITS

        self._fork_x_span_list = ProductionParams.instance().cross_width_y_list + 2 * \
                                 (
                                         ProductionParams.instance().xmon_fork_gnd_gap_list +
                                         ProductionParams.instance().fork_metal_width_list)

        for (origin_on_line,
             to_line,
             above_line,
             trans,
             L_coupling,
             L0,
             L1,
             tail_segment_length,
             tail_turn_angle,
             tail_trans_in,
             fork_metal_width,
             fork_x_span,
             fork_y_span,
             fork_gnd_gap) in zip(ProductionParams.instance().resonator_origin_on_line_list,
                                  ProductionParams.instance().resonator_to_line_list,
                                  ProductionParams.instance().resonator_above_line_list,
                                  ProductionParams.instance().resonator_trans_list,
                                  ProductionParams.instance().resonator_L_coupling_list,
                                  ProductionParams.instance().resonator_L0_list,
                                  ProductionParams.instance().resonator_L1_list,
                                  self._tail_segment_lengths_list,
                                  self._tail_turn_angles_list,
                                  self._tail_trans_in_list,
                                  ProductionParams.instance().fork_metal_width_list,
                                  self._fork_x_span_list,
                                  ProductionParams.instance().fork_y_span_list,
                                  ProductionParams.instance().fork_gnd_gap_list):
            self.resonators.append(
                EMResonatorTL3QbitWormRLTailXmonFork(
                    Z0=DefaultParams.Z_res,
                    start=origin_on_line + DVector(0, to_line if above_line else -to_line),
                    L_coupling=L_coupling,
                    L0=L0,
                    L1=L1,
                    r=DefaultParams.res_r,
                    N=DefaultParams.n_coils,
                    tail_shape=DefaultParams.res_tail_shape,
                    tail_turn_radiuses=DefaultParams.res_r,
                    tail_segment_lengths=tail_segment_length,
                    tail_turn_angles=tail_turn_angle,
                    tail_trans_in=tail_trans_in,
                    fork_x_span=fork_x_span,
                    fork_y_span=fork_y_span,
                    fork_metal_width=fork_metal_width,
                    fork_gnd_gap=fork_gnd_gap,
                    trans_in=trans
                )
            )

    def _init_xmons(self):
        logging.debug("Initializing X-mons")

        xmon_params_packed = [
            {
                "sideX_length": cross_len_x,
                "sideX_width": cross_width_x,
                "sideX_gnd_gap": cross_gnd_gap_x,
                "sideY_length": cross_len_y,
                "sideY_width": cross_width_y,
                "sideY_gnd_gap": cross_gnd_gap_y,
                "sideX_face_gnd_gap": cross_gnd_gap_x,
                "sideY_face_gnd_gap": cross_gnd_gap_face_y
            }
            for
            (cross_len_x,
             cross_width_x,
             cross_gnd_gap_x,
             cross_len_y,
             cross_width_y,
             cross_gnd_gap_y,
             cross_gnd_gap_face_y)
            in
            zip(ProductionParams.instance().cross_len_x_list,
                ProductionParams.instance().cross_width_x_list,
                ProductionParams.instance().cross_gnd_gap_x_list,
                ProductionParams.instance().cross_len_y_list,
                ProductionParams.instance().cross_width_y_list,
                ProductionParams.instance().cross_gnd_gap_y_list,
                ProductionParams.instance().cross_gnd_gap_face_y_list,
                )
        ]

        data_iter = zip(self.resonators,
                        ProductionParams.instance().resonator_trans_list,
                        ProductionParams.instance().resonator_above_line_list,
                        ProductionParams.instance().xmon_res_d_list,
                        ProductionParams.instance().fork_metal_width_list,
                        ProductionParams.instance().fork_gnd_gap_list,
                        xmon_params_packed)
        self.xmons = []

        for resonator, trans, above_line, xmon_res_d, fork_metal_width, fork_gnd_gap, xmon_params_dict in data_iter:
            # d = self.xmon_res_d_list[res_idx] + xmon_params_dict['sideY_length']
            #      + xmon_params_dict['sideX_width'] / 2 + \
            #     self.fork_metal_width_list[res_idx] + self.fork_gnd_gap
            d = xmon_res_d + xmon_params_dict['sideY_length'] + xmon_params_dict[
                'sideX_width'] / 2 + fork_metal_width + fork_gnd_gap
            if above_line:
                xmon_center = resonator.end + DVector(0, -d)
            else:
                xmon_center = resonator.end + DVector(0, d)

            self.xmons.append(
                XmonCross(
                    origin=xmon_center,
                    trans_in=trans,
                    **xmon_params_dict
                )
            )

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

        self.qubit_origins = [((xmon.cpw_bempt.end + xmon.cpw_bempt.start) / 2 if not above_line
                              else xmon.cpw_tempt.end + xmon.cpw_tempt.start) / 2
                              for above_line, xmon
                              in zip(ProductionParams.instance().resonator_above_line_list,
                                     self.xmons)]

        self.qubit_trans = [trans * Trans.R180 if above_line
                            else trans
                            for above_line, trans
                            in zip(ProductionParams.instance().resonator_above_line_list,
                                   ProductionParams.instance().resonator_trans_list)]

    def _init_qubits(self):
        logging.debug("Initializing qubits")
        self.qubits = [Kinemon(origin, kinemon_params, trans)
                       for (origin, kinemon_params, trans) in
                       zip(self.qubit_origins, self.qubit_params_list, self.qubit_trans)]

    def draw(self, design_params=None):
        logging.debug("Drawing started")
        self.draw_chip()

        self.draw_contact_pads()

        self.draw_readout_lines()

        self.draw_resonators()

        self.draw_xmons()

        # Resonators are need to be drawn again as xmons' CPWs can cut parts of resonators' forks
        logging.debug("Redrawing resonators")
        self.draw_resonators()

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
        logging.debug("Drawing readout lines")
        for ro_line in self.readout_lines:
            ro_line.place(self.region_ph)

    def draw_resonators(self):
        logging.debug("Drawing resonators")
        for resonator in self.resonators:
            resonator.place(self.region_ph)

    def draw_xmons(self):
        logging.debug("Drawing X-mons")
        for xmon in self.xmons:
            xmon.place(self.region_ph)

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
    ProductionParams.instance().check_assertions()
    if ProductionParams.instance().start_mode == StartMode.SHOW:
        logging.info("Executing script in showing mode")
        design = DesignKmon()
        design.draw()
        design.show()
