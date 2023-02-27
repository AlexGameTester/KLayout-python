import pya
from pya import Point, DPoint, Vector, DVector, DSimplePolygon, SimplePolygon, DPolygon, Polygon, Region, Box, DBox
from pya import Trans, DTrans, CplxTrans, DCplxTrans, ICplxTrans

from importlib import reload

import classLib

reload(classLib)
from classLib.chipDesign import ChipDesign
from classLib.shapes import TmonT, Rectangle, XmonCross, Circle
from classLib.coplanars import CPWRLPath, CPWParameters, CPW
from classLib.baseClasses import ElementBase, ComplexBase
from classLib.chipTemplates import CHIP_16p5x16p5_20pads
import numpy as np
from collections import OrderedDict
import os
import csv
import shutil
from pathlib import Path
import sonnetSim

reload(sonnetSim)
from sonnetSim import SonnetLab, SonnetPort, SimulationBox

refractive_index = np.sqrt(6.26423)
PROJECT_DIR = os.path.dirname(__file__)


class CircleResonator(ComplexBase):

    def __init__(self, origin, inner_r, outer_r, outer_gap, z, L, r, N, n_pts=50, trans_in=None):
        self.origin = origin
        self.inner_r = inner_r
        self.outer_r = outer_r
        self.outer_gap = outer_gap
        self.z = z
        self.L = L
        self.r = r
        self.N = N
        self.n_pts = n_pts

        super().__init__(origin, trans_in)

    def init_primitives(self):
        origin = DPoint(0, 0)
        if self.N == 0:
            shape = 'L'
            lengths = [self.inner_r + self.outer_r + self.outer_gap]
            rads = []
            angles = []
        else:
            shape = 'LRL' + 'RRLRRL' * self.N + 'RRLRL'
            inner_subwire_len = self.inner_r - self.r * (2*self.N + 1)
            lengths = [
                self.outer_gap + self.outer_r - self.inner_r + inner_subwire_len,
                self.L/2
            ] + [self.L] * 2 * self.N + [
                self.L/2,
                inner_subwire_len
            ]
            rads = [self.r] * 4*(self.N + 1)
            angles = [-np.pi/2, np.pi/2, np.pi/2] + [-np.pi/2, -np.pi/2, np.pi/2, np.pi/2] * self.N + [-np.pi/2]

        # placing circle r_out with dr clearance from ground polygon
        self.empt_circle = Circle(origin, self.outer_r + self.outer_gap, n_pts=self.n_pts, inverse=True)
        self.outer_circle = Circle(origin, self.outer_r, n_pts=self.n_pts, inverse=False)
        self.inner_circle = Circle(origin, self.inner_r, n_pts=self.n_pts, inverse=True)
        self.empt_circle.empty_region -= self.outer_circle.metal_region
        self.empt_circle.empty_region += self.inner_circle.metal_region
        self.primitives["empt_circle"] = self.empt_circle
        self.primitives["outer_circle"] = self.outer_circle
        self.primitives["inner_circle"] = self.inner_circle

        self.inductive_cop = CPWRLPath(
            DPoint(0, -self.outer_r - self.outer_gap), shape, self.z, rads, lengths, angles, trans_in=Trans.R90
        )
        self.primitives["inductive_cop"] = self.inductive_cop


class CircleResonatorDesign(ChipDesign):
    center: DPoint
    cop_pars: CPWParameters
    resonator: CircleResonator

    design_dict = OrderedDict()

    last_draw = None

    def __init__(self, cell_name="testScript"):
        super().__init__(cell_name)

        self.center = DPoint(16500e3 / 2, 16500e3 / 2)
        self.main_axe = self.center.x

        self.cop_rad = 100e3

        self.ind_cop_pars = CPWParameters(4e3, 4e3)
        self.cop_pars = CPWParameters(20e3, 12e3)
        self.cap_plate_pars = CPWParameters(10e3, 6e3)

        self.chip = CHIP_16p5x16p5_20pads
        self.chip_box: pya.DBox = self.chip.box

        self.res_rad = 160e3
        self.res_gap = 20e3
        self.res_in_rad = 100e3

        self.cap_plt_len = 60e3

    def set_parts(self):
        self.resonator = CircleResonator(
            self.center, self.res_in_rad, self.res_rad, self.res_gap, self.ind_cop_pars, 80e3, 10e3, 3
        )

        self.line_1 = CPW(
            start=self.center - DVector(self.res_rad - 1e3, 0),
            end=self.center - DVector(self.res_rad + 100e3, 0),
            cpw_params=self.cop_pars
        )
        # left capacitor
        self.cap_plates = [
            CPW(
                start=self.center + DVector((self.res_rad + self.res_in_rad) / 2, self.cap_plt_len / 2),
                end=self.center + DVector((self.res_rad + self.res_in_rad) / 2, self.cap_plt_len / 2 + self.cap_plate_pars.gap),
                width=0, gap=self.cap_plate_pars.width/2 + self.cap_plate_pars.gap
            ),
            CPW(
                start=self.center + DVector((self.res_rad+self.res_in_rad)/2, self.cap_plt_len/2),
                end=self.center + DVector((self.res_rad+self.res_in_rad)/2, -self.cap_plt_len/2),
                cpw_params=self.cap_plate_pars
            ),
            CPW(
                start=self.center + DVector((self.res_rad + self.res_in_rad) / 2, -self.cap_plt_len / 2),
                end=self.center + DVector((self.res_rad + self.res_in_rad) / 2,
                                          -self.cap_plt_len / 2 - self.cap_plate_pars.gap),
                width=0, gap=self.cap_plate_pars.width / 2 + self.cap_plate_pars.gap
            )
        ]

        self.line_2 = CPW(
            start=self.center + DVector((self.res_rad+self.res_in_rad)/2 + self.cap_plate_pars.width/2, 0),
            end=self.center + DVector(self.res_rad + 100e3, 0),
            cpw_params=self.cop_pars
        )

    def draw(self, parts=None):
        self.set_parts()
        self.region_ph.insert(self.chip_box)

        self.resonator.place(self.region_ph)
        self.line_1.place(self.region_ph)
        for cpw in self.cap_plates:
            cpw.place(self.region_ph)
        self.line_2.place(self.region_ph)


def simulate_Cfork():
    resolution_dx = 2e3
    resolution_dy = 2e3

    ### DRAWING SECTION START ###
    design = CircleResonatorDesign("testScript")

    design.draw({'transmon', 'resonator'})

    worm = design.resonator
    xmonCross = design.transmon
    worm_start = list(worm.primitives.values())[0].start

    # draw open end at the resonators start
    p1 = worm_start - DVector(design.cop_pars.b / 2, 0)
    rec = Rectangle(p1, design.cop_pars.b, design.cop_pars.b / 2,
                    inverse=True)
    rec.place(design.region_ph)

    if worm_start.x < xmonCross.center.x:
        dr = (worm_start - xmonCross.cpw_r.end)
    else:
        dr = (worm_start - xmonCross.cpw_l.end)
    dr.x = abs(dr.x)
    dr.y = abs(dr.y)

    box_side_x = 4 * xmonCross.sideX_length
    box_side_y = 4 * xmonCross.sideY_length
    dv = DVector(box_side_x / 2, box_side_y / 2)

    crop_box = pya.Box().from_dbox(pya.Box(
        xmonCross.center + dv,
        xmonCross.center + (-1) * dv
    ))
    design.crop(crop_box)
    dr = DPoint(0, 0) - crop_box.p1

    # finding the furthest edge of cropped resonator`s central line polygon
    # sonnet port will be attached to this edge
    reg1 = worm.metal_region & Region(crop_box)
    reg1.merge()
    max_distance = 0
    port_pt = None
    for poly in reg1.each():
        for edge in poly.each_edge():
            edge_center = (edge.p1 + edge.p2) / 2
            dp = edge_center - xmonCross.cpw_b.end
            d = max(abs(dp.x), abs(dp.y))
            if d > max_distance:
                max_distance = d
                port_pt = edge_center
    design.sonnet_ports.append(port_pt)
    design.sonnet_ports.append(xmonCross.cpw_l.end)

    design.transform_region(design.region_ph, DTrans(dr.x, dr.y),
                            trans_ports=True)

    design.show()
    design.lv.zoom_fit()
    ### DRAWING SECTION END ###

    for prt in design.sonnet_ports:
        print(prt)

    simulate_C12(crop_box, design, filename="ps_Cqr_results.csv")


def simulate_C12(crop_box, design, filename='c12.csv', resolution_dx=2e3, resolution_dy=2e3):
    ### SIMULATION SECTION START ###
    ml_terminal = SonnetLab()
    # print("starting connection...")
    from sonnetSim.cMD import CMD

    ml_terminal._send(CMD.SAY_HELLO)
    ml_terminal.clear()
    simBox = SimulationBox(
        crop_box.width(),
        crop_box.height(),
        crop_box.width() / resolution_dx,
        crop_box.height() / resolution_dy
    )

    ml_terminal.set_boxProps(simBox)
    # print("sending cell and layer")
    from sonnetSim.pORT_TYPES import PORT_TYPES

    ports = [
        SonnetPort(design.sonnet_ports[0], PORT_TYPES.AUTOGROUNDED),
        SonnetPort(design.sonnet_ports[1], PORT_TYPES.AUTOGROUNDED)
    ]
    # for sp in ports:
    #     print(sp.point)
    ml_terminal.set_ports(ports)

    ml_terminal.send_polygons(design.cell, design.layer_ph)
    ml_terminal.set_linspace_sweep(7, 7, 1)
    print("simulating...")
    result_path = ml_terminal.start_simulation(wait=True)
    ml_terminal.release()

    ### SIMULATION SECTION END ###

    ### CALCULATE C_QR CAPACITANCE SECTION START ###
    C12 = None
    with open(result_path.decode("ascii"), "r") as csv_file:
        data_rows = list(csv.reader(csv_file))
        ports_imps_row = data_rows[6]
        R = float(ports_imps_row[0].split(' ')[1])
        data_row = data_rows[8]
        freq0 = float(data_row[0])

        s = [[0, 0], [0, 0]]  # s-matrix
        # print(data_row)
        for i in range(0, 2):
            for j in range(0, 2):
                s[i][j] = complex(float(data_row[1 + 2 * (i * 2 + j)]),
                                  float(data_row[
                                            1 + 2 * (i * 2 + j) + 1]))
        for i in range(0, 2):
            for j in range(0, 2):
                s[i][j] = complex(
                    float(data_row[1 + 2 * (i * 2 + j)]),
                    float(data_row[1 + 2 * (i * 2 + j) + 1]))

        import math
        delta = (1 + s[0][0]) * (1 + s[1][1]) - s[0][1] * s[1][0]
        y11 = 1 / R * ((1 - s[0][0]) * (1 + s[1][1]) + s[0][1] * s[1][0]) / delta
        y22 = 1 / R * ((1 - s[1][1]) * (1 + s[0][0]) + s[0][1] * s[1][0]) / delta
        C1 = -1e15 / (2 * math.pi * freq0 * 1e9 * (1 / y11).imag)
        C2 = -1e15 / (2 * math.pi * freq0 * 1e9 * (1 / y22).imag)
        # formula taken from https://en.wikipedia.org/wiki/Admittance_parameters#Two_port
        y21 = -2 * s[1][0] / delta * 1 / R
        C12 = 1e15 / (2 * math.pi * freq0 * 1e9 * (1 / y21).imag)
        C1 -= C12
        C2 -= C12

        ### CALCULATE C_QR CAPACITANCE SECTION START ###

        ### SAVING REUSLTS SECTION START ###
        design.layout.write(
            str(PROJECT_DIR / "purcell_filter" / "final_design" / "capacitance" / (filename[:-4] + '.gds')))
        output_filepath = PROJECT_DIR / "purcell_filter" / "final_design" / "capacitance" / filename
        if os.path.exists(str(output_filepath)):
            # append data to file
            with open(str(output_filepath), "a", newline='') as csv_file:
                writer = csv.writer(csv_file)
                writer.writerow(
                    [C12, C1, C2]
                )
        else:
            # create file, add header, append data
            with open(str(output_filepath), "w", newline='') as csv_file:
                writer = csv.writer(csv_file)
                # create header of the file
                writer.writerow(
                    ["C12, fF", "C1, fF", "C2, fF"])
                writer.writerow(
                    [C12, C1, C2]
                )

        design.save_logdata(PROJECT_DIR / "purcell_filter" / "final_design" / "capacitance" / (filename[:-4] + '.xml'),
                            comment=f'Capacitance simulation results saved at {filename}.')

        ### SAVING RESULTS SECTION END ###


def simulate_full_S_pars(filename='cs_S12.csv', min_freq=6.8, max_freq=7.8):
    ### DRAWING SECTION START ###
    design = CircleResonatorDesign("testScript")

    design.draw()

    p1 = design.line_1.end
    p2 = design.line_2.end

    box_side_x = (p1-p2).x
    box_side_y = box_side_x
    dv = DPoint(box_side_x / 2, box_side_y / 2)

    crop_box = pya.Box().from_dbox(pya.Box(
        design.center + dv,
        design.center + (-1) * dv
    ))
    design.crop(crop_box)
    dr = DPoint(0, 0) - crop_box.p1

    design.sonnet_ports.append(p1)
    design.sonnet_ports.append(p2)

    design.transform_region(design.region_ph, DTrans(dr.x, dr.y),
                            trans_ports=True)

    design.show()
    design.lv.zoom_fit()
    ### DRAWING SECTION END ###

    for prt in design.sonnet_ports:
        print(prt)

    simulate_S_pars(design, crop_box, filename, min_freq=min_freq, max_freq=max_freq)


def simulate_S_pars(design, crop_box, filename, min_freq=6.5, max_freq=7.5, resolution_dx=1e3, resolution_dy=1e3):
    ### SIMULATION SECTION START ###
    ml_terminal = SonnetLab()
    from sonnetSim.cMD import CMD

    ml_terminal._send(CMD.SAY_HELLO)
    ml_terminal.clear()
    simBox = SimulationBox(
        crop_box.width(), crop_box.height(),
        crop_box.width() / resolution_dx,
        crop_box.height() / resolution_dy
    )

    ml_terminal.set_boxProps(simBox)
    from sonnetSim.pORT_TYPES import PORT_TYPES

    ports = [
        SonnetPort(prt, PORT_TYPES.BOX_WALL) for prt in design.sonnet_ports
    ]
    ml_terminal.set_ports(ports)
    ml_terminal.send_polygons(design.cell, design.layer_ph)
    ml_terminal.set_ABS_sweep(min_freq, max_freq)
    # print(f"simulating...{resonator_idx}")
    result_path = ml_terminal.start_simulation(wait=True)
    ml_terminal.release()

    all_params = design.get_geometry_parameters()

    # creating directory with simulation results
    results_dirpath = PROJECT_DIR / "circular_resonator"

    shutil.copy(
        result_path.decode("ascii"),
        str(results_dirpath / filename)
    )

    design.layout.write(str(results_dirpath / (filename[:-4] + '.gds')))
    design.save_logdata(results_dirpath / (filename[:-4] + '.xml'),
                        comment=f'S-parameters simulation results saved at {filename}.')

    ### RESULT SAVING SECTION END ###


### MAIN FUNCTION ###
if __name__ == "__main__":
    # my_design = CircleResonatorDesign("testScript")
    # my_design.draw()
    # my_design.show()

    for i in range(3, 11):
        simulate_full_S_pars(filename=f'cs_S12_{i}-{i+1}_GHz.csv', min_freq=float(i), max_freq=float(i+1))
