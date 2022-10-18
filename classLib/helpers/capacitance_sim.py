import os
import csv

from sonnetSim import SonnetLab, SonnetPort

from pya import Region

def simulate_cij(
        design, reg, subreg1, subreg2=None, resolution=(5e3,5e3)
):
    '''
    Measures capacity matrix between polygons represented by `Region()`
    class in 1to1 correspondence
    1 obj - scalar, assuming all other polygons grounded.
    2 obj - Cij Kirchgoff matrix, assuming all other polygons grounded.
    Only <=2 polygons are currently supported.

    Parameters
    ----------
    reg : Region
        simulation environment region
    subreg1 : Region
        region of 1st metal struct
    subreg2 : DPoint
        region of 2nd metal struct

    resolution : tuple(float,float)
        dx and dy mesh step for Sonnet

    Returns
    -------
    None
    '''
    resolution_dx, resolution_dy = resolution

    ''' CROP NECESSARY DESIGN PART '''
    bbox1, bbox2 = subreg1.bbox(), subreg2.bbox()
    crop_box = Region((subreg1 + subreg2).bbox())
    crop_box = crop_box.sized(
        3*(bbox1.width() + bbox2.width()),
        3*(bbox1.height() + bbox2.height())
    )
    design.crop(crop_box)

    ''' PLACE SONNET PORTS '''
    from itertools import product
    edgeCenter_cr1_best, edgeCenter_cr2_best = None, None
    max_distance = 0
    edge_centers_it = product(
        subreg1.edges().centers(0, 0).each(),
        subreg2.edges().centers(0, 0).each()
    )
    edge_centers_it = map(
        lambda edge_tuple: (edge_tuple[0].p1, edge_tuple[1].p1),
        edge_centers_it
    )
    for edgeCenter_cr1, edgeCenter_cr2 in edge_centers_it:
        centers_d = edgeCenter_cr1.distance(edgeCenter_cr2)
        if centers_d > max_distance:
            edgeCenter_cr1_best, edgeCenter_cr2_best = \
                edgeCenter_cr1, edgeCenter_cr2
            max_distance = centers_d
        else:
            continue

    design.sonnet_ports.append(edgeCenter_cr1_best)
    design.sonnet_ports.append(edgeCenter_cr2_best)

    # x_distance_dx_list = [0]
    for dl in dl_list:
        ''' DRAWING SECTION START '''
        design = Design8QStair("testScript")
        q1 = design.qubits[q1_idx]
        q2 = design.qubits[q2_idx]
        design.N_coils = [1] * design.NQUBITS

        q1.qubit_params.qubit_cap_params.disk_r += dl
        print("disk_r = ", q1.qubit_params.qubit_cap_params.disk_r)

        design.draw()
        design.show()
        design.layout.write(
            os.path.join(PROJECT_DIR, f"Cqq_{q1_idx}_{q2_idx}_"
                                      f"{dl:.3f}_.gds")
        )


        from itertools import product

        edgeCenter_cr1_best, edgeCenter_cr2_best = None, None
        max_distance = 0
        edge_centers_it = product(
            cross1.metal_region.edges().centers(0, 0).each(),
            cross2.metal_region.edges().centers(0, 0).each()
        )
        edge_centers_it = map(
            lambda edge_tuple: (edge_tuple[0].p1, edge_tuple[1].p1),
            edge_centers_it
        )
        for edgeCenter_cr1, edgeCenter_cr2 in edge_centers_it:
            centers_d = edgeCenter_cr1.distance(edgeCenter_cr2)
            if centers_d > max_distance:
                edgeCenter_cr1_best, edgeCenter_cr2_best = \
                    edgeCenter_cr1, edgeCenter_cr2
                max_distance = centers_d
            else:
                continue

        design.sonnet_ports.append(edgeCenter_cr1_best)
        design.sonnet_ports.append(edgeCenter_cr2_best)

        dr = DPoint(0, 0) - crop_box.p1
        design.transform_region(design.region_ph, DTrans(dr.x, dr.y),
                                trans_ports=True)

        design.show()
        design.lv.zoom_fit()
        '''DRAWING SECTION END'''

        '''SIMULATION SECTION START'''
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
        ml_terminal.set_linspace_sweep(0.01, 0.01, 1)
        print("simulating...")
        result_path = ml_terminal.start_simulation(wait=True)
        ml_terminal.release()

        ### SIMULATION SECTION END ###

        ### CALCULATE CAPACITANCE SECTION START ###
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
            import math

            delta = (1 + s[0][0]) * (1 + s[1][1]) - s[0][1] * s[1][0]
            y11 = 1 / R * ((1 - s[0][0]) * (1 + s[1][1]) + s[0][1] * s[1][
                0]) / delta
            y22 = 1 / R * ((1 - s[1][1]) * (1 + s[0][0]) + s[0][1] * s[1][
                0]) / delta
            C1 = -1e15 / (2 * math.pi * freq0 * 1e9 * (1 / y11).imag)
            C2 = -1e15 / (2 * math.pi * freq0 * 1e9 * (1 / y22).imag)
            # formula taken from https://en.wikipedia.org/wiki/Admittance_parameters#Two_port
            y21 = -2 * s[1][0] / delta * 1 / R
            C12 = 1e15 / (2 * math.pi * freq0 * 1e9 * (1 / y21).imag)

        print("C_12 = ", C12)
        print("C1 = ", C1)
        print("C2 = ", C2)
        print()
        '''CALCULATE CAPACITANCE SECTION END'''

        '''SAVING REUSLTS SECTION START'''
        geometry_params = design.get_geometry_parameters()
        output_filepath = os.path.join(PROJECT_DIR, "Xmon_Cqq_results.csv")
        if os.path.exists(output_filepath):
            # append data to file
            with open(output_filepath, "a", newline='') as csv_file:
                writer = csv.writer(csv_file)
                writer.writerow(
                    [q1_idx, q2_idx, *list(geometry_params.values()),
                     design.xmon_x_distance / 1e3,
                     C1, C12]
                )
        else:
            # create file, add header, append data
            with open(output_filepath, "w", newline='') as csv_file:
                writer = csv.writer(csv_file)
                # create header of the file
                writer.writerow(
                    ["q1_idx", "q2_idx", *list(geometry_params.keys()),
                     "xmon_x_distance, um",
                     "C1, fF", "C12, fF"])
                writer.writerow(
                    [q1_idx, q2_idx, *list(geometry_params.values()),
                     design.xmon_x_distance / 1e3,
                     C1, C12]
                )
        '''SAVING REUSLTS SECTION END'''