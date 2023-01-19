# import built-ins
import os
import csv
from typing import List, Tuple

# import good 3rd party
import numpy as np

# import project specific 3rd party
from pya import Region, DPoint, DTrans, Box

# import project lib
from sonnetSim import SonnetLab, SonnetPort, SimulationBox
from classLib.chipDesign import ChipDesign


# TODO: `reg` and `layer` will be no longer needed both after introducing
#  region_ids global dictionary
def simulate_cij(
    design, layer, subregs, env_reg=None, resolution=(5e3, 5e3),
    print_values=True
):
    '''
    Measures capacity matrix between polygons supplied in `subregs`
    in metal environment supplied by `reg`.
    1 obj - scalar, assuming all other polygons grounded.
    2 obj - Cij Kirchgoff matrix, assuming all other polygons grounded.
    Only <=2 polygons are currently supported.

    Parameters
    ----------
    design : ChipDesign
        design instance
    subregs : List[Region]
        list of different metal polygons to calculate capacitance between.
    env_reg : Region
        simulation environment region.
        If not specified will be derived from `subregs` parameter
        by enlarging it in 2 to 3 times
    resolution : tuple(float,float)
        dx and dy mesh step for Sonnet
    print_values : bool
        Whether to print capacitance values in KLayout console

    Returns
    -------
    Tuple[float,float,float]
        C1, C12, C2
    '''

    resolution_dx, resolution_dy = resolution
    if env_reg is None:
        crop_box = sum(subregs, Region()).bbox()
        # enlarge box around its center
        dr = crop_box.p2 - crop_box.p1
        center = crop_box.center()
        crop_box = Box(center - 3 / 2 * dr, center + 3 / 2 * dr)
    else:
        crop_box = env_reg.bbox()

    ''' SONNET PORTS POSITIONS SECTION START '''
    n_terminals = len(subregs)
    if n_terminals == 1:
        # find smallest width terminal
        edge_center_best = None
        min_width = 1e9
        edge_centers_it = subregs[0].edges().each()
        for edge in edge_centers_it:
            if edge.length() < min_width:
                # AUTOGROUNDED PORTS CANNOT BE DIAGONAL (TODO: CHANGE LATER)
                if any(
                        [
                            (edge.d().y == 0),
                            (edge.d().x == 0)
                        ]
                ):
                    min_width = edge.length()
                    edge_center_best = (edge.p1 + edge.p2) / 2
        design.sonnet_ports.append(edge_center_best)
    elif n_terminals == 2:
        # TODO: not tested
        from itertools import product
        edge1_best, edge2_best = None, None
        max_distance = 0
        edge_centers_it = product(
            subregs[0].edges().centers(0, 0).each(),
            subregs[1].edges().centers(0, 0).each()
        )
        edge_centers_it = map(
            lambda edge_tuple: (edge_tuple[0].p1, edge_tuple[1].p1),
            edge_centers_it
        )
        for edge1, edge2 in edge_centers_it:
            centers_d = edge1.distance(edge2)
            if all(
                    [
                        centers_d > max_distance,
                        any([(edge1.d().y == 0), (edge1.d().x == 0)]),
                        any([(edge1.d().y == 0), (edge1.d().x == 0)])
                    ]
            ):
                edge1_best, edge2_best = edge1, edge2
                max_distance = centers_d
            else:
                continue

        design.sonnet_ports.append(edge1_best)
        design.sonnet_ports.append(edge2_best)
    # [print(port.x, port.y) for port in design.sonnet_ports]
    ''' SONNET PORTS POSITIONS SECTION END '''

    ''' CROP + SNAP TO ORIGIN SECTION START '''
    dr = DPoint(0, 0) - crop_box.p1
    print(crop_box.p1.x, crop_box.p1.y)
    print(crop_box.p2.x, crop_box.p2.y)
    design.crop(crop_box, design.region_ph)
    print("sonnet ports positions BEFORE transform:")
    for sp in design.sonnet_ports:
        print(sp.x, sp.y)
    design.transform_region(
        design.region_ph, DTrans(dr.x, dr.y),
        trans_ports=True
        )
    design.layout.clear_layer(layer)
    design.show()
    ''' CROP + SNAP TO ORIGIN SECTION END '''

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
    print(simBox.y, simBox.x, simBox.x_n, simBox.x_n)
    ml_terminal.set_boxProps(simBox)
    # print("sending cell and layer")
    from sonnetSim.pORT_TYPES import PORT_TYPES

    ports = []
    for sonnet_port in design.sonnet_ports:
        ports.append(
            SonnetPort(
                sonnet_port,  # `nm` to `um`
                PORT_TYPES.AUTOGROUNDED
            )
        )

    print("sonnet ports positions after transform:")
    for sp in ports:
        print(sp.point)
    ml_terminal.set_ports(ports)

    ml_terminal.send_polygons(design.cell, layer_i=layer)
    ml_terminal.set_linspace_sweep(0.01, 0.01, 1)
    print("simulating...")
    result_path = ml_terminal.start_simulation(wait=True)
    ml_terminal.release()

    """ SIMULATION SECTION END """

    """ CALCULATE CAPACITANCE SECTION START """
    C1 = None
    C12 = None
    C2 = None
    n_terminals = len(design.sonnet_ports)
    with open(result_path.decode("ascii"), "r") as csv_file:
        data_rows = list(csv.reader(csv_file))
        ports_imps_row = data_rows[6]
        R = float(ports_imps_row[0].split(' ')[1])
        data_row = data_rows[8]
        freq0 = float(data_row[0])

        # lowest frequency S-matrix will be stored here
        s = np.zeros((n_terminals, n_terminals), dtype=complex)
        # print(data_row)
        for i in range(0, n_terminals):
            for j in range(0, n_terminals):
                s[i, j] = complex(
                    float(data_row[1 + 2 * (i * 2 + j)]),
                    float(
                        data_row[
                            1 + 2 * (i * 2 + j) + 1]
                        )
                    )
        import math
        if n_terminals == 1:
            y11 = 1 / R * (1 - s[0, 0]) / (1 + s[0, 0])
            C1 = -1e15 / (2 * math.pi * freq0 * 1e9 * (1 / y11).imag)
            if print_values:
                print("C1 = ", C1)
        elif n_terminals == 2:
            delta = (1 + s[0, 0]) * (1 + s[1, 1]) - s[0, 1] * s[1, 0]
            y11 = 1 / R * (
                    (1 - s[0, 0]) * (1 + s[1, 1]) +
                    s[0, 1] * s[1, 0]
            ) / delta
            y22 = 1 / R * (
                    (1 - s[1, 1]) * (1 + s[0, 0]) +
                    s[0, 1] * s[1, 0]
            ) / delta
            C1 = -1e15 / (2 * math.pi * freq0 * 1e9 * (1 / y11).imag)
            C2 = -1e15 / (2 * math.pi * freq0 * 1e9 * (1 / y22).imag)
            # formula taken from https://en.wikipedia.org/wiki/Admittance_parameters#Two_port
            y21 = -2 * s[1][0] / delta * 1 / R
            C12 = 1e15 / (2 * math.pi * freq0 * 1e9 * (1 / y21).imag)
            if print_values:
                print("C1 = ", C1)
                print("C_12 = ", C12)
                print("C2 = ", C2)
    '''CALCULATE CAPACITANCE SECTION END'''

    return C1, C12, C2


def save_sim_results(output_filepath, design, additional_pars):
    """
    Saves simulation results. Associates single filename with single
    simulation type (no matter how many points acquired).

    Parameters
    ----------
    output_filepath : str
        path to the output file
    design : ChipDesign
        design instance with simulated parameters
    additional_pars : Dict[str,Any]
        dedicated parameters exctracted from simulation

    Returns
    -------
    None
    """
    design_pars = design.get_geometry_parameters()
    if os.path.exists(output_filepath):
        # append data to file
        with open(output_filepath, "a", newline='') as csv_file:
            writer = csv.writer(csv_file)
            writer.writerow(
                list(additional_pars.values()) +
                list(design_pars.values())
            )
    else:
        # create file, add header, append data
        with open(output_filepath, "w", newline='') as csv_file:
            writer = csv.writer(csv_file)
            # create header of the file
            writer.writerow(
                list(additional_pars.keys()) +
                list(design_pars.keys())
            )
            writer.writerow(
                list(additional_pars.values()) +
                list(design_pars.values())
            )
