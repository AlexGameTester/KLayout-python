from importlib import reload

from classLib.helpers import polygon_splitting
reload(polygon_splitting)
split_polygons = polygon_splitting.split_polygons

from classLib.helpers import pinning_grid
reload(pinning_grid)
fill_holes = pinning_grid.fill_holes

from classLib.helpers import region_manipulation
reload(region_manipulation)
extended_region = region_manipulation.extended_region
rotate_around = region_manipulation.rotate_around

from classLib.helpers import capacitance_sim
reload(capacitance_sim)
simulate_cij = capacitance_sim.simulate_cij
save_sim_results = capacitance_sim.save_sim_results


class FABRICATION:
    """
        Metal polygons edges are overetched by this value expressen in nm.
        Overetching results in broadening of gaps between polygons.
        In other words, every polygon edge is shifted along direction
    perpendicular to the edge itself from empty space
    to polygon's body by FABRICATION.OVERETCHING distance in nm.
        To account for overetching polygons has to be constructed
    in width way that results in software design with polygons "widened" by
    FABRICATIO.OVERETCHING value. For e.g. witdth of the coplanar
    waveguide central conductor has to be "widened" by 2*FABRICATION.OVERETCHING
    while preseving symmetry along center of the wavegiude.

        Correponding adjustments have to be made to every design element
    that undergoues overetching during fabrication.
    In addition, different areas of the sample can undergo different
    overetching, depending on the design and fabrication process.
    """
    # 0.0 - for development
    # 0.8e3 - estimation for fabrication by Bolgar photolytography etching
    # recipe
    OVERETCHING = 0.0e3