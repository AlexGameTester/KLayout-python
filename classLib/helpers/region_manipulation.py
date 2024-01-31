# import built-ins
import logging

# import good 3rd party

# import project specific 3rd party
import pya
from pya import DVector
from pya import Region
from pya import ICplxTrans, DCplxTrans
logging.debug("classlib.helpers.region_manipulation loaded")

def extended_region(reg, extension=0):
    """
    extends region in outer direction by `extension` value.
    extension is orthogonal to region`s polygons edges.
    Parameters
    ----------
    reg : Region
        pya.Region instance
    extension : float
        extension in nm
    Returns
    -------
    Region
        extended version of region
    """
    tmp_reg = Region()
    ep = pya.EdgeProcessor()
    for poly in reg.each():
        tmp_reg.insert(
            ep.simple_merge_p2p(
                [
                    poly.sized(
                        extension,
                        extension,
                        2
                    )
                ],
                False,
                False,
                1
            )
        )
    return tmp_reg


def rotate_around(primitive, rotation_center, angle_deg):
    rotation_center = rotation_center.dup()
    trans = DCplxTrans(1, 0, False, DVector(rotation_center)) * \
            DCplxTrans(1, angle_deg, False, 0, 0) * \
            DCplxTrans(1, 0, False, -DVector(rotation_center))
    if isinstance(primitive, Region):
        primitive.transform(ICplxTrans(trans))
    elif hasattr(primitive, "make_trans"):
        primitive.make_trans(trans)
    else:
        logging.error(
            "unknown datastructure passed to classlib.helpers.rotate_around"
            f"object type is f{type(primitive)}\n"
        )
