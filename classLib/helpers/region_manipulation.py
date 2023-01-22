import pya
from pya import DVector
from pya import Region
from pya import DCplxTrans


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
    primitive.make_trans(DCplxTrans(1, 0, False, DVector(-rotation_center)))
    primitive.make_trans(DCplxTrans(1, angle_deg, False, 0, 0))
    primitive.make_trans(DCplxTrans(1, 0, False, DVector(rotation_center)))
