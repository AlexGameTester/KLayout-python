from pya import DPoint, DVector, DPolygon, DCplxTrans

from Projects.Dmon.Design import RFSquidParams, DesignDmon
from classLib import ElementBase, CPW
from classLib.josJ import AsymSquid


def check_positive(val, full_name, short_name):
    if val < 0:
        raise ValueError(f"{full_name} must be positive but {short_name}={val:.1e} was given")


class MeanderParams:
    def __init__(self, dr: DPoint,
                 line_gap: float,
                 line_squares_n: int,
                 line_width_dx: float,
                 line_width_dy: float,
                 line_length: float,
                 add_dx_mid: float,
                 ):
        self.dr = dr

        check_positive(line_gap, "Line gap", "line_gap")
        self.line_gap = line_gap

        check_positive(line_squares_n, "Line squares number", "line_squares_n")
        self.line_squares_n = line_squares_n

        check_positive(line_width_dx, "Line width along Ox axis", "line_width_dx")
        self.line_width_dx = line_width_dx

        check_positive(line_width_dy, "Line width along Oy axis", "line_width_dy")
        self.line_width_dy = line_width_dy

        check_positive(line_length, "Total meander length", "line_length")
        self.line_length = line_length

        check_positive(add_dx_mid, "Additional Ox shift", "add_dx_mid")
        self.add_dx_mid = add_dx_mid


class KinIndMeander(ElementBase):
    # TODO: Возможно, нужно будет задавать trans_in, чтобы правильно распологать меандр относительно других элементов
    def __init__(self, meander_params: MeanderParams, trans_in=None, region_id="default"):
        self.meander_params = meander_params
        self.start = None
        self.end = None
        super().__init__(DPoint(0,0), trans_in=trans_in, region_id=region_id)

    def place(self, dest, layer_i=-1, region_id="default", merge=False):
        super().place(dest, layer_i, region_id, merge)

    def init_regions(self):
        self.dx = self.meander_params.dr.x
        self.dy = self.meander_params.dr.y 

        if abs(self.dy) < self.meander_params.line_gap:
            raise ValueError(f"Meander line gap must be less than its total height dy={abs(self.dy)} but line_gap={self.meander_params.line_gap} was given")

        self.n_periods = (self.dy - self.meander_params.line_gap) // \
                         (2 * self.meander_params.line_gap)
        self.n_periods = int(self.n_periods)

        self.dx_step = self.dx / (self.n_periods + 1)
        self.dy_step = self.dy / (2 * self.n_periods + 1)  # >= `line_gap`

        self.s = (self.meander_params.line_length - self.dy + self.dx - 2 * self.meander_params.add_dx_mid) / (
                    self.n_periods + 1) / 2

        poly = self.construct_poly()
        self.metal_region.insert(poly)

    def construct_poly(self):
        # Exterior points
        ext_points = []
        dy_hw = self.meander_params.line_width_dy / 2
        dx_hw = self.meander_params.line_width_dx / 2
        p1 = DPoint(0, -dy_hw)
        p2 = p1 + DVector(self.s / 2 + dx_hw, 0)
        p3 = p2 + DVector(0, self.dy_step + 2 * dy_hw)
        p4 = p3 + DVector(-self.s + self.dx_step, 0)
        ext_points += [p1, p2, p3, p4]
        for i in range(self.n_periods):
            p1 = ext_points[-1]
            p2 = p1 + DVector(0, self.dy_step - 2 * dy_hw)
            p3 = p2 + DVector(self.s, 0)
            p4 = p3 + DVector(0, self.dy_step + 2 * dy_hw)
            if i == self.n_periods - 1:
                p5 = p4 + DVector(-self.s / 2 + self.dx_step, 0)
            else:
                p5 = p4 + DVector(-self.s + self.dx_step, 0)

            ext_points += [p1, p2, p3, p4, p5]

        # Interior points
        int_points = []
        p1 = DPoint(0, dy_hw)
        p2 = p1 + DVector(self.s / 2 - dx_hw, 0)
        p3 = p2 + DVector(0, self.dy_step - 2 * dy_hw)
        p4 = p3 + DVector(-self.s + self.dx_step, 0)
        int_points += [p1, p2, p3, p4]

        for i in range(self.n_periods):
            p1 = int_points[-1]
            p2 = p1 + DVector(0, self.dy_step + 2 * dy_hw)
            p3 = p2 + DVector(self.s, 0)
            p4 = p3 + DVector(0, self.dy_step - 2 * dy_hw)

            # half of the width for the last point
            if i == self.n_periods - 1:
                p5 = p4 + DVector(-self.s / 2 + self.dx_step, 0)
            else:
                p5 = p4 + DVector(-self.s + self.dx_step, 0)

            int_points += [p1, p2, p3, p4, p5]

        # Last point correction
        int_points[-1] += DVector(dx_hw, 0)
        ext_points[-1] += DVector(-dx_hw, 0)

        int_points.reverse()
        points = ext_points + int_points

        return DPolygon(points)


class KinemonParams(RFSquidParams):
    def __init__(self, rf_sq_params: RFSquidParams, meander_params: MeanderParams, area_ratio=1 / 2,
                 MC_dy=None,
                 MC_dx=None,
                 KI_bridge_width=1e3,
                 KI_bridge_height=4e3,
                 KI_pad_y_offset=0.2e3,
                 KI_pad_width=3e3,
                 KI_ledge_y_offset=0.3e3):
        """
        @param rf_sq_params:
        @param meander_params: Object that contains parameters of kinetic inductance meander
        @param area_ratio:
        @param MC_dy:
        @param MC_dx:
        @param KI_bridge_width: Width of a thin bridge that connects kinetic inductance meander to its large contact pad
        @param KI_bridge_height: Height of a thin bridge that connects kinetic inductance meander to its large contact pad
        @param KI_pad_y_offset: Offset of end of kinetic inductance contact pad relatively to contact pads of josephson junctions
        @param KI_pad_width: Width of kinetic inductance contact pad
        """
        self.__dict__.update(rf_sq_params.__dict__)
        self.meander_params = meander_params
        # TODO: Implement variable area_ratio
        self.area_ratio = area_ratio

        if MC_dy is None:
            self.MC_dy = self.TC_dy
        else:
            self.MC_dy = MC_dy

        if MC_dx is None:
            self.MC_dx = self.TC_dx
        else:
            self.MC_dx = MC_dx

        self.KI_bridge_width = KI_bridge_width
        self.KI_bridge_height = KI_bridge_height
        self.KI_pad_y_offset = KI_pad_y_offset
        self.KI_pad_width = KI_pad_width
        self.KI_ledge_y_offset = KI_ledge_y_offset


class Kinemon(AsymSquid):
    def __init__(self, origin: DPoint, squid_params: KinemonParams, trans_in=None):
        self.center = origin

        self.r_curve = max(squid_params.line_width_dx,
                           squid_params.line_width_dy)
        self.TCBC_round_r = 500  # nm
        self.TCBC_round_n = 50
        # declare for proper code completion
        self.kin_ind_meander: KinIndMeander = None
        self.squid_params: KinemonParams = None
        super().__init__(origin=origin, params=squid_params,
                         trans_in=trans_in)

    def init_kin_ind(self):
        # Making meander
        assert len(self.BC_list) == 2
        self.BC_list.sort(key=lambda bc: bc.start.x)
        left_bc_object = self.BC_list[0]
        right_bc_object = self.BC_list[1]
        tc_object = self.TC

        tc = tc_object.end
        left_bc = left_bc_object.start
        right_bc = right_bc_object.start

        r1 = tc - left_bc
        r2 = tc - right_bc
        print(f'r1 = {r1}, r2 = {r2}')
        r0 = (r1 + r2) / 2 + DVector(0, -2 * self.squid_params.KI_bridge_height)
        # r0 = DVector(10e3, 10e3)

        start = tc - r0 + DVector(0, -self.squid_params.KI_bridge_height)
        end = start + r0

        self.squid_params.meander_params.dr = r0
        trans = DCplxTrans(start)
        self.kin_ind_meander = KinIndMeander(self.squid_params.meander_params, trans_in=trans, region_id="kinInd")
        self.kin_ind_meander.start = start
        self.kin_ind_meander.end = end
        self.primitives["kinIndMeander"] = self.kin_ind_meander

        # Making meander contacts
        # top kinetic inductance bridge
        self.TKIB = CPW(start=self.kin_ind_meander.end + DVector(0, - self.squid_params.meander_params.line_width_dy / 2),
                        end=self.kin_ind_meander.end + DVector(0, self.squid_params.KI_bridge_height),
                        width=self.squid_params.KI_bridge_width, gap=0,
                        region_id="kinInd")
        self.primitives["TKIB"] = self.TKIB

        self.BKIB = CPW(start=self.kin_ind_meander.start + DVector(0, + self.squid_params.meander_params.line_width_dy / 2),
                        end=self.kin_ind_meander.start + DVector(0, -self.squid_params.KI_bridge_height),
                        width=self.squid_params.KI_bridge_width, gap=0,
                        region_id="kinInd")
        self.primitives["BKIB"] = self.BKIB

        # top kinetic inductance pad
        self.TKIP = CPW(start=self.TKIB.end,
                        end=self.TC.start + DVector(0, -self.squid_params.KI_pad_y_offset),
                        width=self.squid_params.KI_pad_width,
                        gap=0,
                        region_id="kinInd")
        self.TKIP.metal_region.round_corners(0, self.TCBC_round_r, self.TCBC_round_n)
        self.primitives["TKIP"] = self.TKIP

        bkip_end = DVector((self.BC_list[0].start.x + self.BC_list[1].start.x) / 2, self.BC_list[0].end.y + self.squid_params.KI_pad_y_offset)
        self.BKIP = CPW(start=self.BKIB.end,
                        end=bkip_end,
                        width=self.squid_params.KI_pad_width,
                        gap=0,
                        region_id="kinInd")
        self.BKIP.metal_region.round_corners(0, self.TCBC_round_r, self.TCBC_round_n)
        self.primitives["BKIP"] = self.BKIP
        self.BC_list.append(self.BKIP)

        # # botttom kinetic inductance ledge
        # self.BKIL = CPW(start=self.BKIB.end,
        #                 end=bkip_end + DVector(0, -self.squid_params.KI_pad_y_offset + self.squid_params.KI_ledge_y_offset),
        #                 width=self.squid_params.BC_dx[0],
        #                 gap=0,
        #                 region_id="photo")
        # self.primitives["BKIL"] = self.BKIL


    # def init_kin_ind_contact(self):
    #     print("Kinemon init_kin_ind_contact")
    #     # Creating a ledge in isolator across conductance cross
    #     # TODO: this

        # contact_start = self.kin_ind_meander.start + DVector(0, -self.squid_params.MC_dy)
        # contact_end = self.kin_ind_meander.start
        # self.BC_KI = CPW(start=contact_start, end=contact_end, width=self.squid_params.MC_dx / 2, gap=0)
        # self.MC = CPW(start=contact_start, end=contact_end + DVector(0, self.squid_params.MC_dy / 2), width=self.squid_params.MC_dx, gap=0)
        # # self.MC = CPW(start=contact_start, end=contact_end, width=self.squid_params.MC_dx, gap=0)
        # self.primitives["MC"] = self.MC
        # self.primitives["BC_KI"] = self.BC_KI
        #
        # meander_finish = self.kin_ind_meander.start + self.kin_ind_meander.meander_params.dr
        # self.TC_KI = CPW(start=meander_finish, end=meander_finish + DVector(0, self.squid_params.MC_dy / 2),
        #                  width=self.squid_params.MC_dx / 2, gap=0)


    def init_primitives(self):
        super().init_primitives()

        self.init_kin_ind()

class DesignKinemon(DesignDmon):
    def __init__(self, cell_name="testScript"):
        super().__init__(cell_name)

    def draw(self):
        self.draw_chip()



    def draw_kin_ind(self):
        pass


if __name__ == "__main__":
    kmon = DesignKinemon()
