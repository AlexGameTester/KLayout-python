import numpy as np
from pya import DPoint, DVector, DPolygon

from Projects.Dmon.Design import RFSquidParams, DPathCPWStraight, DesignDmon
from classLib import CPWParameters, ChipDesign, ElementBase
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
        super().__init__(DPoint(0,0), trans_in=trans_in, region_id=region_id)

    def init_regions(self):
        self.dx = self.meander_params.dr.x
        self.dy = self.meander_params.dr.y 

        if self.dy < self.meander_params.line_gap:
            print(
                "Kinetic inductance meander drawing error:"
                f"impossible to draw line with meander step "
                f"{self.meander_params.line_gap:.0f} nm"
                f"total dy of meander is too small: dy = {self.dy:.f} nm\n"
            )
            return
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
        p2 = p1 + DVector(self.s + dx_hw, 0)
        p3 = p2 + DVector(0, self.dy_step + 2 * dy_hw)
        p4 = p3 + DVector(-self.s + self.dx_step, 0)
        ext_points += [p1, p2, p3, p4]
        for i in range(self.n_periods):
            p1 = ext_points[-1]
            p2 = p1 + DVector(0, self.dy_step - 2 * dy_hw)
            p3 = p2 + DVector(self.s, 0)
            p4 = p3 + DVector(0, self.dy_step + 2 * dy_hw)
            p5 = p4 + DVector(-self.s + self.dx_step, 0)
            ext_points += [p1, p2, p3, p4, p5]

        # Interior points
        int_points = []
        p1 = DPoint(0, dy_hw)
        p2 = p1 + DVector(self.s - dx_hw, 0)
        p3 = p2 + DVector(0, self.dy_step - 2 * dy_hw)
        p4 = p3 + DVector(-self.s + self.dx_step, 0)
        int_points += [p1, p2, p3, p4]

        for i in range(self.n_periods):
            p1 = int_points[-1]
            p2 = p1 + DVector(0, self.dy_step + 2 * dy_hw)
            p3 = p2 + DVector(self.s, 0)
            p4 = p3 + DVector(0, self.dy_step - 2 * dy_hw)
            p5 = p4 + DVector(-self.s + self.dx_step, 0)
            int_points += [p1, p2, p3, p4, p5]

        # Last point correction
        int_points[-1] += DVector(dx_hw, 0)
        ext_points[-1] += DVector(-dx_hw, 0)

        int_points.reverse()
        points = ext_points + int_points

        return DPolygon(points)


class KinemonParams(RFSquidParams):
    def __init__(self, rf_sq_params: RFSquidParams):
        self.__dict__.update(rf_sq_params.__dict__)


class Kinemon(AsymSquid):
    def __init__(self, origin: DPoint, squid_params: KinemonParams, trans_in=None):
        self.center = origin

        self.r_curve = max(squid_params.line_width_dx,
                           squid_params.line_width_dy)
        self.TCBC_round_r = 500  # nm
        # declare for proper code completion
        self.squid_params: RFSquidParams = None
        super().__init__(origin=origin, params=squid_params,
                         trans_in=trans_in)

    def init_primitives(self):
        super().init_primitives()


class DesignKinemon(DesignDmon):
    def __init__(self, cell_name="testScript"):
        super().__init__(cell_name)

    def draw(self):
        self.draw_chip()

    def draw_kin_ind(self):
        pass


if __name__ == "__main__":
    kmon = DesignKinemon()
