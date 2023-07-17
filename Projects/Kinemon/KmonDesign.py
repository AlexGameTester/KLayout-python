from classLib import ComplexBase
from classLib.josJ import AsymSquidParams, AsymSquid
from Projects.Dmon.Design import RFSquidParams, DPathCPWStraight
from pya import DPoint


class MeanderParams:
    def __init__(self, dr: DPoint, line_gap: float, line_squares_n: int, line_width_dy: float):
        self.dr = dr

        if line_gap < 0:
            raise ValueError(f"Line gap must be positive but line_gap={line_gap:.1e} was given")
        self.line_gap = line_gap

        if line_squares_n < 0:
            raise ValueError(f"Line squares number must be positive but line_squares_n={line_squares_n} was given")
        self.line_squares_n = line_squares_n

        if line_width_dy < 0:
            raise ValueError(f"Line width must be positive but line_width_dy={line_width_dy} was given")
        self.line_width_dy = line_width_dy




class KinIndMeander(DPathCPWStraight):
    # TODO: Возможно, нужно будет задавать trans_in, чтобы правильно распологать меандр относительно других элементов
    def __init__(self, meander_params: MeanderParams, trans_in=None, region_id="default"):
        self.meander_params = meander_params
        points, cpw_pars = self._initialize_dpath()
        super().__init__(points, cpw_pars, trans_in, region_id)

    def _initialize_dpath(self):
        self.dx = self.meander_params.dr.dx
        self.dy = self.meander_params.dr.dy

        if self.dy < self.meander_params.line_gap:
            print(
                "Kinetic inductivity meander drawing error:"
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

        y_squares_n = self.dy/self.meander_params.line_width_dy
        if y_squares_n > self.meander_params.line_squares_n:
            print("Error, to little squares for fixed inductor height.")



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



class DesignKinemon:
    pass

if __name__ == "__main__":
    kmon = DesignKinemon()

