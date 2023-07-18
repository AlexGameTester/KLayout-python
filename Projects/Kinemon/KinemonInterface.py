from Projects.Dmon import Design as dmon
from Projects.Kinemon.CalculationIntefrace import CalculationInterface
from contextlib import redirect_stdout
import logging

# Variables dynamically assigned in script arguments
# mode = ""
# pts = ""
LOG_PATH = "./log.txt"
logging.basicConfig(filename=LOG_PATH, level=logging.DEBUG)


class KinemonInterface(CalculationInterface):
    def __init__(self):
        # self.dmon = dmon.DesignDmon('testScript')
        # self.stdout_f = open('out.txt', 'w')
        # redirect_stdout(self.stdout_f)
        pass

    def process(self):
        global pts, mode
        logging.info('Process started')
        if mode == "Cqr":
            pts = int(pts)
            dmon.simulate_Cqr(resolution=(4e3, 4e3), mode="Cqr", pts=pts, par_d=10e3)
        else:
            logging.error('Incorrect mode')
        logging.info(f"pts={pts}, mode={mode}")


if __name__ == "__main__":
    logging.info('File execution started')
    kinemon_interface = KinemonInterface()
    kinemon_interface.process()