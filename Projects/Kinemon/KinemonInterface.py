from importlib import reload

import Projects.Dmon.Design
from Projects.Kinemon.CalculationIntefrace import CalculationInterface
from contextlib import redirect_stdout
import logging

# Variables dynamically assigned in script arguments
# mode = ""
# pts = ""
# par_d = ""
LOG_PATH = "./log.txt"
LOG_FORMAT = "%(asctime)s - %(level)s - %(message)s"
LOG_DATE = "%d.%m.%Y-%H:%M:%S"
logging.basicConfig(filename=LOG_PATH, level=logging.DEBUG)
# , format=LOG_FORMAT, datefmt=LOG_DATE
class KinemonInterface(CalculationInterface):
    def __init__(self, output_fname):
        self.output_fname = output_fname

    def process(self):
        global pts, mode, par_d
        logging.info('Process started')
        if mode == "Cqr":
            pts = int(pts)
            par_d = float(par_d)
            Projects.Dmon.Design.simulate_Cqr(resolution=(4e3, 4e3), mode="Cqr", pts=pts, par_d=par_d, output_fname=self.output_fname)
        else:
            logging.error('Incorrect mode')
        logging.info(f"pts={pts}, mode={mode}")


if __name__ == "__main__":
    reload(Projects.Dmon.Design)
    logging.info('File execution started')
    kinemon_interface = KinemonInterface('test_out.csv')
    kinemon_interface.process()