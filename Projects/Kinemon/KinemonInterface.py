from importlib import reload

import Projects.Dmon.Design
from Projects.Kinemon.CalculationIntefrace import CalculationInterface
from contextlib import redirect_stdout
import logging

# Variables dynamically assigned in script arguments
# mode = ""
# pts = ""
# par_d = ""
# output_file = ""
LOG_PATH = "./log.txt"
LOG_FORMAT = "%(asctime)s - %(level)s - %(message)s"
LOG_DATE = "%d.%m.%Y-%H:%M:%S"
logging.basicConfig(filename=LOG_PATH, level=logging.DEBUG)
# , format=LOG_FORMAT, datefmt=LOG_DATE
class KinemonInterface(CalculationInterface):
    def __init__(self, output_fname):
        self.output_fname = output_fname

    def process(self):
        global pts, mode, par_d, output_file
        design = Projects.Dmon.Design.DesignDmon('testScript')
        Pe


if __name__ == "__main__":
    reload(Projects.Dmon.Design)
    logging.info('File execution started')
    kinemon_interface = KinemonInterface('test_out.csv')
    kinemon_interface.process()