#from Projects.SFS_transmon.MyDesign import My_Design

import pya
from math import sqrt, cos, sin, tan, atan2, pi, copysign
from pya import Point, DPoint, Vector, DVector, DSimplePolygon, SimplePolygon, DPolygon, Polygon, Region
from pya import Trans, DTrans, CplxTrans, DCplxTrans, ICplxTrans

from importlib import reload
from ClassLib import *
reload(ClassLib)
from ClassLib import *

class Test_Squid(Complex_Base):
    """ @brief:     class represents a rectangular capacitor with a dc-SQUID between its plates
        @params:    DPoint origin - position of the center of a structure
                    params{} - a dictionary with geometric parameters of a capacitor
                    squid_params - a list with dc-SQUID parameters
                    Trans trans_in - initial transformation (None by default)
    """
    def __init__(self, origin, params, squid_params, trans_in=None):
        self.width = params['width']
        self.height = params['height']
        self.innergap = params['innergap']
        self.outergap = params['outergap']
        self.squid_params = squid_params
        super().__init__(origin, trans_in)

    def init_primitives(self):
        origin = DPoint(0, 0)
        self.primitives['empty_rect'] = Rectangle(origin - DPoint(self.width/2 + self.outergap, self.height + self.innergap / 2 + self.outergap),
                                    self.width + 2 * self.outergap,
                                    2 * self.height + 2 * self.outergap + self.innergap,
                                    inverse=True)  
        self.primitives['top_rect'] = Rectangle(origin + DPoint(-self.width/2, self.innergap/2),
                                    self.width,
                                    self.height)
        self.primitives['bottom_rect'] = Rectangle(origin - DPoint(self.width/2, self.height + self.innergap / 2),
                                    self.width,
                                    self.height)  
        self.squid = Squid(origin, self.squid_params)
        self.primitives['qubit'] = self.squid

    def place(self, dest, layer_ph=-1, layer_el=-1):
        for prim_name in list(self.primitives.keys())[:-1]:
            self.primitives[prim_name].place(dest, layer_ph)
        self.squid.place(dest, layer_el)  

class My_Design(Chip_Design):

    origin = DPoint(0, 0)
    Z = CPWParameters(20e3, 10e3) # normal CPW
    Z_narrow = CPWParameters(10e3,7e3) # narrow CPW
    cpw_curve = 200e3 # Curvature of CPW angles
    chip = None

    # Call other methods drawing parts of the design from here
    def draw(self):
        self.draw_chip()
        self.draw_mark()
        self.draw_single_photon_source()
        self.draw_mixing_qubit()
        self.draw_resonators()
        self.draw_test_squids()
    
    def draw_chip(self):
        Z_params = [self.Z_narrow] + [self.Z]*7
        self.chip = Chip5x10_with_contactPads(self.origin, Z_params)
        self.chip.place(self.cell, self.layer_ph)  
    
    def draw_mark(self):
        # Placing the mark
        mark = Mark2(DPoint(500e3, self.chip.chip_y - 500e3))
        mark.place(self.cell, self.layer_ph)   
    
    def draw_single_photon_source(self):
        # Ports to which a single photon source is connected
        p1 = self.chip.connections[0]
        p2 = self.chip.connections[5]

        # Single photon source and dc-SQUID parameters
        pars = self.get_sps_params()
        pars_squid = self.get_dc_squid_params()
        
        # Drawing a qubit
        p = p1 + DPoint((p2-p1).x, (p2-p1).y/2) # position of a qubit
        sfs = SFS_Csh_emb(p, pars, pars_squid)
        sfs.place(self.cell, self.layer_ph, self.layer_el)
        
        # Input line
        sfs_line_in = CPW_RL_Path(p1, "LRL", pars['Z1'], self.cpw_curve, [(p2 - p1).x, -(p2 - p1).y / 2 - pars['r_out']], -pi / 2)
        sfs_line_in.place(self.cell, self.layer_ph)
        
        # Output line
        sfs_line_out = CPW(pars['Z2'].width, pars['Z2'].gap, sfs.connections[1], p2)
        sfs_line_out.place(self.cell, self.layer_ph)
    
    def draw_mixing_qubit(self):
        # Ports to which a maxing qubit is connected
        p1 = self.chip.connections[3]
        p2 = self.chip.connections[4]

        # Mixing qubit and dc-SQUID parameters
        pars = self.get_mixing_qubit_params()
        pars_squid = self.get_dc_squid_params()
    
        # Drawing the probe line
        probe_line = CPW_RL_Path(p1, "LRLRL", self.Z, self.cpw_curve, [self.chip.chip_y/10, (p2-p1).x, self.chip.chip_y/10], [pi/2, pi/2], trans_in=Trans.R270)
        probe_line.place(self.cell, self.layer_ph)
        
        # Drawing the qubit near the probe line
        to_line = 40e3 # Distance between the qubit and the probe line
        p = DPoint((p2.x + p1.x)/2, p1.y - abs(self.chip.chip_y/10) - to_line - pars['r_out']) # Position of the qubit
        sfs = SFS_Csh_emb(p, pars, pars_squid, squid_pos=-1, trans_in=Trans.R180)
        sfs.place(self.cell, self.layer_ph, self.layer_el)

        # Drawing a flux bias line
        p3 = self.chip.connections[2] # port connected to the flux bias line
        conn_len = 50e3 # length of a connection between two parts of the line
        dx = (p-p3).x - pars['r_out'] - 40e3 - conn_len
        dy = p3.y - p.y + pars['r_out']
        flux_line1 = CPW_RL_Path(p3, 'LRL', self.Z, self.cpw_curve, [dy, dx], pi/2, trans_in=Trans.R270)
        flux_line1.place(self.cell, self.layer_ph)
        Z_end = CPWParameters(5e3, 5e3) # parameters of a CPW at the end of the line
        p_start = p3 + DPoint(dx, -dy)
        p_end = p3 + DPoint(dx + conn_len, -dy)
        flux_conn = CPW2CPW(self.Z, Z_end, p_start, p_end)
        flux_conn.place(self.cell, self.layer_ph)
        flux_line2 = CPW_RL_Path(p_end, 'RL', Z_end, 50e3, [pars['r_out'] + 40e3], pi/2, trans_in=None)
        flux_line2.place(self.cell, self.layer_ph)
    
    def draw_resonators(self):
        # drawing the line
        p1 = self.chip.connections[1]
        p2 = self.chip.connections[7]
        probe_line = CPW_RL_Path(p1, "LRL", self.Z, self.cpw_curve, [(p1-p2).x, (p1-p2).y], pi / 2, trans_in = Trans.R180)
        probe_line.place(self.cell, self.layer_ph)

        # drawing resonators
        to_line = 30e3
        self.draw_one_resonator(DPoint(p2.x + to_line + 2 * (self.Z.width/2 + self.Z.gap), p1.y * 1/3 + p2.y * 2/3), freq=5, trans_in=Trans.R90)
        self.draw_one_resonator(DPoint(p1.x * 2/5 + p2.x * 3/5 , p1.y - to_line - 2 * (self.Z.width/2 + self.Z.gap)), freq=5.1)
        self.draw_one_resonator(DPoint(p1.x * 4/5 + p2.x * 1/5, p1.y - to_line - 2 * (self.Z.width/2 + self.Z.gap)), freq=4.9)
    
    def draw_one_resonator(self, pos, freq, trans_in=None):
        turn_radius = 50e3
        eps = 11.45
        wavelength_fraction = 1/4
        coupling_length = 300e3
        meander_periods = 4
        neck_length = 900e3
        worm = CPWResonator2(pos, self.Z, turn_radius, freq, eps, wavelength_fraction, coupling_length, meander_periods, neck_length, trans_in)
        worm.place(self.cell, self.layer_ph)

    def draw_test_squids(self):
        pars_probe = {'width': 300e3, 'height': 200e3, 'innergap': 30e3, 'outergap': 30e3}
        pars_squid = self.get_dc_squid_params()
        Test_Squid(DPoint(1e6, 1e6), pars_probe, pars_squid).place(self.cell, self.layer_ph, self.layer_el)
        Test_Squid(DPoint(self.chip.chip_x - 1e6, self.chip.chip_y - 1e6), pars_probe, pars_squid).place(self.cell, self.layer_ph, self.layer_el)

    def get_sps_params(self):
        pars = {'r_out'	:	175e3, # Radius of an outer ring including the empty region
                'dr'	:	25e3, # Gap in the outer ring
                'n_semiwaves'	:	2,
                's'	:	10e3, # Gap between two pads of a central capacitor
                'alpha'	:	pi/4, # period of a gap zigzag
                'r_curve'	:	30e3, # curvature of the roundings at the edges of a zigzag
                'n_pts_cwave'	:	200, # number of points for drawing a wave gap between to conductors
                'Z1'	:	self.Z_narrow, # Parameters of a top CPW
                'd_alpha1'	:	0, # width of a tip  of a central conductor of the top CPW
                'width1'	:	0, # width of a conductor in the top semiring
                'gap1'	:	25e3 - 1.33e3, # gap between the top semiring and the central capacitor
                'Z2'	:	self.Z, # Paramters of a bottom CPW
                'd_alpha2'	:	2 / 9 * pi, # length of a circumference covered by the bottom semiring
                'width2'	:	25e3/3, # width of a conductor in the bottom semiring
                'gap2'	:	25e3/3, # gap between the bottom semiring and the central capacitor
                'n_pts_arcs'	:	 50, # number of points for drawing a circle
                }
        return pars

    def get_mixing_qubit_params(self):
        pars = {'r_out'	:	175e3, # Radius of an outer ring including the empty region
                'dr'	:	25e3, # Gap in the outer ring
                'n_semiwaves'	:	2,
                's'	:	10e3, # Gap between two pads of a central capacitor
                'alpha'	:	pi/4, # period of a gap zigzag
                'r_curve'	:	30e3, # curvature of the roundings at the edges of a zigzag
                'n_pts_cwave'	:	200, # number of points for drawing a wave gap between to conductors
                'Z1'	:	self.Z_narrow, # Parameters of a top CPW
                'd_alpha1'	:	0, # width of a tip  of a central conductor of the top CPW
                'width1'	:	0, # width of a conductor in the top semiring
                'gap1'	:	25e3, # gap between the top semiring and the central capacitor
                'Z2'	:	self.Z, # Paramters of a bottom CPW
                'd_alpha2'	:	0, # length of a circumference covered by the bottom semiring
                'width2'	:	0, # width of a conductor in the bottom semiring
                'gap2'	:	25e3, # gap between the bottom semiring and the central capacitor
                'n_pts_arcs'	:	 50, # number of points for drawing a circle
                }
        return pars

    def get_dc_squid_params(self):
        pad_side = 5e3
        pad_r = 1e3
        pads_distance = 30e3
        p_ext_width = 3e3
        p_ext_r = 0.5e3
        sq_len = 7e3
        sq_area = 15e6
        j_width = 0.2e3
        low_lead_w = 0.5e3
        b_ext =   0.9e3
        j_length =  0.1e3
        n = 7
        bridge = 0.2e3
        return [pad_side, pad_r, pads_distance, p_ext_width,
                p_ext_r, sq_len, sq_area, j_width, low_lead_w,
                b_ext, j_length, n,bridge]

### MAIN FUNCTION ###
if __name__ == "__main__":
    my_design = My_Design('testScript')
    my_design.show()
    # my_design.save_as_gds2(r'C:\Users\andre\Documents\chip_designs\chip_design.gds2')