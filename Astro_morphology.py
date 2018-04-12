from sys import exit
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

from Finite_Difference import *
from multi_comp_param import *
from multi_comp_init import *

class CalcMorpholy(object):
    def calc_volume(self):
        self.volume = (self.diameters/2.)**2 * self.length_comp * np.pi

    def calc_shell_surface(self):
        self.surface = self.diameters * self.length_comp * np.pi

    def calc_SVR(self):
        self.SVR = self.surface/self.volume

    def add_boundaries(self):
        # add outer boundaries to cross_section areas
        self.cross_plus_boundary = np.zeros(len(self.cross_area)+2)
        self.cross_plus_boundary[1:-1] = self.cross_area
        self.cross_plus_boundary[0] = self.cross_area[0]
        self.cross_plus_boundary[-1] = self.cross_area[-1]

class Process(CalcMorpholy):
    def __init__(self, diameter = None, length = None, length_comp = None, end_condition = None):
        '''
        Process of an astrocyte with constant diameter.
        Parameters
        ----------
        diameter :
            Diameter of the whole process. Default value is set to 1 * 1e-6 metre.
        length :
            Length of the process. Default value is set to 40 * 1e-6 metre.
        length_comp :
            Length of a single compartments within the process. Default value is set to
            0.5 * 1e-6.
        '''

        #default values
        self.diameter = 1. * 1e-6 #meter
        self.length = 40. * 1e-6 #meter
        self.length_comp = 0.5 * 1e-6 #meter, length of single compartment

        if diameter is not None:
            self.diameter = diameter

        if length is not None:
            self.length = length

        if length_comp is not None:
            self.length_comp = length_comp

        if end_condition is not None:
            self.end_condition = end_condition

        if end_condition == 'sealed_end':
            self.N = int(self.length/self.length_comp) # number of compartments
        elif end_condition == 'open_end':
            self.N = int(self.length/self.length_comp) #+ 2 # number of compartments

        self.cylinder_to_rectangles()
        self.calc_volume()
        self.calc_cross_area_cylinder()
        self.calc_shell_surface()
        self.calc_SVR()
        self.add_boundaries()

    def cylinder_to_rectangles(self):
        self.diameters = np.ones(self.N)*self.diameter

    def calc_cross_area_cylinder(self):
        radius = self.diameters[:-1]/2.
        self.cross_area = 2*(radius**2)*np.pi


class Soma(CalcMorpholy):
    def __init__(self, diameter = None, length = None, length_comp = None):
        '''
        Soma of an astrocyte consisting of cylinders with increasing/decreasing diameter.
        Parameters
        ----------
        diameter :
            Diameter of the whole process. Default value is set to 10 * 1e-6 metre.
        length :
            Length of the process. Default value is set to 10 * 1e-6 metre.
        length_comp :
            Length of a single compartments within the process. Default value is set to
            0.5 * 1e-6.
        '''

        #default values
        self.diameter = 10. * 1e-6 #meter
        self.length = 10. * 1e-6 #meter
        self.length_comp = 0.2 * 1e-6 #meter, length of singl compartment

        if diameter is not None:
            self.diameter = diameter

        if length is not None:
            self.length = length

        if length_comp is not None:
            self.length_comp = length_comp

        self.N = int(self.length/self.length_comp) # number of compartments

        self.x_y_coord_soma()
        self.soma_to_rectangles()
        self.calc_volume()
        self.calc_cross_area_soma()
        self.calc_shell_surface()
        self.calc_SVR()
        self.add_boundaries()

    def x_y_coord_soma(self):
        angle = np.arange(0,25)/12. * np.pi #angles described with radian
        self.x_soma = np.cos(angle)*self.diameter/2. # x coordinates of circle
        self.y_soma = np.sin(angle)*self.diameter/2. # y coordinates of circle

    def soma_to_rectangles(self):
        steps = self.N
        x_min = np.min(self.x_soma)
        x_max = np.max(self.x_soma)
        x_steps = np.linspace(x_min, x_max, steps+1)
        offset = (x_steps[1]-x_steps[0])/2.
        x_offset = np.linspace(x_min+offset,x_max-offset,steps)
        y_offset = np.sin(np.arccos(x_offset/(self.diameter/2.)))*self.diameter/2.

        self.x_coord = []
        self.y_coord = []
        self.diameters = []
        for i in range(len(x_offset)):
            x0 = x_steps[i]
            x1 = x_steps[i+1]
            y0 = y_offset[i]
            self.x_coord.append([x0, x1, x1, x0, x0])
            self.y_coord.append([-y0, -y0, y0, y0, -y0])
            self.diameters.append(2*y0)
        self.diameters = np.array(self.diameters)

    def calc_cross_area_soma(self):
        # check if it is symmetrical
        if len(self.diameters)%2 == 0:
            middle = len(self.diameters)/2
            start_from_middle = middle
        elif len(self.diameters)%2 == 1:
            middle = len(self.diameters)/2
            start_from_middle = middle+1
        if np.mean(np.round(np.sort(self.diameters[:middle]) - np.sort(self.diameters[start_from_middle:]),15)) == 0:
            radius = self.diameters[:middle]/2.
            cross_area_to_mid = 2*(radius**2)*np.pi
            if len(self.diameters)%2 == 0:
                self.cross_area = np.append(cross_area_to_mid, cross_area_to_mid[:middle-1][::-1])
            elif len(self.diameters_soma)%2 == 1:
                self.cross_area = np.append(cross_area_to_mid, cross_area_to_mid[::-1])
        else:
            print 'Not symmetrical morphology'


class Connect(object):
    def __init__(self, subcell_comp, connections = None):
        '''
        Connections between single subcellular compartments
        Parameters
        ----------
        subcell_comp :
            List of subcellular compartments (Endfeet, Process, Soma),
            which shall be connected.
        connections :
            List of lists of positions, which define the connections between the
            subcellular compartments.
        '''

        self.subcell_comp = subcell_comp

        if connections is not None:
            self.connections = connections

    def no_connection(self):
        self.N = self.subcell_comp[0].N
        self.length_comp = self.subcell_comp[0].length_comp
        self.cross_plus_boundary = self.subcell_comp[0].cross_plus_boundary
        self.volume = self.subcell_comp[0].volume
        self.SVR = self.subcell_comp[0].SVR
        self.length = self.subcell_comp[0].length
        self.diameters = self.subcell_comp[0].diameters
        self.end_condition = self.subcell_comp[0].end_condition

        conn = Finite_difference(N = self.N, length_comp = self.length_comp,
                       cross_plus_boundary = self.cross_plus_boundary,
                       volume = self.volume, end_condition = self.end_condition)
        conn.connectivity_matrix_simple()
        self.conn_matrix = conn.conn_matrix

    def add_branch(self):
        '''
        Function to add a branch consisting of two subcellular compartments
        (two processes, two endfeet) to one subcellular compartment. Checks for the
        cross section areas between the three subcellular compartments.
        '''
        if len(self.subcell_comp) != 3:
            print "Not enough subcellular compartments for a branching point."
            exit

        if len(self.subcell_comp) == 3:

            if self.subcell_comp[0].cross_area[0] > \
                self.subcell_comp[1].cross_area[0] + self.subcell_comp[2].cross_area[0]:
                self.N = 0
                self.length = 0
                self.cross_plus_boundary = []
                self.volume = []
                self.diameters = []
                self.SVR = []
                self.start_comp = []
                self.end_comp = []
                for ix, sc in enumerate(self.subcell_comp):
                    self.N += sc.N
                    self.length += sc.length
                    self.volume.append(sc.volume)
                    self.diameters.append(sc.diameters)
                    self.SVR.append(sc.SVR)
                    if ix == 0:
                        self.cross_plus_boundary.append(sc.cross_plus_boundary)
                        self.start_comp.append(0)
                        self.end_comp.append(sc.N-1)
                    elif ix > 0:
                        self.cross_plus_boundary.append(sc.cross_plus_boundary[1:])
                        self.start_comp.append(sc.N*ix)
                        self.end_comp.append(sc.N-1+sc.N*ix)
                self.cross_plus_boundary = np.concatenate(self.cross_plus_boundary)
                self.cross_plus_boundary[self.subcell_comp[0].N] = self.subcell_comp[1].cross_area[0] + \
                                                self.subcell_comp[2].cross_area[0]
                self.volume = np.concatenate(self.volume)
                self.diameters = np.concatenate(self.diameters)
                self.SVR = np.concatenate(self.SVR)
                self.length_comp = sc.length_comp
                self.end_condition = sc.end_condition

                conn = Finite_difference(N = self.N, length_comp = self.length_comp,
                               cross_plus_boundary = self.cross_plus_boundary, volume = self.volume,
                               end_condition = end_condition,
                               start_comp = start_comp, end_comp = self.end_comp)
                conn.connectivity_matrix_branch()
                self.conn_matrix = conn.conn_matrix


    def add_subcell_comp(self):
        '''
        Function to add a another subcellular compartment to an existing
         subcellular comprtment.
        '''

        self.N = 0
        self.length = 0
        self.cross_plus_boundary = []
        self.volume = []
        self.SVR = []
        self.start_comp = []
        self.end_comp = []
        for ix, sc in enumerate(self.subcell_comp):
            self.N += sc.N
            self.length += sc.length
            self.volume.append(sc.volume)
            self.SVR.append(sc.SVR)
            if ix == 0:
                self.cross_plus_boundary.append(sc.cross_plus_boundary)
                self.start_comp.append(0)
                self.end_comp.append(sc.N-1)
            elif ix > 0 and len(self.subcell_comp) == 2:
                self.cross_plus_boundary.append(sc.cross_plus_boundary[1:])
                self.start_comp.append(sc.N*ix)
                self.end_comp.append(sc.N-1+sc.N*ix)
            elif ix > 0 and ix < len(self.subcell_comp)-1:
                self.cross_plus_boundary.append(sc.cross_plus_boundary[1:-1])
                self.start_comp.append(sc.N*ix)
                self.end_comp.append(sc.N-1+sc.N*ix)
            elif ix == len(self.subcell_comp)-1:
                self.cross_plus_boundary.append(sc.cross_plus_boundary)
                self.start_comp.append(sc.N*ix)
                self.end_comp.append(sc.N-1+sc.N*ix)
        self.cross_plus_boundary = np.concatenate(self.cross_plus_boundary)
        self.volume = np.concatenate(self.volume)
        self.SVR = np.concatenate(self.SVR)
        self.length_comp = sc.length_comp

        conn = Finite_difference(N = self.N, length_comp = self.length_comp,
                       cross_plus_boundary = self.cross_plus_boundary,
                       volume = self.volume, end_condition = self.end_condition)
        self.conn_matrix = conn.conn_matrix








