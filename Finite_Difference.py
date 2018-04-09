#from Astro_morphology import Astrocyte_morphology
from Parameters import p

import numpy as np
import matplotlib.pyplot as plt

class Finite_difference(object):

    def __init__(self, N = None, length_comp = None,
                 cross_plus_boundary = None, volume = None, end_condition = None,
                 start_comp = None, end_comp = None):

        if N is not None:
            self.N = N
        if length_comp is not None:
            self.length_comp = length_comp
        if cross_plus_boundary is not None:
            self.cross_plus_boundary = cross_plus_boundary
        if volume is not None:
            self.volume = volume
        if start_comp is not None:
            self.start_comp = start_comp
        if end_comp is not None:
            self.end_comp = end_comp

        self.end_condition = end_condition
        #if self.end_condition == 'open_end':
             #self.N_max = self.N
             #self.N = self.N - 2

            # #volume for end compartments
            # vol = np.zeros(len(self.volume)+2)
            # vol[1:-1] = self.volume
            # vol[0] = self.volume[0]
            # vol[-1] = self.volume[-1]
            # self.volume = vol
            #
            # cross = np.zeros(len(self.cross_plus_boundary)+2)
            # cross[1:-1] = self.cross_plus_boundary
            # cross[0] = self.cross_plus_boundary[0]
            # cross[-1] = self.cross_plus_boundary[-1]
            # self.cross_plus_boundary = cross

        self.diffusion_coeff_scaled_right = self.volume/self.cross_plus_boundary[1:]
        self.diffusion_coeff_scaled_left = self.volume/self.cross_plus_boundary[:-1]


    def forward_difference(self):
        self.for_diff = np.zeros((self.N,self.N)) # matrix for forward difference

        #input from left
        for i in range(self.N):
            self.for_diff[i,i] = -1
            self.for_diff[i,i] *= self.diffusion_coeff_scaled_left[i] * \
                                  self.cross_plus_boundary[i]/self.volume[i]

        #input from right
        for i in range(self.N-1):
            self.for_diff[i,i+1] = 1
            self.for_diff[i,i+1] *= self.diffusion_coeff_scaled_right[i] * \
                                    self.cross_plus_boundary[i+1]/self.volume[i]

        #sealed end cond for last comp
        if self.end_condition == 'sealed_end':
            self.for_diff[-1,-1] = 0

        #divide by length of single compartment
        self.for_diff /= self.length_comp


    def backward_difference(self):
        self.back_diff = np.zeros((self.N,self.N)) # matrix for backward difference

        for i in range(self.N):
            self.back_diff[i,i] = 1

        for i in range(1,self.N):
            self.back_diff[i,i-1] = -1

        #divide by length of single compartment
        self.back_diff /= self.length_comp


    def connectivity_matrix_simple(self):

        # set up connectivity matrix
        self.forward_difference()
        self.backward_difference()
        self.conn_matrix = self.back_diff.dot(self.for_diff)

        if self.end_condition == 'open_end':
            conn = np.zeros((self.N, self.N+2))
            conn[:,1:-1] = self.conn_matrix
            conn[0,0] = conn[0,2]
            conn[-1,-1] = conn[-1,-3]
            conn[0,1] *= 2
            self.conn_matrix = conn


    def connectivity_matrix_branch(self):

        # set up connectivity matrix
        self.forward_difference()
        self.backward_difference()
        self.conn_matrix = self.back_diff.dot(self.for_diff)

        # set connection between first and third compartment
        # (based on connection between first and second cmpartment)
        self.conn_matrix[self.end_comp[0],self.start_comp[2]] = \
            self.conn_matrix[self.end_comp[0],self.start_comp[1]]

        self.conn_matrix[self.start_comp[2],self.end_comp[0]] = \
            self.conn_matrix[self.end_comp[0],self.start_comp[1]]

        # condition to meet the Kirchhoffs law for the branching compartment
        self.conn_matrix[self.end_comp[0],self.end_comp[0]] *= 1.5

        # block conn between second and third compartment
        self.conn_matrix[self.end_comp[1],self.end_comp[1]] *= 0.5
        self.conn_matrix[self.end_comp[1],self.start_comp[2]] = 0
        self.conn_matrix[self.start_comp[2],self.end_comp[1]] = 0

