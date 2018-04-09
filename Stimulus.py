import numpy as np

class Stimulus(object):

    def __init__(self, t, dt, tstart, tstop, N, stimulus):
	
        self.t = t
        self.dt = dt
        self.tstart = tstart
        self.tstop = tstop
        self.N = N
        self.stimulus = stimulus
	
	
    def apply_stimulus(self):
	try:
		self.glut_input = self.stimulus[:,int(self.t/self.dt)]
	except IndexError:
		self.glut_input = self.stimulus[:,-1] #plus extra input for euler integration
