import numpy as np
from scipy import signal

class Stimulus_gen(object):

    def __init__(self, N, time, dt, tstart, tstop, input_length, length, glut_conc,
                 stimulus_type = 'constant', num_pulses = None, comp_start = None):
        self.N = N
        self.time = time
        self.dt = dt
        self.tstart = tstart
        self.tstop = tstop
        self.input_length = input_length
        self.length = length
        self.glut_conc = glut_conc
        self.stimulus_type = stimulus_type
        if num_pulses is not None:
            self.num_pulses = num_pulses
        self.comp_start = comp_start


    def generate_stimulus(self):
        self.input_zone = int(np.ceil(self.input_length/(self.length/self.N))) #calculates number of compartments in input zone
        self.stimulus = np.zeros((self.N, int(self.time/self.dt)))
        if self.comp_start is not None:
            self.comp_start = self.comp_start
        elif self.comp_start is None:
            self.comp_start = int((self.N - self.input_zone)/2)
        self.comp_stop = self.comp_start + self.input_zone
        self.stimulus[self.comp_start:self.comp_stop,int(self.tstart/self.dt):int(self.tstop/self.dt)] = self.glut_conc

        if self.stimulus_type == 'constant':
            self.glut_stimulus = self.stimulus
        elif self.stimulus_type == 'pulse':
            #pulses = (np.arange(int((self.tstop - self.tstart)/self.dt)) % self.pulse_period < self.pulse_dur)*self.glut_conc
            t = np.linspace(0, 1, int((self.tstop - self.tstart)/self.dt), endpoint=False)
            pulses = (signal.square(2*np.pi*self.num_pulses*t)+1.)*((1./2)*self.glut_conc)
            pulses_ar = np.ones((self.input_zone, len(pulses)))*pulses
            self.stimulus[self.comp_start:self.comp_stop,int(self.tstart/self.dt):int(self.tstop/self.dt)] = pulses_ar
            self.glut_stimulus = self.stimulus
        elif self.stimulus_type == 'constant_plus_noise':
            noise = np.round(np.random.normal(0,0.01,int((self.tstop - self.tstart)/self.dt)),4)
            noise_ar = np.ones((self.input_zone, len(noise)))*noise
            self.stimulus[self.comp_start:self.comp_stop,int(self.tstart/self.dt):int(self.tstop/self.dt)] += noise_ar
            self.glut_stimulus = self.stimulus
