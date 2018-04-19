__author__ = 'franziska'

import matplotlib.pyplot as plt
from matplotlib import rcParams
import numpy as np
import os

#rcParams.update(params)

from Analysis import Astro_analysis
from Astro import Astro_multi_compartment
from Astro_morphology import Process, Soma, Connect
from helper_plots import *
from Parameters import p
import Stimulus
from Stimulus_gen import Stimulus_gen

# ratio of internal calcium stores
p['ratio'] = 0.15 * np.exp(-0.002*(1**2.32))

# define morphology
p1 = Process(end_condition = 'sealed_end')
p2 = Process(diameter = 0.5 * 1e-6, end_condition = 'sealed_end')
p3 = Process(diameter = 0.5 * 1e-6, end_condition = 'sealed_end')

morpho = Connect([p1, p2, p3])
morpho.add_branch()

p['time'] = 300
p['tstart'] = 100
p['tstop'] = 200

# generate stimulus
glut_conc = 0.1
stim = Stimulus_gen(morpho.N, p['time'], p['dt'],
                          p['tstart'], p['tstop'],
                          p['input_length'], morpho.length, glut_conc, comp_start = 70)
stim.generate_stimulus()

# simulate system
p['I_GluT_max'] = 4
p['P_max'] = 0.4 * 1.52/p['F']
p['D_C'] = 0.223 * 1e-11 #meter**2 * second **-1

astro = Astro_multi_compartment(params = p, model_type='NKV_Calcium',
                                stimulus = stim.glut_stimulus, morpho = morpho)

plt.figure(figsize = (3.087,3.087))
ax = plt.subplot(111)
plt.plot(np.linspace(-p1.length, p2.length, p1.N + p2.N)*1e6, astro.Na[2000,:160])
plt.xlabel(r'$\mathsf{x\/\/[\mu m]}$')
plt.ylabel(r'$\mathsf{Na^{+}\/\/[mM]}$')
plt.ylim(14,33)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.savefig('figures/branching_stimulation_bigger_branch.svg')