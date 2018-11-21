import numpy as np
import matplotlib.pyplot as plt

from Parameters import p

C=.000073
CER=.025
h=.7892
IP3=.15659

p['SVR'] = 1.
p['glut_input'] = 1.

r_er = .15 * np.exp(-.002*(p['SVR'])**2.32)
# membrane currents at the ER
ICERleak = ((p['F'])/(p['SVR']*np.sqrt(r_er))) * p['rl'] * (CER - C)
ISerca = ((p['F'])/(p['SVR']*np.sqrt(r_er))) * p['ver'] * C ** 2 / (C ** 2 + p['Ker'] ** 2)

m_infty = IP3 / (IP3 + p['d1'])
n_infty = C / (C + p['d5'])
IIP3R = ((p['F'])/(p['SVR']*np.sqrt(r_er))) * p['rc'] * m_infty ** 3 * n_infty ** 3 * h ** 3 \
		* (CER - C)


print "ICERleak="+str(ICERleak)
print "ISerca="+str(ISerca)
print "ICERleak="+str(ICERleak)

# IP3
K_gamma = p['Kr'] * (1 + p['Kp'] / p['Kr'] * C / (C + p['Kpi']))
v_glu = p['vb'] * p['glut_input'] ** 0.7 / (p['glut_input'] ** 0.7 + K_gamma ** 0.7)
v_3K = p['v3k'] * C ** 4 / (C ** 4 + p['KD'] ** 4) * IP3 / (IP3 + p['K3'])
v_delta = p['vd'] / (1 + IP3/p['kd']) \
		* C ** 2 / (C ** 2 + p['Kplcd'] ** 2)
prod_degr_IP3 = v_glu + v_delta - v_3K - p['r5p'] * IP3

# h
Q_2 = p['d2'] * (IP3 + p['d1']) / (IP3 + p['d3'])
h_infty = Q_2 / (Q_2 + C)
tau_h = 1 / (p['a2'] * (Q_2 + C))
dhdt = (h_infty - h) / tau_h

# transmembrane flux densities
J_CER_m = (- ISerca  + ICERleak + IIP3R)

dC = (p['SVR']*np.sqrt(r_er)/(p['F'])) * J_CER_m 


dCER = ((p['SVR']*np.sqrt(r_er))/(p['F']*r_er)) * (-J_CER_m) 





