# generate dictionary for parameters
# parameter values from Oschmann et al., 2017
p = {}

# simulation duration and timestep
p['time'] = 300. #seconds
p['dt'] = .1 #seconds
p['tstart'] = 100.
p['tstop'] = 110.
p['input_length'] = 1. * 1e-6 #length of input zone

# initial concentrations
p['Na_0']    =  15.# + 0.189 #mole/meter**3
p['Na_o_0']          = 150.# - 0.378  #mole/meter**3
p['K_0']          = 100.#99.79088145 #mole/meter**3
p['K_o_0']          = 3.#41823711 #mole/meter**3
#p['V_0'] = (-85. - 1.4) * 1e-3 #V
p['C_0']          = 7.3 * 1e-5 #mole/meter**3
p['C_o_0']          = 1.8 #mole/meter**3


# physical constants
p['F'] = 96500.0 # coulomb / mole
p['R'] =  8.315 # joule / (mole*kelvin)
p['T'] = 298. # kelvin
p['psifac'] = ((p['R']*p['T'])/p['F'])

# geometrical parameters
#p['N'] = 101
#p['length'] = 2.5 * 1e-4 #m
#p['h'] = p['length']/p['N'] #length of each compartment
#p['SVR'] = 1. * 1e6 #meter**-1, surface volume ratio of the astrocyte
p['a_i'] = 0.4 # Tissue volume fraction being astrocytes
#p['O_m'] = p['SVR'] * p['a_i']#a_i/(5.00 * 1e-8) # meter, Astrocytic membrane area per tissue volume
p['a_o'] = 0.2 #(Tissue volume fraction being ECS)
#p['SVR_i'] = p['O_m']/p['a_i']
#p['SVR_o'] = p['O_m']/p['a_o']

# membrane parameters
#p['gK'] = 16.96 # S/m**2
#p['gN'] = 1. # S/m**2
p['C_m'] = 1.0 *1e-2 #farad/meter**2
p['rL'] = 5 * 1e15

# diffusion coefficients
p['D_Na'] = 0.6 * 1e-9 #meter**2 * second **-1 1.33*1e-9
p['D_C'] = 0.223 * 1e-9 #meter**2 * second **-1
p['D_K'] = 1. * 1e-9 #meter**2 * second **-1 1.96*1e-9
p['D_IP3'] = 0.28 * 1e-9 #meter**2 * second **-1

# tortuosity
p['lamb_intra'] = 3.2
p['lamb_extra'] = 1.6

# valence
p['z_C'] = 2.
p['z_Na'] = 1.
p['z_K'] = 1.

#iglu
p['K_mN_glu'] = 15. #mole/meter**3
p['K_mK_glu'] = 5. #mole/meter**3
p['K_mg'] = 34. * 1e-3 #mole/meter**3
p['I_GluT_max'] = .68 # amp/meter**2

#INCX
p['I_NCX_max']    =  0.0001 # amp/meter**2
p['K_mN']     = 87.5 #mole/meter**3
p['K_mC']     =  1.380 #mole/meter**3
p['eta']      =     0.35
p['k_sat']    =     0.1

# P (NKA)
p['K_mN_NKA'] = 10. #mole/meter**3
p['K_mK'] = 1.5 #mole/meter**3
p['P_max'] = 1.52/p['F'] #mol/(meter**2 * second) #1.115*1e-6

# i_serca
p['Ker']      = 0.1 * 1e-3 #mole/meter**3
p['ver']      =  5*0.9*1e-3#0.9 * 1e-3 #0.9 * 1e-3 #mole/(meter**3 * second) #44

# i_ip3 - CICR through I_IP3
p['d1']       = 0.13 * 1e-3 #mole/meter**3, IP3 dissociation constant
p['d5']       = 0.08234 * 1e-3 #mole/meter**3, Ca2+ activation dissociation constant
p['rc']       = 0.5*6. # second **-1

# i_ca_er-leak
p['rl']    = 0.5*0.11 # second **-1 0.55

# h variable
p['d2']    = 1.049  * 1e-3 #mole/meter**3, Ca2+ inactivation dissociation constant
p['d3']    = 0.9434 * 1e-3 #mole/meter**3, IP3 dissociation constant
p['a2']    = 0.2 * 1e3 #1/(msecond), IP3R binding rate for Ca2+ inhibition

# IP3 variable
# IP3 production - PLC_beat
p['vb']    =  0.1 * 1e-3 #mole/(meter**3 * second)
p['Kr']    =  1.3 * 1e-3 #mole/meter**3
p['Kp']    = 10.0 * 1e-3 #mole/meter**3
p['Kpi']   =  0.6 * 1e-3 #mole/meter**3
# IP3 production - PLC_delta
p['vd']    = 0.02 * 1e-3 #mole/(meter**3 * second)
p['kd']    = 1.5  * 1e-3 #mole/meter**3, Inhibition constant of PLCdelta activity
p['Kplcd'] = 0.1  * 1e-3 #mole/meter**3, Ca2+ affinity of PLCdetla
# IP3 degradation - IP3-3K
p['v3k']   = 2.0 * 1e-3 #mole/(meter**3 * second), Maximal rate of degradation by IP3-3K
p['K3']    = 1.0 * 1e-3 #mole/meter**3, IP3 affinity of IP3-3K
p['KD']    = 0.7 * 1e-3 #mole/meter**3, Ca2+ affinity of IP3-3K
# IP3 degradation- IP-5P
p['r5p']   = 0.04 # second **-1


# external bath
p['k_dec_N'] = 0*1e-7 #Output factor (m/s)
p['k_dec_K'] = 0*1e-7 #Output factor (m/s)










