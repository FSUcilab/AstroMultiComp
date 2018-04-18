import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

from Astro_morphology import Process
from Astro_morphology import Connect
from Init import Init
from Parameters import p
from Finite_Difference import Finite_difference

class Astro_multi_compartment(object):

    def __init__(self, params, model_type, stimulus, morpho, diff_range = None):

        '''
        Simulation of an astrocytic subcellular compartment consisting of a multi-compartment model
        ----------
        params :
            Dictionary with parameters.
        model_type:
            Specifies the considered currents of the model.
            NKV: Membrane currents affecting the concentrations of Na+ and K+ and
                the membrane voltage.
            Calcium: Ca2+ currents at the internal Ca2+ store.
            NKV_Calcium: All currents affecting the concentration of Na+, Ca2+, K+ and the membrane voltage.
        stimulus:
            Matrix containing the concentration on the stimulus in space and time.
        morpho:
            Object produced by the Astro_morphology class.
        '''

        # converts dict with parameters to variables
        self.__dict__.update(params)

        self.model_type = model_type

        self.stimulus = stimulus

        self.conn_matrix = morpho.conn_matrix

        # morphology
        self.SVR = morpho.SVR
        self.N = morpho.N
        self.diameters = morpho.diameters
        self.length_comp = morpho.length_comp
        self.end_condition = morpho.end_condition

        # duration and time step of simulation
        self.tspan = np.arange(0,self.time, self.dt)

        # calculate initial values
        self.init = Init(p)
        self.init.initialize()

        # duration and time step of simulation
        self.tspan = np.arange(0,self.time, self.dt)

        # set stimulus duration and length of input zone
        self.comp_start = int(self.N * 0.49)
        self.comp_end = int(self.N * 0.51)

        # parameters for diff_range
        if diff_range is not None:
            self.d_cond, self.d_start, self.d_end, self.d_coeff = ['step', 70, 79, self.D_C*0.1]

        if self.model_type == 'NKV':

            #calculate initial values
            self.V_0 = self.init.V_0
            self.gN = self.init.gN
            self.gK = self.init.gK

            # initialize the system
            init_Na = self.Na_0 * np.ones(self.N) # mM
            init_K = self.K_0 * np.ones(self.N) # mM
            init_C = self.C_0 * np.ones(self.N) # Mol
            init_Na_o = self.Na_o_0 * np.ones(self.N) # mM
            init_K_o = self.K_o_0 * np.ones(self.N) # mM
            init_C_o = self.C_o_0 * np.ones(self.N) # Mol

            init = np.hstack((init_Na, init_K, init_C, init_Na_o, init_K_o, init_C_o))

            # simulate spatial astro
            self.sol = odeint(self.spatial_astro_NKV, init, self.tspan, tcrit=[self.tstart, self.tstop])

            # transfer solution of ode system into single variables
            self.Na = self.sol[:,:self.N]
            self.K = self.sol[:,self.N:2*self.N]
            self.C = self.sol[:,2*self.N:3*self.N]
            self.Na_o = self.sol[:,3*self.N:4*self.N]
            self.K_o = self.sol[:,4*self.N:5*self.N]
            self.C_o = self.sol[:,5*self.N:6*self.N]
            #self.V = self.sol[:,6*self.N:7*self.N]

        elif self.model_type == 'Calcium':

            #calculate initial values
            self.IP3_0 = self.init.IP3_0
            self.h_0 = self.init.h_0
            self.CER_0 = self.init.CER_0

            # initialize the system
            init_C = self.C_0 * np.ones(self.N) # Mol
            init_IP3 = self.IP3_0 * np.ones(self.N) # mM
            init_h = self.h_0 * np.ones(self.N) # mM
            init_CER = self.CER_0 * np.ones(self.N) # mM
            init = np.hstack((init_C, init_IP3, init_h, init_CER))

            # simulate spatial astro
            self.sol = odeint(self.spatial_astro_Ca, init, self.tspan, tcrit=[self.tstart, self.tstop])

            # transfer solution od ode system into single variables
            self.C = self.sol[:,:self.N]
            self.IP3 = self.sol[:,self.N:2*self.N]
            self.h = self.sol[:,2*self.N:3*self.N]
            self.CER = self.sol[:,3*self.N:4*self.N]

        elif self.model_type == 'NKV_Calcium':

            # calculate initial values
            self.V_0 = self.init.V_0
            self.gN = self.init.gN
            self.gK = self.init.gK
            self.IP3_0 = self.init.IP3_0
            self.h_0 = self.init.h_0
            self.CER_0 = self.init.CER_0

            # initialize the system
            init_Na = self.Na_0 * np.ones(self.N) # mM
            init_K = self.K_0 * np.ones(self.N) # mM
            init_C = self.C_0 * np.ones(self.N) # Mol

            init_IP3 = self.IP3_0 * np.ones(self.N) # mM
            init_h = self.h_0 * np.ones(self.N) # mM
            init_CER = self.CER_0 * np.ones(self.N) # mM

            init_Na_o = self.Na_o_0 * np.ones(self.N) # mM
            init_K_o = self.K_o_0 * np.ones(self.N) # mM
            init_C_o = self.C_o_0 * np.ones(self.N) # Mol

            init = np.hstack((init_Na, init_K, init_C, init_IP3, init_h, init_CER,
                              init_Na_o, init_K_o, init_C_o))

            self.diffusion_range(condition=self.d_cond, start = self.d_start, end = self.d_end, d_coeff=self.d_coeff)

            # simulate spatial astro
            self.sol = odeint(self.spatial_astro_NKV_Ca, init, self.tspan, tcrit=[self.tstart, self.tstop])

            # transfer solution od ode system into single variables
            self.Na = self.sol[:,:self.N]
            self.K = self.sol[:,self.N:2*self.N]
            self.C = self.sol[:,2*self.N:3*self.N]

            self.IP3 = self.sol[:,3*self.N:4*self.N]
            self.h = self.sol[:,4*self.N:5*self.N]
            self.CER = self.sol[:,5*self.N:6*self.N]

            self.Na_o = self.sol[:,6*self.N:7*self.N]
            self.K_o = self.sol[:,7*self.N:8*self.N]
            self.C_o = self.sol[:,8*self.N:9*self.N]


    def open_end(self, conc, init_conc):
        tmp = np.zeros(len(conc)+2)
        tmp[1:-1] = conc
        tmp[0] = tmp[-1] = init_conc
        return tmp

    def diffusion_range(self, condition = 'equal', start = None, end = None, d_coeff = None):
        if condition == 'equal':
            self.D_Ci = self.D_C
        elif condition == 'step':
            self.D_Ci = self.D_C * np.ones(self.N)
            self.D_Ci[start:end] = d_coeff

    def spatial_astro_NKV(self, state, tspan):

        # initial values
        Na = state[:self.N]
        K = state[self.N:2*self.N]
        C = state[2*self.N:3*self.N]
        Nao = state[3*self.N:4*self.N]
        Ko = state[4*self.N:5*self.N]
        Co = state[5*self.N:6*self.N]
        #V = state[6*self.N:7*self.N]

        # input
        if tspan < self.time:
            self.glut_input = self.stimulus[:,int(int(tspan)/self.dt)]#self.input(tspan)
        elif tspan >= self.time:
            self.glut_input = self.stimulus[:,-1]

        # Astrocytic membrane area per tissue volume
        self.O_m = self.SVR * self.a_i

        # membrane potential
        X_oZ_o = -((self.O_m*self.C_m*self.V_0)/(self.a_o)) - self.F*(self.z_K*self.K_o_0 + self.z_Na*self.Na_o_0 + self.z_C*self.C_o_0)
        V = -((self.a_o)/(self.C_m * self.O_m))*((self.F*(self.z_K * Ko + self.z_Na * Nao + self.z_C * Co)) + X_oZ_o)

        # transmembrane currents
        IGluT = self.I_GluT_max * ((K)/(K + self.K_mK_glu)) * ((Nao)**3/(((Nao)**3) + self.K_mN_glu**3)) * (self.glut_input/(self.glut_input + self.K_mg))
        INCX = self.I_NCX_max * (Nao)**3/(self.K_mN**3+(Nao)**3) * (Co/(self.K_mC+Co)) * \
               (np.exp((self.eta) * V*self.F/(self.R*self.T))*((Na)**3/(Nao)**3) - np.exp((self.eta - 1) * V*self.F/(self.R*self.T))*(C/(Co)))/(1. + self.k_sat*np.exp((self.eta - 1)*V*self.F/(self.R*self.T)))
        INKA = self.P_max * ((Na**1.5)/(Na**1.5 + self.K_mN_NKA**1.5)) * ((Ko)/(Ko + self.K_mK))
        IKleak = ((self.gK) * (V - self.psifac*np.log(Ko/K)))
        INleak = (((self.gN) * (V - self.psifac*np.log(Nao/Na))))

        # transmembrane flux densities
        J_Ca_m = (-1.*INCX)/(self.F)
        J_K_m = (IKleak - (2.*INKA) + (1.*IGluT/self.F))
        J_Na_m = (INleak + (3.*INKA) - (3.*IGluT/self.F) + (3.*INCX/self.F))

        if self.end_condition == 'open_end':
            C = self.open_end(C, self.C_0)
            Co = self.open_end(Co, self.C_o_0)
            K = self.open_end(K, self.K_0)
            Ko = self.open_end(Ko, self.K_o_0)
            Na = self.open_end(Na, self.Na_0)
            Nao = self.open_end(Nao, self.Na_o_0)

        # diffusive flux
        J_CiD = self.conn_matrix.dot(-(self.D_C/(self.lamb_intra**2)) * C)
        J_CoD = self.conn_matrix.dot(-(self.D_C/(self.lamb_extra**2)) * Co)
        J_KiD = self.conn_matrix.dot(-(self.D_K/(self.lamb_intra**2)) * K)
        J_NaiD = self.conn_matrix.dot(-(self.D_Na/(self.lamb_intra**2)) * Na)
        J_KoD = self.conn_matrix.dot(-(self.D_K/(self.lamb_extra**2)) * Ko)
        J_NaoD = self.conn_matrix.dot(-(self.D_Na/(self.lamb_extra**2)) * Nao)

        #define differential equations
        dKdt = -(self.O_m/self.a_i)*(J_K_m) - J_KiD
        dNadt = -(self.O_m/self.a_i)*(J_Na_m) - J_NaiD
        dKodt = (self.O_m/self.a_o)*(J_K_m) - J_KoD
        dNaodt = (self.O_m/self.a_o)*(J_Na_m) - J_NaoD
        dCdt = -(self.O_m/self.a_i) * J_Ca_m - J_CiD
        dCodt = (self.O_m/self.a_o) * J_Ca_m - J_CoD


        return np.hstack((dNadt, dKdt, dCdt, dNaodt, dKodt, dCodt))


    def spatial_astro_Ca(self, state, tspan):

        # initial values
        C = state[:self.N]
        IP3 = state[self.N:2*self.N]
        h = state[2*self.N:3*self.N]
        CER = state[3*self.N:4*self.N]

        # Astrocytic membrane area per tissue volume
        self.O_m = self.SVR * self.a_i

        # input
        if tspan < self.time:
            self.glut_input = self.stimulus[:,int(int(tspan)/self.dt)]#self.input(tspan)
        elif tspan >= self.time:
            self.glut_input = self.stimulus[:,-1]

        # membrane currents at the ER
        ICERleak = ((self.F)/(self.SVR*np.sqrt(self.ratio))) * self.rl * (CER - C)
        ISerca = ((self.F)/(self.SVR*np.sqrt(self.ratio))) * self.ver * C ** 2 / (C ** 2 + self.Ker ** 2)

        m_infty = IP3 / (IP3 + self.d1)
        n_infty = C / (C + self.d5)
        IIP3R = ((self.F)/(self.SVR*np.sqrt(self.ratio))) * self.rc * m_infty ** 3 * n_infty ** 3 * h ** 3 \
                * (CER - C)

        # IP3
        K_gamma = self.Kr * (1 + self.Kp / self.Kr * C / (C + self.Kpi))
        v_glu = self.vb * self.glut_input ** 0.7 / (self.glut_input ** 0.7 + K_gamma ** 0.7)
        v_3K = self.v3k * C ** 4 / (C ** 4 + self.KD ** 4) * IP3 / (IP3 + self.K3)
        v_delta = self.vd / (1 + IP3/self.kd) \
                * C ** 2 / (C ** 2 + self.Kplcd ** 2)
        prod_degr_IP3 = v_glu + v_delta - v_3K - self.r5p * IP3

        # h
        Q_2 = self.d2 * (IP3 + self.d1) / (IP3 + self.d3)
        h_infty = Q_2 / (Q_2 + C)
        tau_h = 1 / (self.a2 * (Q_2 + C))
        dhdt = (h_infty - h) / tau_h

        # transmembrane flux densities
        J_CER_m = (- ISerca  + ICERleak + IIP3R)

        if self.end_condition == 'open_end':
            C = self.open_end(C, self.C_0)
            CER = self.open_end(CER, self.C_ER_0)
            IP3 = self.open_end(IP3, self.IP3_0)

        # diffusive flux
        J_CiD = self.conn_matrix.dot(-(self.D_C/(self.lamb_intra**2)) * C)
        J_CERiD = self.conn_matrix.dot(-(self.D_C/(self.lamb_intra**2)) * CER)
        J_IP3iD = self.conn_matrix.dot(-(self.D_IP3/(self.lamb_intra**2)) * IP3)


        dCdt = ((self.SVR*np.sqrt(self.ratio))/(self.F)) * J_CER_m - J_CiD
        dIP3dt = prod_degr_IP3 - J_IP3iD
        dCERdt = ((self.SVR*np.sqrt(self.ratio))/(self.F*self.ratio)) * (-J_CER_m) - J_CERiD

        return np.hstack((dCdt, dIP3dt, dhdt, dCERdt))

    def spatial_astro_NKV_Ca(self, state, tspan):

        # initial values
        Na = state[:self.N]
        K = state[self.N:2*self.N]
        C = state[2*self.N:3*self.N]

        IP3 = state[3*self.N:4*self.N]
        h = state[4*self.N:5*self.N]
        CER = state[5*self.N:6*self.N]

        Nao = state[6*self.N:7*self.N]
        Ko = state[7*self.N:8*self.N]
        Co = state[8*self.N:9*self.N]

        # Astrocytic membrane area per tissue volume
        self.O_m = self.SVR * self.a_i

        # membrane potential
        X_oZ_o = -((self.O_m*self.C_m*self.V_0)/(self.a_o)) - self.F*(self.z_K*self.K_o_0 + self.z_Na*self.Na_o_0 + self.z_C*self.C_o_0)
        V = -((self.a_o)/(self.C_m * self.O_m))*((self.F*(self.z_K * Ko + self.z_Na * Nao + self.z_C * Co)) + X_oZ_o)

        # input
        if tspan < self.time:
            self.glut_input = self.stimulus[:,int(int(tspan)/self.dt)]#self.input(tspan)
        elif tspan >= self.time:
            self.glut_input = self.stimulus[:,-1]

        # transmembrane currents
        IGluT = self.I_GluT_max * ((K)/(K + self.K_mK_glu)) * ((Nao)**3/(((Nao)**3) + self.K_mN_glu**3)) * (self.glut_input/(self.glut_input + self.K_mg))
        INCX = self.I_NCX_max * (Nao)**3/(self.K_mN**3+(Nao)**3) * (Co/(self.K_mC+Co)) * \
               (np.exp((self.eta) * V*self.F/(self.R*self.T))*((Na)**3/(Nao)**3) - np.exp((self.eta - 1) * V*self.F/(self.R*self.T))*(C/(Co)))/(1. + self.k_sat*np.exp((self.eta - 1)*V*self.F/(self.R*self.T)))
        INKA = self.P_max * ((Na**1.5)/(Na**1.5 + self.K_mN_NKA**1.5)) * ((Ko)/(Ko + self.K_mK))
        E_K = self.psifac*np.log(Ko/K)
        E_Na = self.psifac*np.log(Nao/Na)

        # membrane currents at the ER
        ICERleak = ((self.F)/(self.SVR*np.sqrt(self.ratio))) * self.rl * (CER - C)
        ISerca = ((self.F)/(self.SVR*np.sqrt(self.ratio))) * self.ver * C ** 2 / (C ** 2 + self.Ker ** 2)

        m_infty = IP3 / (IP3 + self.d1)
        n_infty = C / (C + self.d5)
        IIP3R = ((self.F)/(self.SVR*np.sqrt(self.ratio))) * self.rc * m_infty ** 3 * n_infty ** 3 * h ** 3 \
                * (CER - C)

        # IP3
        K_gamma = self.Kr * (1 + self.Kp / self.Kr * C / (C + self.Kpi))
        v_glu = self.vb * self.glut_input ** 0.7 / (self.glut_input ** 0.7 + K_gamma ** 0.7)
        v_3K = self.v3k * C ** 4 / (C ** 4 + self.KD ** 4) * IP3 / (IP3 + self.K3)
        v_delta = self.vd / (1 + IP3/self.kd) \
                * C ** 2 / (C ** 2 + self.Kplcd ** 2)
        prod_degr_IP3 = v_glu + v_delta - v_3K - self.r5p * IP3

        # h
        Q_2 = self.d2 * (IP3 + self.d1) / (IP3 + self.d3)
        h_infty = Q_2 / (Q_2 + C)
        tau_h = 1 / (self.a2 * (Q_2 + C))
        dhdt = (h_infty - h) / tau_h

        # transmembrane flux densities
        J_Na_m = (((self.gN) * (V - E_Na)) + (3.*INKA) - (3.*IGluT/self.F) + (3.*INCX/self.F))
        J_K_m = (((self.gK) * (V - E_K)) - (2.*INKA) + (1.*IGluT/self.F))
        J_Ca_m = (-1.*INCX)/(self.F)
        J_CER_m = (- ISerca  + ICERleak + IIP3R)

        if self.end_condition == 'open_end':
            C = self.open_end(C, self.C_0)
            Co = self.open_end(Co, self.C_o_0)
            K = self.open_end(K, self.K_0)
            Ko = self.open_end(Ko, self.K_o_0)
            Na = self.open_end(Na, self.Na_0)
            Nao = self.open_end(Nao, self.Na_o_0)
            CER = self.open_end(CER, self.C_ER_0)
            IP3 = self.open_end(IP3, self.IP3_0)

        # diffusive flux
        J_CiD = self.conn_matrix.dot(-(self.D_Ci/(self.lamb_intra**2)) * C)
        J_CoD = self.conn_matrix.dot(-(self.D_C/(self.lamb_extra**2)) * Co)
        J_KiD = self.conn_matrix.dot(-(self.D_K/(self.lamb_intra**2)) * K)
        J_NaiD = self.conn_matrix.dot(-(self.D_Na/(self.lamb_intra**2)) * Na)
        J_KoD = self.conn_matrix.dot(-(self.D_K/(self.lamb_extra**2)) * Ko)
        J_NaoD = self.conn_matrix.dot(-(self.D_Na/(self.lamb_extra**2)) * Nao)
        J_CERiD = self.conn_matrix.dot(-(self.D_C/(self.lamb_intra**2)) * CER)
        J_IP3iD = self.conn_matrix.dot(-(self.D_IP3/(self.lamb_intra**2)) * IP3)

        #define differential equations
        dKdt = -(self.O_m/self.a_i)*(J_K_m) - J_KiD
        dNadt = -(self.O_m/self.a_i)*(J_Na_m) - J_NaiD
        dKodt = (self.O_m/self.a_o)*(J_K_m) - J_KoD
        dNaodt = (self.O_m/self.a_o)*(J_Na_m) - J_NaoD
        dCdt = -(self.O_m/self.a_i) * J_Ca_m + ((self.SVR*np.sqrt(self.ratio))/(self.F)) * J_CER_m - J_CiD
        dCodt = (self.O_m/self.a_o) * J_Ca_m - J_CoD
        dIP3dt = prod_degr_IP3 - J_IP3iD
        dCERdt = ((self.SVR*np.sqrt(self.ratio))/(self.F*self.ratio)) * (-J_CER_m) - J_CERiD

        return np.hstack((dNadt, dKdt, dCdt, dIP3dt, dhdt, dCERdt, dNaodt, dKodt, dCodt))


