import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

from Astro_morphology import Astrocyte_morphology
from Astro_morphology import Process
from Astro_morphology import Connect
from Init import Init
from Parameters import p
from Finite_Difference import Finite_difference

class Astro_multi_compartment(object):

    def __init__(self, params, model_type, stimulus, morpho):

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

        # if self.end_condition == 'open_end':
        #     self.SVR = self.SVR[1:-1]
        #     self.glut_input = self.glut_input[1:-1]

        self.O_m = self.SVR * self.a_i#a_i/(5.00 * 1e-8) # meter, Astrocytic membrane area per tissue volume

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

        J_V = (2*IGluT - INKA*self.F - INCX - INleak * self.F - IKleak*self.F)

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
        #J_VD = self.conn_matrix.dot(-((1.)/(self.rL * self.length_comp**2)) * V)

        #define differential equations
        dKdt = -(self.O_m/self.a_i)*(J_K_m) - J_KiD
        dNadt = -(self.O_m/self.a_i)*(J_Na_m) - J_NaiD
        dKodt = (self.O_m/self.a_o)*(J_K_m) - J_KoD
        dNaodt = (self.O_m/self.a_o)*(J_Na_m) - J_NaoD
        dCdt = -(self.O_m/self.a_i) * J_Ca_m - J_CiD
        dCodt = (self.O_m/self.a_o) * J_Ca_m - J_CoD
        #dVdt = (1./self.C_m)*(J_V)# - J_VD)


        return np.hstack((dNadt, dKdt, dCdt, dNaodt, dKodt, dCodt))


    def spatial_astro_Ca(self, state, tspan):

        # initial values
        C = state[:self.N]
        IP3 = state[self.N:2*self.N]
        h = state[2*self.N:3*self.N]
        CER = state[3*self.N:4*self.N]

        self.O_m = self.SVR * self.a_i#a_i/(5.00 * 1e-8) # meter, Astrocytic membrane area per tissue volume


        #CER = (self.Ctot_0 - C)*(self.ratio_inv)

        # # finite difference method for calculation of pde: matrix for first derivative
        # I_1a = np.zeros((self.N,self.N))
        # for i in range(0,len(I_1a)-1):
        #     I_1a[i,i] = -1.
        #     I_1a[i,i+1] = 1.
        # I_1a /= (self.h)
        # I_1a[-1,-1] = -1./(-self.h)
        # I_1a[-1,-2] = 1./(-self.h)
        #
        # I_1b = np.zeros((self.N,self.N))
        # for i in range(0,len(I_1b)-1):
        #     I_1b[i,i] = -1.
        #     I_1b[i,i-1] = 1.
        # I_1b /= (-self.h)
        # I_1b[0,0] = -1./(self.h)
        # I_1b[0,1] = 1./(self.h)
        # I_1b[-1,-1] = -1./(-self.h)
        # I_1b[-1,-2] = 1./(-self.h)

        # membrane potential
        #X_oZ_o = -((self.O_m*self.C_m*self.V_0)/(self.a_o)) - self.F*(self.z_K*self.K_o_0 + self.z_Na*self.Na_o_0 + self.z_C*self.C_o_0)
        #V = -((self.a_o)/(self.C_m * self.O_m))*((self.F*(self.z_K * Ko + self.z_Na * Nao + self.z_C * Co)) + X_oZ_o)
        #X_oZ_o = -((self.O_m*self.C_m*self.V_0)/(self.a_o)) - self.F*(self.K_o_0 + self.Na_o_0)
        #V = -((self.a_o)/(self.C_m * self.O_m))*((self.F*(Ko + Nao)) + X_oZ_o)

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

        # diffusive flux
        J_CiD = self.conn_matrix.dot(-(self.D_C/(self.lamb_intra**2)) * C)
        J_CERiD = self.conn_matrix.dot(-(self.D_C/(self.lamb_intra**2)) * CER)
        J_IP3iD = self.conn_matrix.dot(-(self.D_IP3/(self.lamb_intra**2)) * IP3)
        # J_CD = -(self.D_C/(self.lamb_intra**2)) * I_1a.dot(C)
        # J_CERD = -(self.D_C/(self.lamb_intra**2)) * I_1a.dot(CER)
        # J_IP3D = -(self.D_IP3/(self.lamb_extra**2)) * I_1a.dot(IP3)

        # #axial flux densities
        # J_C = J_CD
        # J_CER = J_CERD
        # J_IP3 = J_IP3D
        #
        # # set boundary conditions
        # J_C[0] = J_C[-1] = 0.
        # J_CER[0] = J_CER[-1] = 0.
        # J_IP3[0] = J_IP3[-1] = 0.

        #define differential equations
        # dKdt = -(self.O_m/self.a_i)*(J_K_m) - I_1b.dot(J_Ki)
        # dNadt = -(self.O_m/self.a_i)*(J_Na_m) - I_1b.dot(J_Nai)
        # dKodt = (self.O_m/self.a_o)*(J_K_m) - I_1b.dot(J_Ko)
        # dNaodt = (self.O_m/self.a_o)*(J_Na_m) - I_1b.dot(J_Nao)
        # dCdt = -(self.O_m/self.a_i) * J_Ca_m - I_1b.dot(J_Ci)
        # dCodt = (self.O_m/self.a_o) * J_Ca_m - I_1b.dot(J_Co)
        #dCdt = ((self.SVR*np.sqrt(self.ratio))/(self.F)) * J_CER_m - I_1b.dot(J_C)
        #dIP3dt = prod_degr_IP3 - I_1b.dot(J_IP3)
        #dCERdt = ((self.SVR*np.sqrt(self.ratio))/(self.F*self.ratio)) * (-J_CER_m) - I_1b.dot(J_CER)

        dCdt = ((self.SVR*np.sqrt(self.ratio))/(self.F)) * J_CER_m - J_CiD
        dIP3dt = prod_degr_IP3 - J_IP3iD
        dCERdt = ((self.SVR*np.sqrt(self.ratio))/(self.F*self.ratio)) * (-J_CER_m) - J_CERiD

        #return np.hstack((dNadt, dKdt, dNaodt, dKodt))
        return np.hstack((dCdt, dIP3dt, dhdt, dCERdt))#, dVdt))

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

        self.O_m = self.SVR * self.a_i#a_i/(5.00 * 1e-8) # meter, Astrocytic membrane area per tissue volume


        # # finite difference method for calculation of pde: matrix for first derivative
        # self.I_1a = np.zeros((self.N,self.N))
        # for i in range(0,len(self.I_1a)-1):
        #     self.I_1a[i,i] = -1.
        #     self.I_1a[i,i+1] = 1.
        # self.I_1a[-1,-1] = 1.
        # self.I_1a[-1,-2] = -1.
        # self.I_1a /= self.h
        #
        # self.I_1b = np.zeros((self.N,self.N))
        # for i in range(0,len(self.I_1b)-1):
        #     self.I_1b[i,i] = 1.
        #     self.I_1b[i,i-1] = -1.
        # self.I_1b[0,0] = -1.
        # self.I_1b[0,1] = 1.
        # self.I_1b[-1,-1] = 1.
        # self.I_1b[-1,-2] = -1.
        # self.I_1b /= self.h

        # membrane potential
        X_oZ_o = -((self.O_m*self.C_m*self.V_0)/(self.a_o)) - self.F*(self.z_K*self.K_o_0 + self.z_Na*self.Na_o_0 + self.z_C*self.C_o_0)
        V = -((self.a_o)/(self.C_m * self.O_m))*((self.F*(self.z_K * Ko + self.z_Na * Nao + self.z_C * Co)) + X_oZ_o)
        #X_oZ_o = -((self.O_m*self.C_m*self.V_0)/(self.a_o)) - self.F*(self.K_o_0 + self.Na_o_0)
        #V = -((self.a_o)/(self.C_m * self.O_m))*((self.F*(Ko + Nao)) + X_oZ_o)

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

        # # diffusive flux
        # J_NaiD = -(self.D_Na/(self.lamb_intra**2)) * self.I_1a.dot(Na)
        # J_KiD = -(self.D_K/(self.lamb_intra**2)) * self.I_1a.dot(K)
        # J_CiD = -(self.D_C/(self.lamb_intra**2)) * self.I_1a.dot(C)
        #
        # J_CERD = -(self.D_C/(self.lamb_intra**2)) * self.I_1a.dot(CER)
        # J_IP3D = -(self.D_IP3/(self.lamb_extra**2)) * self.I_1a.dot(IP3)
        #
        # J_NaoD = -(self.D_Na/(self.lamb_extra**2)) * self.I_1a.dot(Nao)
        # J_KoD = -(self.D_K/(self.lamb_extra**2)) * self.I_1a.dot(Ko)
        # J_CoD = -(self.D_C/(self.lamb_extra**2)) * self.I_1a.dot(Co)

        # diffusive flux
        J_CiD = self.conn_matrix.dot(-(self.D_C/(self.lamb_intra**2)) * C)
        J_CoD = self.conn_matrix.dot(-(self.D_C/(self.lamb_extra**2)) * Co)
        J_KiD = self.conn_matrix.dot(-(self.D_K/(self.lamb_intra**2)) * K)
        J_NaiD = self.conn_matrix.dot(-(self.D_Na/(self.lamb_intra**2)) * Na)
        J_KoD = self.conn_matrix.dot(-(self.D_K/(self.lamb_extra**2)) * Ko)
        J_NaoD = self.conn_matrix.dot(-(self.D_Na/(self.lamb_extra**2)) * Nao)
        J_CERiD = self.conn_matrix.dot(-(self.D_C/(self.lamb_intra**2)) * CER)
        J_IP3iD = self.conn_matrix.dot(-(self.D_IP3/(self.lamb_intra**2)) * IP3)




        # # intra- and extracellular resistivity
        # r_o = (self.psifac*(self.lamb_extra**2))/(self.F*(self.D_Na*(self.z_Na**2)*Nao+self.D_K*(self.z_K**2)*Ko+self.D_C*(self.z_C**2)*Co))
        # r_i = (self.psifac*(self.lamb_intra**2))/(self.F*(self.D_Na*(self.z_Na**2)*Na+self.D_K*(self.z_K**2)*K+self.D_C*(self.z_C**2)*C))
        # #r_o = (self.psifac*(self.lamb_extra**2))/(self.F*(self.D_Na*Nao+self.D_K*Ko))
        # #r_i = (self.psifac*(self.lamb_intra**2))/(self.F*(self.D_Na*Na+self.D_K*K))
        #
        # # current densities due to diffusion
        # #i_odiff = self.F*(self.z_K*J_KoD + self.z_Na*J_NaoD)
        # #i_idiff = self.F*(self.z_K*J_KiD + self.z_Na*J_NaiD)
        # i_odiff = self.F*(self.z_K*J_KoD + self.z_Na*J_NaoD + self.z_C*J_CoD)
        # i_idiff = self.F*(self.z_K*J_KiD + self.z_Na*J_NaiD + self.z_C*J_CiD)
        #
        # # calculate intra- and extracellular membrane voltage
        # dVidx = (self.I_1a.dot(V) + ((r_o*self.a_i*i_idiff)/self.a_o) + (r_o*i_odiff))*((1. + ((r_o*self.a_i)/(r_i*self.a_o)))**-1)
        # dVodx = (-self.I_1a.dot(V) + ((r_i*self.a_o*i_odiff)/self.a_i) + (r_i*i_idiff))*((1. + ((r_i*self.a_o)/(r_o*self.a_i)))**-1)
        #
        # # field flux
        # J_KiV = -((self.D_K*self.z_K)/((self.lamb_intra**2) * self.psifac)) * (K*dVidx)
        # J_NaiV = -((self.D_Na*self.z_Na)/((self.lamb_intra**2) * self.psifac)) * (Na*dVidx)
        # J_KoV = -((self.D_K*self.z_K)/((self.lamb_extra**2) * self.psifac)) * (Ko*dVodx)
        # J_NaoV = -((self.D_Na*self.z_Na)/((self.lamb_extra**2) * self.psifac)) * (Nao*dVodx)
        # J_CiV = -((self.D_C*self.z_C)/((self.lamb_intra**2) * self.psifac)) * (C*dVidx)
        # J_CoV = -((self.D_C*self.z_C)/((self.lamb_extra**2) * self.psifac)) * (Co*dVodx)

        # #axial flux densities
        # J_Nai = J_NaiD #+ J_NaiV
        # J_Ki = J_KiD #+ J_KiV
        # J_Ci = J_CiD #+ J_CiV
        #
        # J_IP3 = J_IP3D
        # J_CER = J_CERD
        #
        # J_Nao = J_NaoD #+ J_NaoV
        # J_Ko = J_KoD #+ J_KoV
        # J_Co = J_CoD #+ J_CoV
        #
        # # set boundary conditions: sealed end
        # J_Nai[0] = J_Nai[-1] = 0.
        # J_Ki[0] = J_Ki[-1] = 0.
        # J_Ci[0] = J_Ci[-1] = 0.
        #
        # J_IP3[0] = J_IP3[-1] = 0.
        # J_CER[0] = J_CER[-1] = 0.
        #
        # J_Nao[0] = J_Nao[-1] = 0.
        # J_Ko[0] = J_Ko[-1] = 0.
        # J_Co[0] = J_Co[-1] = 0.

        # external bath
        #bath_K = - self.k_dec*(Ko - self.K_o_0)
        #bath_Na = - self.k_dec*(Nao - self.Na_o_0)

        # input
        #self.Imax = self.input(tspan)

        # #define differential equations
        # dNadt = -(self.O_m/self.a_i)*(J_Na_m) - self.I_1b.dot(J_Nai)
        # dKdt = -(self.O_m/self.a_i)*(J_K_m) - self.I_1b.dot(J_Ki)
        # dCdt = -(self.O_m/self.a_i) * J_Ca_m + ((self.SVR*np.sqrt(self.ratio))/(self.F)) * J_CER_m - self.I_1b.dot(J_Ci)
        #
        # dIP3dt = prod_degr_IP3 - self.I_1b.dot(J_IP3)
        # dCERdt = ((self.SVR*np.sqrt(self.ratio))/(self.F*self.ratio)) * (-J_CER_m) - self.I_1b.dot(J_CER)
        #
        # dNaodt = (self.O_m/self.a_o)*(J_Na_m) - self.I_1b.dot(J_Nao)
        # dKodt = (self.O_m/self.a_o)*(J_K_m) - self.I_1b.dot(J_Ko)
        # dCodt = (self.O_m/self.a_o) * J_Ca_m - self.I_1b.dot(J_Co)

        #define differential equations
        dKdt = -(self.O_m/self.a_i)*(J_K_m) - J_KiD
        dNadt = -(self.O_m/self.a_i)*(J_Na_m) - J_NaiD
        dKodt = (self.O_m/self.a_o)*(J_K_m) - J_KoD
        dNaodt = (self.O_m/self.a_o)*(J_Na_m) - J_NaoD
        dCdt = -(self.O_m/self.a_i) * J_Ca_m + ((self.SVR*np.sqrt(self.ratio))/(self.F)) * J_CER_m - J_CiD
        dCodt = (self.O_m/self.a_o) * J_Ca_m - J_CoD
        dIP3dt = prod_degr_IP3 - J_IP3iD
        dCERdt = ((self.SVR*np.sqrt(self.ratio))/(self.F*self.ratio)) * (-J_CER_m) - J_CERiD

        #return np.hstack((dNadt, dKdt, dNaodt, dKodt))
        return np.hstack((dNadt, dKdt, dCdt, dIP3dt, dhdt, dCERdt, dNaodt, dKodt, dCodt))#, dVdt))


