# Yusuf's qutip simulator. Thanks Yusuf!
import numpy as np
import qutip as qt
import matplotlib.pyplot as plt
# import plotly.express as px
# import plotly.graph_objects as go
# import plotly.io as pio
import configparser

# Define a class named Pulse to create laser pulses 
class Pulse:
    def __init__(self, tau, e_start, w_gain=0, t0=0, e0=1, phase=0, polar_x=1):
        self.tau = tau  # in ps (tau is a time constant)
        self.e_start = e_start  # in meV (e_start is the starting energy)
        self.w_start = e_start / HBAR  # in 1 / ps (calculate angular frequency)
        self.w_gain = float(w_gain)  # in 1/ps^2 (rate of change of angular frequency)
        self.t0 = t0  # initial time
        self.e0 = e0  # initial energy
        self.phase = phase  # phase of the pulse
        self.freq = None  # frequency , initialized to None
        self.phase_ = None  # phase_ , initialized to None
        self.polar_x = polar_x  # x-component of polarization
        self.polar_y = np.sqrt(1 - polar_x ** 2)  # calculate y-component of polarization

    def __repr__(self):
        
        return "%s(tau=%r, e_start=%r, w_gain=%r, t0=%r, e0=%r)" % (
            self.__class__.__name__, self.tau, self.e_start, self.w_gain, self.t0, self.e0
        )

    def get_envelope(self, t):
        # Calculate the envelope of the pulse at a given time t
        return self.e0 * np.exp(-0.5 * ((t - self.t0) / self.tau) ** 2) / (np.sqrt(2 * np.pi) * self.tau)
    
    def set_frequency(self, f):
        #Set the time-dependent frequency using a lambda function f that takes time t.
        
        self.freq = f

    def get_frequency(self, t):
        #Calculate the frequency (omega) for a given time t.
        
        if self.freq is not None:
            return self.freq(t)
        return self.w_start + self.w_gain * (t - self.t0)

    def get_full_phase(self, t):
        # Calculate the full phase of the pulse at a given time t
        return self.w_start * (t - self.t0) + 0.5 * self.w_gain * ((t - self.t0) ** 2) + self.phase
    
    def get_energies(self):
        #Get the energy difference of +- tau for a chirped pulse.
        
        low = self.get_frequency(-self.tau)
        high = self.get_frequency(self.tau)
        energy_range = np.abs(high - low) * HBAR  # Calculate energy range in meV
        return energy_range

    def get_total(self, t):
        # Calculate pulse envelope at a given time t
        return self.get_envelope(t) * np.exp(-1j * self.get_full_phase(t))


# Chirping the laser pulses
class ChirpedPulse(Pulse):
    def __init__(self, tau_0, e_start, alpha=0, t0=0, e0=1 * np.pi, polar_x=1):
        self.tau_0 = tau_0  # Initial tau value
        self.alpha = alpha  # Chirp Value
        # Initialize the ChirpedPulse using parameters and super() to call the parent class's constructor
        super().__init__(tau=np.sqrt(alpha**2 / tau_0**2 + tau_0**2), e_start=e_start, w_gain=alpha / (alpha**2 + tau_0**4), t0=t0, e0=e0, polar_x=polar_x)
    
    def get_parameters(self):
        """
        Returns tau and chirp parameter.
        """
        return "tau: {:.4f} ps , a: {:.4f} ps^-2".format(self.tau, self.w_gain)

    def get_envelope(self, t):
        # Calculate the pulse envelope  at a given time t
        return self.e0 * np.exp(-0.5 * ((t - self.t0) / self.tau) ** 2) / (np.sqrt(2 * np.pi * self.tau * self.tau_0))

    def get_ratio(self):
        """
        Returns the ratio of pulse area chirped/unchirped: tau / sqrt(tau * tau_0).
        """
        return np.sqrt(self.tau / self.tau_0)



# Biexciton Function - Run this function
HBAR = 0.6582173  # meV*ps

# delta_b: binding energy between exciton and biexciton
def energies(delta_b=4., delta_0=0.):
    # energy levels of the system
    E_X = -delta_0/2
    E_Y =  delta_0/2
    E_B = -delta_b
    return E_X, E_Y, E_B

def read_calibration_file(calibration_file):

    # reads in experimentally aquired quantum dot parameters 
    config = configparser.ConfigParser()
    config.read(calibration_file)

    # read the calibration file
    central_wavelength = float(config['EMISSION']['exciton_wavelength']) #nm
    biexciton_wavelength = float(config['EMISSION']['biexciton_wavelength'])
    dark_wavelength = float(config['EMISSION']['dark_wavelength']) 

    fss_bright = float(config['SPLITTING']['fss_bright'])*1e-3 #meV
    fss_dark = float(config['SPLITTING']['fss_dark']) *1e-3 # meV 

    lifetime_exciton = float(config['LIFETIMES']['exciton']) #ps
    lifetime_biexciton = float(config['LIFETIMES']['biexciton'])
    #lifetime_dark = float(config['LIFETIMES']['dark']) 

    g_ex = float(config['G_FACTORS']['g_ex'])
    g_hx = float(config['G_FACTORS']['g_hx'])
    g_ez = float(config['G_FACTORS']['g_ez'])
    g_hz = float(config['G_FACTORS']['g_hz'])

    exciton_meV = 1239.8*1e3/central_wavelength #meV
    biexciton_meV = 1239.8*1e3/biexciton_wavelength
    dark_meV = 1239.8*1e3/dark_wavelength

    exciton_x_energy = fss_bright/2
    exciton_y_energy = -fss_bright/2
    binding_energy = -(exciton_meV - biexciton_meV) # negatively defined
    dark_energy = (dark_meV-exciton_meV)
    dark_x_energy = dark_energy + fss_dark/2
    dark_y_energy = dark_energy - fss_dark/2 

    gamma_e = 1/lifetime_exciton
    gamma_b = 1/(lifetime_biexciton*2)
    #gamma_d = 1/lifetime_dark

    return exciton_x_energy, exciton_y_energy, dark_x_energy, dark_y_energy, binding_energy, gamma_e, gamma_b, g_ex, g_hx, g_ez, g_hz

def fourlevel_system(collapse= False, tau1=3, tau2=3, area1=1*np.pi, area2=0, det1=0, det2=0, alpha1=0, alpha2=0, prob_b=0,
    pol1_x=1, pol2_x=1, delay=0, delta_b=4, delta_0=0.0, gamma_e=1/100, gamma_b=1/100, epsilon=0.01, timeAxis = "auto", t_user=None,timeAxis_smart = False,
    dt_1=0.1, dt_2=0.1, options=qt.Options(atol=1e-7), mode="population", calibration_file=None, pulse_x = None, pulse_y = None,t0 = 0, tend = 200, dt=0.01):
    """
    In qutip, every energy has to be provided in 1/ps
    Here, a rotating frame with the unsplit exciton energy is chosen. 
    collapse: choose "decay" if you wanna include decay otherwise choose "nodecay" to leave collapse operators empty. 
    tau1/2: pulse 1/2 duration in ps
    area1/2: pulsearea of pulse 1/2
    det1/2: detuning of pulse 1/2 to unsplit exciton energy in meV
    alpha1/2: chirp of pulse1/2 in ps^2
    pol1/2_x: x polarization component of pulse 1/2. possible options = 0,...,1
    delay: delay of pulse 2 to pulse 1 in ps
    delta_b: biexciton binding in meV
    delta_0: exciton X/Y splitting in meV (FSS)
    gamma_e: inverse exciton lifetime in 1/ps
    gamma_b: inverse biexciton lifetime in 1/ps
    epsilon: exponential decay, until epsilon is reached
    dt_1: timestep during pulse (0,..,8tau)
    dt_2: timestep after the pulse, during the decay
    mode: "pop" for population. This is for various possible modes in the future. 
    """
    
    delta_b = delta_b / HBAR  # delta_b in 1/ps
    delta_0 = delta_0 / HBAR

    # system states
    g =qt.basis(4,0)
    x = qt.basis(4,1)
    y = qt.basis(4,2)
    b = qt.basis(4,3)
    
    # Initialize the ground state. For example, if you want to start from the ground state, set prob_b to zero.
    # "prob_b" represents the probability amplitude of the biexciton state. If you want to start with a superposition of the ground state and exciton, then change "np.sqrt(prob_b) * b" to "x" or "y".
    # If you want to have a superposition of the ground state, exciton, and biexciton, then you can change "gxbas" as you wish.
    
    gxbas =np.sqrt(1-prob_b)*g + np.sqrt(prob_b)*b

    # number operators
    n_g = g * g.dag()
    n_x = x * x.dag()
    n_y = y * y.dag()
    n_b = b * b.dag()

    # transition operators / polarizations
    p_gx = g * x.dag()
    p_gy = g * y.dag()
    p_xb = x * b.dag()
    p_yb = y * b.dag()
    # this one is not needed for the hamiltonian
    p_gb = g * b.dag()

    # Check the value of the collapse parameter
    
        
    # Read values from calibration file
    if calibration_file is not None:
        E_X, E_Y, _, _, E_B, gamma_e, gamma_b, _, _, _, _ = read_calibration_file(calibration_file)
        E_X = E_X / HBAR
        E_Y = E_Y / HBAR
        E_B = E_B / HBAR
    else:
        E_X, E_Y, E_B = energies(delta_b=delta_b, delta_0=delta_0) # note they are already divided by HBAR
    # Check the value of the collapse parameter
   
    if not collapse:
        gamma_b = 0
        gamma_e = 0 

    c_ops = [np.sqrt(gamma_e) * p_gx, np.sqrt(gamma_e) * p_gy, np.sqrt(gamma_b) * p_xb, np.sqrt(gamma_b) * p_yb]

    if not collapse:
        c_ops = []
    # system Hamiltonian
    H_sys = E_X * n_x + E_Y * n_y + E_B * n_b
    #print(H_sys)

    # pulse 1 and 2, right now assume delay > 0
    tau11=np.sqrt(alpha1**2 / tau1**2 + tau1**2)
    tau22=np.sqrt(alpha2**2 / tau2**2 + tau2**2)
    # choose the longer of the two
    t_start1 = 4*tau11 #if tau11 > tau22 else 4*tau22
    # further delay pulse 2
    t_start2 = t_start1 + delay
    pulse1 = ChirpedPulse(tau1, det1, alpha1, t0=t_start1, e0=area1, polar_x=pol1_x)
    pulse2 = ChirpedPulse(tau2, det2, alpha2, t0=t_start2, e0=area2, polar_x=pol2_x)

    # excitation Hamiltonians (daggered, as expressed by polarization operators)
    H_x_dag = -0.5 * (p_gx + p_xb)  # this has to be paired with the conjugated total x-field 
    H_y_dag = -0.5 * (p_gy + p_yb)
    # print(H_x_dag.dag())
    if (pulse_x is None) or (pulse_y is None):
        H = [H_sys, [H_x_dag, lambda t,args : np.conj(pulse1.polar_x*pulse1.get_total(t) + pulse2.polar_x*pulse2.get_total(t))],
                    [H_y_dag, lambda t,args : np.conj(pulse1.polar_y*pulse1.get_total(t) + pulse2.polar_y*pulse2.get_total(t))],
                    [H_x_dag.dag(), lambda t,args : pulse1.polar_x*pulse1.get_total(t) + pulse2.polar_x*pulse2.get_total(t)],
                    [H_y_dag.dag(), lambda t,args : pulse1.polar_y*pulse1.get_total(t) + pulse2.polar_y*pulse2.get_total(t)]]
    else:
        H = [H_sys, [H_x_dag, lambda t,args : np.conj(pulse_x(t))],
                    [H_y_dag, lambda t,args : np.conj(pulse_y(t))],
                    [H_x_dag.dag(), lambda t,args : pulse_x(t)],
                    [H_y_dag.dag(), lambda t,args : pulse_y(t)]]
        
        # def pulse_x_func(t,args):
        #     return pulse_x(t).astype(complex)
        
        # def pulse_x_func_conj(t,args):
        #     return np.conj(pulse_x(t)).astype(complex)
        
        # def pulse_y_func(t,args):
        #     return pulse_y(t).astype(complex)
        
        # def pulse_y_func_conj(t,args):
        #     return np.conj(pulse_y(t)).astype(complex)    
        
        
        # H = qt.Qobj([H_sys, [H_x_dag*qt.coefficient(pulse_x_func_conj)],
        #             [H_y_dag*qt.coefficient(pulse_y_func_conj)],
        #             [H_x_dag.dag()*qt.coefficient(pulse_x_func)],
        #             [H_y_dag.dag()*qt.coefficient(pulse_y_func)]])
        
        # H = H_sys + H_x_dag*qt.coefficient(pulse_x_func_conj) + H_y_dag*qt.coefficient(pulse_y_func_conj) + H_x_dag.dag()*qt.coefficient(pulse_x_func) + H_y_dag.dag()*qt.coefficient(pulse_y_func)
        
        # H = qt.QobjEvo([H_sys, [H_x_dag, lambda t,args : np.conj(pulse_x(t))],
        #              [H_y_dag, lambda t,args : np.conj(pulse_y(t))],
        #              [H_x_dag.dag(), lambda t,args : pulse_x(t)],
        #              [H_y_dag.dag(), lambda t,args : pulse_y(t)]])
    
    

    if timeAxis == "manual":
        t_axis= np.arange(0, t_user, dt_1)

    if timeAxis == "auto":

        # time axes. has to start at 0 due to limitations in the function calculating the 2-time quantities
        # two different time steps are used: a small dt_1, during the time the pulses are active, and a larger dt_2 during the decay 
        # time axis during the pulses
        t_off = t_start2 + t_start1  # time window where pulse 1 or 2 is still active
        rate = 2*gamma_b if 2*gamma_b<gamma_e else gamma_e
        t_end = t_off - 1/rate *np.log(epsilon)  # note that log(epsilon) is in general negative
        t_axis1 = np.arange(0, t_off, dt_1)
        t_axis2 = np.arange(t_off, t_end, dt_2)  # note that arange does not include the final value so t_off is not in both arrays
        t_axis = np.append(t_axis1, t_axis2)

    if timeAxis.lower()[0] == 'p': 
        #for useage with PULSE
        if not timeAxis_smart:
            t_axis = np.arange(t0,tend,dt)
        else:
            t_axis = [t0]
            while t_axis[-1] < tend:
                cur_field = np.abs(pulse_x(t_axis[-1])) + np.abs(pulse_y(t_axis[-1]))
                cur_dt = np.min([3, 0.1 * 1/cur_field])
                cur_dt = np.max([dt, cur_dt])
                t_axis.append(t_axis[-1] + cur_dt)   
            t_axis[-1] = tend
            t_axis = np.array(t_axis)
            #print(t_axis)
            
    if mode == "population":
        g_occ, x_occ, y_occ, b_occ, polar_gx, polar_xb, polar_gb = qt.mesolve(H, gxbas, t_axis, c_ops=c_ops, e_ops=[n_g, n_x, n_y, n_b, p_gx, p_xb, p_gb], options=options).expect
        
        # g_occ, x_occ, y_occ, b_occ, polar_gx, polar_xb, polar_gb = qt.mcsolve(H, gxbas, t_axis, c_ops=c_ops, e_ops=[n_g, n_x, n_y, n_b, p_gx, p_xb, p_gb], options=options, ntraj=5).expect
        
        return t_axis, g_occ, x_occ, y_occ, b_occ, polar_gx, polar_xb, polar_gb, t_axis, pulse1, pulse2


# gamma_e=1/100 ; gamma_b=1/100 ; 
# g, x,y,b,gx,xb,gb,t, pulse1,pulse2 = fourlevel_system(collapse="nodecay", tau1=2, dt_1=0.1,dt_2=0.1, prob_b=0,
#             gamma_b=gamma_b, gamma_e=gamma_e, delta_b=4, timeAxis= "manual", t_user=150,
#             area1=4*np.pi, alpha1=0, det1=-2, pol1_x=np.sqrt(0), pol2_x=np.sqrt(0), 
#             tau2=2 ,alpha2=0, area2=0.5*np.pi, det2=-4, delay=15, mode="population"
#             ) 