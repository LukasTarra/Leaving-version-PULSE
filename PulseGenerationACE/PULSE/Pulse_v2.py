import numpy as np
from matplotlib import pyplot as plt
from scipy.io import savemat
import pyaceqd.pulsegenerator as pg
from pyaceqd.six_level_system.linear import sixls_linear
from pyaceqd.six_level_system.linear import energies_linear
from pyaceqd.six_level_system.linear import sixls_linear_dressed_states 
from pyaceqd.four_level_system.linear import biexciton
from pyaceqd.four_level_system.linear import biexciton_dressed_states
from QUtip_four_level_system import fourlevel_system 
from IPython.display import clear_output
import time as time 
from uncertainties import ufloat
from uncertainties import *
from uncertainties import unumpy as unp
from scipy.io import loadmat
from pyaceqd.tools import read_calibration_file
import configparser  
from scipy.optimize import curve_fit
from datetime import datetime
from scipy.optimize import minimize
from scipy.optimize import fsolve
import pickle
import copy
import os
from time import sleep
import csv
import qutip as qt
from pyaceqd.tools import read_calibration_file

import pulse_shaper_control as psc
import spectrometer_control as smc
import attenuator_control as atc
import simulator_control as simc
import time_delay_control as tdc
import powermeter_control as pmc
import wave_plate_control as wpc

HBAR = 0.6582173666

class create_experiment():
    def __init__(self,path = None, name='Experiment', save = False) -> None:
        # check if path exists
        if path is not None:
            self.exp_path = path + '/' + name
        else:
            # set current working directory
            self.exp_path = os.getcwd() + '/' + name
        
        #create folder
        self.pulse_path = self.exp_path + '/Pulse_files'
        self.device_path = self.exp_path + '/Device_files'
        self.temp_path = self.exp_path + '/Temp_files'

        os.chdir(self.exp_path)

        os.makedirs(self.exp_path,exist_ok=True)
        os.makedirs(self.pulse_path,exist_ok=True)
        os.makedirs(self.device_path,exist_ok=True)
        os.makedirs(self.temp_path,exist_ok=True)

        self.name = name
        self.time_of_creation = datetime.now()
        self.suffix = 'ex'
        self.print_info()
        if save:
            save_device(self, self.name,suffix = self.suffix)
        pass

    def print_info(self):
        print('\nExperiment created!')
        print('Experiment name:',self.name)
        print('Experiment time of creation:',self.time_of_creation)
        pass

class time_delay():
    def __init__(self,device = None,name='Time Delay',save = False,experiment = None,time_delay_sleep = 0, scale = 1, offset = 0, num_passes = 1, unit = 'mm') -> None:
        self.name = name
        self.device = device
        self.time_of_creation = datetime.now()
        self.suffix = 'td'
        self.time_delay_sleep = time_delay_sleep
        self.print_info()
        _check_method(device,'set_position')
        _check_method(device,'close')
        
        self.scale = scale # in mm 
        self.offset = offset # in ps 
        self.num_passes = num_passes # number of passes, e.g. backrefector = 2    
        self.unit = unit

        if save:
            save_device(self, self.name,suffix = self.suffix,experiment = experiment)
    
    def convert_ps_to_mm(self):
        c = 0.299792 # mm per ps
        self.motor_delay =c*(self.delay/self.num_passes -self.offset)
        self.motor_delay = self.motor_delay/self.scale
    
    def get_offset(self):
        return self.offset
    
    def excecute(self, excecute = True):
        # meant to be called from saved object and excecute device specific functions
        self.convert_ps_to_mm()
        self.device.set_position(self.motor_delay, excecute = excecute)
        sleep(self.time_delay_sleep)

    def set_delay(self,delay=None,pulse_object = None,excecute=False):
            
            if delay is not None:
                self.delay = delay
    
            if excecute:
                if self.unit == 'mm':
                    self.convert_ps_to_mm()
                elif self.unit == 'ps':
                    self.motor_delay = self.delay/self.scale - self.offset
                self.device.set_position(self.motor_delay, excecute = excecute) 
    
            if pulse_object is not None:
                if type(pulse_object) is str:
                    pulse_object = pg.load_pulse(pulse_object)
                out_pulse_object = pulse_object.copy_pulse()
                out_pulse_object = pulse_time_delay(out_pulse_object,delay)
                return out_pulse_object

    def print_info(self):
        print('\nTime Delay connected!')
        print('Time Delay name:',self.name)
        print('Time Delay device:',self.device)
        print('Time Delay time of creation:',self.time_of_creation)
        pass
    
    def open_control(self,pulse_object = None,open_gui = True,parent_window = None,previous_control = None):
       
        control_object = tdc.time_delay_control(self,pulse_object = pulse_object, open_gui = open_gui,parent_window=parent_window,previous_control= previous_control)
        return control_object
    
    def close(self):
        self.device.close()

class attenuator():
    def __init__(self,device=None,name='Attenuator',save = False,experiment = None,attenuator_sleep = 0, command = 'set_attenuation') -> None:
        self.name = name
        self.device = device
        self.time_of_creation = datetime.now()
        self.suffix = 'at'
        self.attenuator_sleep = attenuator_sleep
        self.current_attenuation = 0
        self.command = command
        
        self.print_info()
        _check_method(device,self.command)
        _check_method(device,'close')

        if save:
            save_device(self, self.name,suffix = 'at',experiment = experiment)
    
    def device_command(self,position):
        com = getattr(self.device,self.command)
        com(position)
    
    def excecute(self):
        # meant to be called from saved object and excecute device specific functions
        #self.device.set_attenuation(self.attenuation)
        self.device_command(self.attenuation)
        sleep(self.attenuator_sleep)

    def set_attenuation(self,attenuation=None,pulse_object = None,excecute=False):
        
        if attenuation is not None:
            self.attenuation = attenuation

        if excecute and attenuation != self.current_attenuation:
            self.current_attenuation = attenuation
            #self.device.set_attenuation(attenuation) 
            self.device_command(attenuation)

        if pulse_object is not None:
            if type(pulse_object) is str:
                 pulse_object = pg.load_pulse(pulse_object)
            out_pulse_object = pulse_object.copy_pulse()
            out_pulse_object.clear_filter()
            out_pulse_object.add_filter_rectangle(transmission = np.sqrt(attenuation),cap_transmission=False)
            out_pulse_object.apply_frequency_filter()
            out_pulse_object.clear_filter()
            return out_pulse_object 
    
    def print_info(self):
        print('\nAttenuator connected!')
        print('Attenuator name:',self.name)
        print('Attenuator device:',self.device)
        print('Attenuator time of creation:',self.time_of_creation)
        pass
    
    def open_control(self,pulse_object = None,open_gui = True,parent_window = None,previous_control = None):
       
        control_object = atc.attenuator_control(self,pulse_object = pulse_object, open_gui = open_gui,parent_window=parent_window,previous_control= previous_control)
        return control_object

    def close(self):
        self.device.close()
        

class motor():
    def __init__(self,device=None,name='Motor',save = False,experiment = None, motor_sleep = 0, min_position = None, max_position = None, command = 'set_position') -> None:
        self.name = name
        self.device = device
        self.time_of_creation = datetime.now()
        self.suffix = 'mo'
        self.motor_sleep = motor_sleep
        self.current_position = 0
        self.command = command
   
        self.print_info()
        _check_method(device,self.command)
        _check_method(device,'close')
        
        if save:
            save_device(self, self.name,suffix = self.suffix,experiment = experiment)
        pass
    
        if min_position is None:
            self.min_position = -np.inf 
        else:
            self.min_position = min_position
            
        if max_position is None:
            self.max_position = np.inf
        else:
            self.max_position = max_position

    def device_command(self,position):
        com = getattr(self.device,self.command)
        com(position)
    
    def set_position(self,position=0,excecute=False):
        
        if position < self.min_position:
            position = self.min_position
        elif position > self.max_position:
            position = self.max_position
        
        if excecute and position != self.current_position:
            self.current_position = position
            #self.device.set_position(position)
            self.device_command(position)
            sleep(self.motor_sleep)
        pass

    def print_info(self):
        print('\nMotor connected!')
        print('Motor name:',self.name)
        print('Motor device:',self.device)
        print('Motor time of creation:',self.time_of_creation)
        pass

    def close(self):
        self.device.close()

# measurement devices 
class power_meter():
    def __init__(self,device=None,name='Power Meter',save = False,experiment = None, command = 'get_power') -> None:
        self.name = name
        self.device = device
        self.time_of_creation = datetime.now()
        self.suffix = 'pm'
        self.print_info()
        self.command = command  
        _check_method(device,self.command)
        _check_method(device,'close')
        self.output = [None, # experiment_power
                       None] # pulse_power
        if save:
            save_device(self, self.name,suffix = self.suffix,experiment = experiment)
        pass
    
    def device_command(self):
        com = getattr(self.device,self.command)
        return com()
    
    def get_power(self,excecute=True, pulse_object = None):
        if excecute:
            if self.name == 'fake':
                self.output[0] = self.device.get_power(pulse_object)
            else:
                #self.output[0] = self.device.get_power()
                self.output[0] = self.device_command()
        
        if pulse_object is not None:
            if type(pulse_object) is str:
                 pulse_object = pg.load_pulse(pulse_object)
            self.output[1] = pulse_object.pulse_power
        return self.output
        
    def print_info(self):
        print('\nPower Meter connected!')
        print('Power Meter name:',self.name)
        print('Power Meter device:',self.device)
        print('Power Meter time of creation:',self.time_of_creation)
        pass
    
    def open_control(self, pulse_object = None,previous_control=[],open_gui = True):
        return pmc.powermeter_control(self,pulse_object,previous_control, open_gui)
    
    def close(self):
        self.device.close()


class spectrometer():
    def __init__(self,device=None,name='Spectrometer',save = False, use_sim_pulse = False, polarisation_angle = None, command = 'get_spectrum') -> None:
        self.name = name
        self.device = device
        self.time_of_creation = datetime.now()
        self.suffix = 'sp'
        self.command = command
        self.print_info()
        _check_method(device,self.command)
        _check_method(device,'close')
        self.output = [[None,None], #exp_wl, exp_counts
                       [None,None], # pul_wl, pul_counts
                        [None,None]] # sim_wl, sim_counts
        #exp_wl, _ =  self.device.get_spectrum()
        exp_wl, _ = self.device_command()
        self.exp_wl_0 = exp_wl
        self.use_sim_pulse = use_sim_pulse
        self.central_wavelength = np.mean(exp_wl)
        self.polarisation_angle = polarisation_angle
        self.lp = linear_polarizer()
        
        if save:
            save_device(self, self.name,suffix = self.suffix)
        pass
    
    def device_command(self):
        com = getattr(self.device,self.command)
        wl, counts = com()
        return wl, counts
    
    def set_polarisation(self,polarisation):
        self.polarisation_angle = polarisation
    
    def get_central_wavelength(self):
        return self.central_wavelength
    
    def get_spectrum(self,pulse_object= None,sim_object = None,excecute=True):
        if excecute:
            #exp_wl, exp_counts =  self.device.get_spectrum()
            exp_wl, exp_counts = self.device_command()
            self.output[0] = [exp_wl, exp_counts]

        else:
            exp_wl = self.exp_wl_0

        if pulse_object is not None:
            if type(pulse_object) is str:
                pulse_object = pg.load_pulse(pulse_object)
            
            if self.polarisation_angle is not None:
                self.lp.set_angle(self.polarisation_angle) 
                pulse_object = self.lp.rotate(pulse_object=pulse_object.copy_pulse())    
            
            pul_counts = _interpolate_pulse_to_spectrum(pulse_object,exp_wl,intensity=True)
            
            self.output[1] = [exp_wl, pul_counts]

        if sim_object is not None:
            if type(sim_object) is not list:
                sim_object = [sim_object]
                
            if not self.use_sim_pulse:
                sim_counts = np.zeros_like(exp_wl)
                for i in range(len(sim_object)):
                    for i_sim, sim_wl in enumerate(sim_object[i].photon_wavelength):
                        if self.polarisation_angle is None:
                            weight = 1
                        else:
                            #print(sim_object[i].photon_polarisation[i_sim])
                            cur_pol = self.lp.rotate(pol_vector=sim_object[i].photon_polarisation[i_sim])
                            weight = np.abs(cur_pol[0])**2 + np.abs(cur_pol[1])**2
                        sim_counts[int(np.argmin(np.abs(exp_wl-sim_wl)))] += sim_object[i].photon_emission[i_sim]*weight
                self.output[2] = [exp_wl, sim_counts]

            else: 
                sim_pulse = sim_object[0].photon_pulse.copy_pulse()
                sim_pulse.clear_filter()
                if self.polarisation_angle is not None:
                    self.lp.set_angle(self.polarisation_angle) 
                    sim_pulse = self.lp.rotate(pulse_object=sim_pulse)
                for i in range(len(sim_object)-1):
                    sim_pulse.merge_pulses(sim_object[i+1].photon_pulse)
                
                sim_counts = _interpolate_pulse_to_spectrum(sim_pulse,exp_wl,intensity=True)
                self.output[2] = [exp_wl, sim_counts]
        return self.output
    

    def print_info(self):
        print('\nSpectrometer connected!')
        print('Spectrometer name:',self.name)
        print('Spectrometer device:',self.device)
        print('Spectrometer time of creation:',self.time_of_creation)
        pass

    def get_slice(self,slice_center = [], slice_width = [], excecute=False, output_kind = None):
        # output kind : which alreade recorded spectrum to use
        # 0 -> lab, 1-> pulse, 2-> sim
        if excecute:
            wavelength, counts = self.get_spectrum()[0]
        else:
            wavelength, counts = self.output[output_kind]
    
        slice_counts_vec = []
        slice_center_vec = []
        for i_slice, slice_c in enumerate(slice_center):
            
            i_slice_width = np.min([i_slice,len(slice_width)-1])
            L_slice = np.where(np.logical_and(wavelength >= slice_c-slice_width[i_slice_width]/2,wavelength <= slice_c+slice_width[i_slice_width]/2))
            slice_counts_vec.append(np.sum(counts[L_slice]))
            slice_center_vec.append(np.mean(wavelength[L_slice]))
            print('Slice center:',slice_center_vec[i_slice])
            print('Slice counts:',slice_counts_vec[i_slice])
        return slice_center_vec, slice_counts_vec

    def open_control(self,pulse_object = None,simulation_control = None,open_gui = True,parent_window = None,previous_control = None):
        
        control_object = smc.spectrometer_control(self,pulse_object = pulse_object,simulation_control= simulation_control, open_gui = open_gui,parent_window=parent_window, previous_control = previous_control)
        return control_object
    
    def close(self):
        self.device.close()

    
class pulse_shaper_obj():
    def __init__(self,device=None,name='Pulse Shaper',calibration_file = None, save = False,motor_sleep = 0) -> None:
        #device is a PULSE motor object
        self.name = name
        self.device = device
        self.time_of_creation = datetime.now()
        self.suffix = 'ps'
        self.calibration_file = calibration_file
        self.motor_sleep = motor_sleep
        self.slit_center_wavelength = None
        self.print_info()
        _check_method(device,'set_position')
        _check_method(device,'close')

        self.set_calibration_file(self.calibration_file)

        if save:
            save_device(self, self.name,suffix = self.suffix)
        pass

    def print_info(self):
        print('\nPulse Shaper connected!')
        print('Pulse Shaper name:',self.name)
        print('Pulse Shaper device:',self.device)
        print('Pulse Shaper time of creation:',self.time_of_creation)
        print('Pulse Shaper calibration file:',self.calibration_file)
        pass
    
    def set_calibration_file(self, calibration_file = None):
        self.calibration_file = calibration_file
        if self.calibration_file is not None:
            pulse_config = configparser.ConfigParser()
            pulse_config.read(self.calibration_file)
            
            self.poly_1 = float(pulse_config['pulse_shaper']['poly_1'])
            self.poly_0 = float(pulse_config['pulse_shaper']['poly_0'])
            self.slit_width = float(pulse_config['pulse_shaper']['slit_width'])
            self.slit_psf = float(pulse_config['pulse_shaper']['psf_width'])
            self.slit_motor_min = float(pulse_config['pulse_shaper']['slit_pos_min'])
            self.slit_motor_max = float(pulse_config['pulse_shaper']['slit_pos_max'])
            self.wavelength_min = self.motor_to_wavelength(self.slit_motor_min)
            self.wavelength_max = self.motor_to_wavelength(self.slit_motor_max)
            self.wavelength_range = np.abs(self.wavelength_max-self.wavelength_min)
        else:
            print('No calibration file provided')
    
    def excecute(self):
        # meant to be called from saved object and excecute device specific functions
        self.move_slit(excecute=True)

    def motor_to_wavelength(self,motor_position):
        return self.poly_1*motor_position+self.poly_0
    
    def wavelength_to_motor(self,wavelength):
        return (wavelength-self.poly_0)/self.poly_1

    def move_slit(self, slit_center_wavelength=None, pulse_object = None, excecute=False):
        if self.calibration_file is None: 
            print('No calibration file provided')
            return None
        
        if slit_center_wavelength is not None:
            self.slit_center_wavelength = slit_center_wavelength

        #self.motor_position = (self.slit_center_wavelength-self.poly_0)/self.poly_1
        self.motor_position = self.wavelength_to_motor(self.slit_center_wavelength)
        if self.motor_position < self.slit_motor_min or self.motor_position > self.slit_motor_max:
            pass
            #print('Caution: position out of confidence bounds.')
        
        if excecute:
            self.device.set_position(self.motor_position,excecute=True)
            sleep(self.motor_sleep)
        
        if pulse_object is not None:
            if type(pulse_object) is str:
                 pulse_object = pg.load_pulse(pulse_object)
            out_pulse_object = pulse_object.copy_pulse()
            out_pulse_object.clear_filter() # Do the error function approach -> erf(center-width/2-wavelengths) - erf(center+width/2-wavelengths))

            # hardcore convoluting -> computational heavy & numerical errors
            # out_pulse_object.add_filter_rectangle(unit = 'nm', central_f = self.slit_center_wavelength, width_f=self.slit_width, merging='+',transmission=1,cap_transmission=False)
            # print(np.sum(out_pulse_object.frequency_filter_x))
            # out_pulse_object.convolute_psf_filter(unit = 'nm', width_f = self.slit_psf, sig_or_fwhm= 'fwhm')

            # analytical solution with double error function 
            out_pulse_object.add_filter_double_erf(unit = 'nm', central_f = self.slit_center_wavelength, width_f=self.slit_width, rise_f = self.slit_psf,merging='+')
            
            self.mask = (np.sqrt(np.abs(out_pulse_object.frequency_filter_x)**2 + np.abs(out_pulse_object.frequency_filter_y)**2)**2)/2 
            
            out_pulse_object.apply_frequency_filter()
            
            self.mask_wavelengths = out_pulse_object.wavelengths
            out_pulse_object.clear_filter()
            
            return out_pulse_object
        else:
            wl_borders = [self.wavelength_min,self.wavelength_max]
            dummy_pulse = pg.PulseGenerator(t0=0,f0=min(wl_borders)-self.wavelength_range,fend = max(wl_borders)+self.wavelength_range,central_wavelength=np.mean(wl_borders))
            dummy_pulse.add_filter_double_erf(unit = 'nm', central_f = self.slit_center_wavelength, width_f=self.slit_width, rise_f = self.slit_psf,merging='+')
            
            self.mask = np.sqrt(np.abs(dummy_pulse.frequency_filter_x)**2 + np.abs(dummy_pulse.frequency_filter_y)**2)**2
            self.mask_wavelengths = dummy_pulse.wavelengths
            dummy_pulse.clear_filter()
            
        
    def calibrate(self, motor_pos_vec, wl_mat, counts_mat, initial_pulse, reference_pulse = None):
        if reference_pulse is not None:
            dummy_env_pulse = reference_pulse.copy_pulse() 
            dummy_slice_pulse = dummy_slice_pulse.copy_pulse() 
        else:
            dummy_env_pulse = pg.PulseGenerator(0,200,0.02,central_wavelength=np.mean(wl_mat))
            dummy_slice_pulse = dummy_env_pulse.copy_pulse()

        dummy_env_pulse.clear_all()
        dummy_slice_pulse.clear_all() 

        #dummy_env_pulse.add_spectrum_frequ(wl_mat[0],initial_pulse,plot=False)
        
        
        def _pulse_shaper_fitting(wavelength, wl_center, slit_width, psf, counts_scaling, background):
            # dummy_slice_pulse.clear_all()
            # dummy_slice_pulse.add_filter_rectangle(unit = 'nm', central_f = wl_center, width_f=slit_width, merging='+',transmission=1,cap_transmission=False)
            # dummy_slice_pulse.convolute_psf_filter(unit = 'nm', width_f = psf, sig_or_fwhm= 'fwhm')
            # dummy_slice_pulse.add_filter_rectangle(transmission=counts_scaling,merging='*',cap_transmission=False)
            # dummy_slice_pulse.add_filter_rectangle(transmission=background,merging='max',cap_transmission=False)
            # dummy_slice_pulse.add_spectrum_frequ(wl_mat[0],np.sqrt(initial_pulse/np.max(initial_pulse)),add_filter=True,plot=False,merging='*')
            # dummy_slice_pulse.add_rectangle_frequ(central_f = 0, width_f = np.inf, hight = 1)
            
            # dummy_slice_pulse.apply_frequency_filter()
            # dummy_slice_pulse_interp = _interpolate_pulse_to_spectrum(dummy_slice_pulse, wavelength,intensity=True)
            dummy_slice_pulse.clear_all()
            dummy_slice_pulse.add_spectrum_frequ(wl_mat[0],(initial_pulse-background)/np.max(initial_pulse),plot=False,power = None)
            # dummy_slice_pulse.add_filter_rectangle(unit = 'nm', central_f = wl_center, width_f=slit_width, merging='+',transmission=counts_scaling,cap_transmission=False)
            # dummy_slice_pulse.convolute_psf_filter(unit = 'nm', width_f = psf, sig_or_fwhm= 'fwhm') 
            ######
            dummy_slice_pulse.add_filter_double_erf(unit = 'nm', central_f = wl_center, width_f=slit_width, rise_f = psf,merging='+', transmission=counts_scaling,cap_transmission=False)
            dummy_slice_pulse.apply_frequency_filter()
            dummy_slice_pulse.add_rectangle_frequ(central_f = 0, width_f = np.inf, hight = np.sqrt(background))

            dummy_slice_pulse_interp = _interpolate_pulse_to_spectrum(dummy_slice_pulse, wavelength,intensity=True)
            return dummy_slice_pulse_interp
        
        wl_center_vec = []
        slit_width_vec = []
        psf_vec = []
        counts_scaling_vec = []
        background_vec = []
        for i_pos, pos in enumerate(motor_pos_vec):
            cur_slice = np.copy(counts_mat[i_pos])/np.max(initial_pulse)
            cur_wl = wl_mat[i_pos]
            
            L_initial_slit = np.where(cur_slice >= 0.2*(np.max(cur_slice)))
            initial_slit_width_guess = np.abs(cur_wl[L_initial_slit[0][0]] - cur_wl[L_initial_slit[0][-1]])
            if i_pos >= 0:
                guess_0 = [np.mean(cur_wl[L_initial_slit]), # wl_center 
                            initial_slit_width_guess, # slitwidth -> mayby like fwhm ?
                            0.1*initial_slit_width_guess, # slit slope
                            1,#np.max(cur_counts), #np.max(cur_counts)/np.max(np.abs(full_pulse_object.frequency_representation_x)**2),
                            np.min(cur_slice)]
            
            else: 
                guess_0 = [np.mean(cur_wl[L_initial_slit]), # wl_center 
                            np.mean(np.array(slit_width_vec)), # slitwidth -> mayby like fwhm ?
                            np.mean(np.array(psf_vec)), # slit slope
                            1,#np.mean(np.array(counts_scaling_vec)),#np.max(cur_counts), #np.max(cur_counts)/np.max(np.abs(full_pulse_object.frequency_representation_x)**2),
                            np.mean(np.array(background_vec))]

            bounds = ([np.min(cur_wl),0,0,0,0],
                        [np.max(cur_wl),np.inf,np.inf,np.inf,np.inf])
            popt, pcov = curve_fit(_pulse_shaper_fitting, cur_wl, cur_slice, p0=guess_0,bounds=bounds)
            wl_center_vec.append(popt[0])
            slit_width_vec.append(popt[1])
            psf_vec.append(popt[2])
            counts_scaling_vec.append(popt[3])
            background_vec.append(popt[4])

            print(str(i_pos+1) + '/' + str(len(motor_pos_vec)) + ' done')


            
            # plt.figure()
            # plt.plot(cur_wl,cur_slice,'b-')
            # plt.plot(cur_wl,_pulse_shaper_fitting(cur_wl,*popt),'r-')
            # plt.plot(cur_wl,_pulse_shaper_fitting(cur_wl,*guess_0),'g-')
            # plt.plot(cur_wl,initial_pulse/np.max(initial_pulse),'k-')
        
        
        poly_1, poly_0 = np.polyfit(motor_pos_vec,wl_center_vec,1)
        
        print('Calibration results:')
        print('Wavelength at motor 0: ' + str(poly_0))
        print('Wavelength per motor step: ' + str(poly_1))
        print('slit_width = ',np.mean(slit_width_vec),' +/- ',np.std(slit_width_vec))
        print('psf = ',np.mean(psf_vec),' +/- ',np.std(psf_vec))
        print('counts_scaling = ',np.mean(counts_scaling_vec),' +/- ',np.std(counts_scaling_vec))
        print('background = ',np.mean(background_vec),' +/- ',np.std(background_vec))
        
        self.calibration_file = 'Just_calibrated'
        self.poly_1 = poly_1
        self.poly_0 = poly_0
        self.slit_width = np.mean(slit_width_vec)
        self.slit_psf = np.mean(psf_vec)
        self.slit_motor_min = np.min(motor_pos_vec)
        self.slit_motor_max = np.max(motor_pos_vec)

        plt.figure()
        plt.title('Calibration results: translation')
        plt.plot(motor_pos_vec,wl_center_vec,'o')
        plt.plot(motor_pos_vec,np.polyval([poly_1,poly_0],motor_pos_vec),'r-')
        plt.xlabel('Motor position')
        plt.ylabel('Wavelength (nm)')
        plt.draw()

        
        plt.figure()
        plt.title('Calibration results: slit filter')
        for i_plot, psf_plot in enumerate(psf_vec):
            dummy_slice_pulse.clear_all()
            dummy_slice_pulse.add_filter_rectangle(unit = 'nm', central_f = wl_center_vec[0], width_f=slit_width_vec[i_plot])
            dummy_slice_pulse.convolute_psf_filter(unit = 'nm', width_f = psf_plot, sig_or_fwhm= 'fwhm')

            plt.plot(dummy_slice_pulse.wavelengths-wl_center_vec[0],np.abs(dummy_slice_pulse.frequency_filter_x),'k-',alpha=0.5)

        dummy_slice_pulse.clear_all()
        dummy_slice_pulse.add_filter_rectangle(unit = 'nm', central_f = wl_center_vec[0], width_f=self.slit_width)
        dummy_slice_pulse.convolute_psf_filter(unit = 'nm', width_f = self.slit_psf, sig_or_fwhm= 'fwhm')
        plt.plot(dummy_slice_pulse.wavelengths-wl_center_vec[0],np.abs(dummy_slice_pulse.frequency_filter_x),'r-')
        plt.xlim([-(self.slit_width+2*self.slit_psf)*1.5,(self.slit_width+2*self.slit_psf)*1.5])
        plt.xlabel('Wavelength (nm)')
        plt.ylabel('Transmission')
        plt.draw()

        plt.show()
        # indx_vec: vector of slit motor positions
        # wl_mat: matrix of wavelengths
        # counts_mat: matrix of counts
        # returns: poly_1, poly_0, slit_width, slit_psf, slit_motor_min, slit_motor_max
        # fit the slit motor positions to the wavelengths
        # fit the counts to the slit motor positions
        # fit the counts to the wavelengths
        # return the parameters
        pass

    def save_calibration(self, file_name = 'calibration_file'):
        pulse_config = configparser.ConfigParser()
        pulse_config['pulse_shaper'] = {             
                                        'poly_1':self.poly_1,
                                        'poly_0':self.poly_0,
                                        'slit_width':self.slit_width,
                                        'psf_width':self.slit_psf,
                                        'slit_pos_min':self.slit_motor_min,
                                        'slit_pos_max':self.slit_motor_max}
        with open(file_name+'.txt', 'w') as configfile:
            pulse_config.write(configfile)

        print('Calibration file saved as:'+str(file_name)+'.txt')
        pass
    
    def open_control(self,pulse_object = None,open_gui = True,parent_window = None,previous_control = None):
       
        control_object = psc.pulse_shaper_control(self,pulse_object = pulse_object, open_gui = open_gui,parent_window=parent_window,previous_control= previous_control)
        return control_object   
        
    def close(self):
        self.device.close()
    

class half_wave_plate():
    def __init__(self,name = 'HWP') -> None:
        self.name = name
        pass
    def rotate(self,pulse_object = None,device = None, angle=0, excecute=False):
        hwp_matrix = 1/2*np.array([[np.cos(2*angle),np.sin(2*angle)],
                               [np.sin(2*angle),-np.cos(2*angle)]])
        
        if excecute:
            print('To do')
            pass

        if pulse_object is not None:
            if type(pulse_object) is str:
                 pulse_object = pg.load_pulse(pulse_object)
            out_pulse_object = pulse_object.copy_pulse()
            out_pulse_object.clear_all()
            
            field_add_x = np.zeros_like(pulse_object.temporal_representation_x).astype(complex)
            field_add_y = np.zeros_like(pulse_object.temporal_representation_y).astype(complex)
            for i_temp, field_x in enumerate(pulse_object.temporal_representation_x):
                field_y = pulse_object.temporal_representation_y[i_temp]
                field_add_x[i_temp] = np.matmul(hwp_matrix,np.array([field_x,field_y]))[0]
                field_add_y[i_temp] = np.matmul(hwp_matrix,np.array([field_x,field_y]))[1]

            out_pulse_object._add_time(field_add_x,field_add_y)
            return out_pulse_object
        
class quarter_wave_plate():
    def __init__(self,name = 'QWP') -> None:
        self.name = name
        pass
    
class generic_wave_plate():
    def __init__(self,device = None, name = 'HWP', phase = np.pi, unit = 'deg') -> None:
        self.device = device
        self.name = name
        self.phase = phase
        self.unit = unit
    
    def set_unit(self,unit):
        self.unit = unit
    
    def get_unit(self):
        return self.unit
    
    def set_phase(self,phase):
        self.phase = phase
        
    def get_phase(self):
        return self.phase
    
    def generate_matrix(self,angle = 0):
        self.wp_mat = np.exp(-1j*self.phase/2)*np.array([[np.cos(angle)**2 + np.exp(1j*self.phase)*np.sin(angle)**2, (1-np.exp(1j*self.phase))*np.cos(angle)*np.sin(angle)],[(1-np.exp(1j*self.phase))*np.cos(angle)*np.sin(angle),np.sin(angle)**2 + np.exp(1j*self.phase)*np.cos(angle)**2]])
    
    def rotate(self,pulse_object = None, angle=0, excecute=False):
        if excecute:
            self.device.set_position(angle, excecute = excecute)
            #print('To do')
            pass
        
        if pulse_object is not None:
            if self.unit == 'deg':
                # range from 0 to 360 degrees
                angle = np.mod(angle,360)
                angle = np.deg2rad(angle)
            elif self.unit == 'rad':
                angle = np.mod(angle,2*np.pi)
            
            self.generate_matrix(angle)
            if type(pulse_object) is str:
                 pulse_object = pg.load_pulse(pulse_object)
            out_pulse_object = pulse_object.copy_pulse()
            out_pulse_object.clear_all()
            
            field_add_x = np.zeros_like(pulse_object.temporal_representation_x).astype(complex)
            field_add_y = np.zeros_like(pulse_object.temporal_representation_y).astype(complex)
            for i_temp, field_x in enumerate(pulse_object.temporal_representation_x):
                field_y = pulse_object.temporal_representation_y[i_temp]
                field_add_x[i_temp] = np.matmul(self.wp_mat,np.array([field_x,field_y]))[0]
                field_add_y[i_temp] = np.matmul(self.wp_mat,np.array([field_x,field_y]))[1]

            out_pulse_object._add_time(field_add_x,field_add_y)
            return out_pulse_object
    
    def open_control(self,pulse_object = None,open_gui = True,parent_window = None,previous_control = None):
       
        control_object = wpc.waveplate_control(self,pulse_object = pulse_object, open_gui = open_gui,parent_window=parent_window,previous_control= previous_control)
        return control_object
    
        
    def set_phase(self,phase):
        self.phase = phase
        
    def close(self):
        pass
        
class linear_polarizer():
    def __init__(self,name = 'LP', angle = 0, unit = 'deg') -> None:
        self.name = name
        self.angle = angle
        self.unit = unit
        pass
    
    def set_angle(self,angle):
        self.angle = angle

    def get_angle(self):
        return self.angle
    
    def rotate(self,pulse_object = None,device = None, angle=None, excecute=False, pol_vector = None):
        if angle is not None:
            self.set_angle(angle)
        if self.unit == 'deg':
            # range from 0 to 360 degrees
            angle = np.mod(self.angle,360)
            angle = np.deg2rad(self.angle)
        else:
            angle = np.mod(self.angle,2*np.pi)
        lp_matrix = np.array([[np.cos(angle)**2,np.cos(angle)*np.sin(angle)],
                               [np.cos(angle)*np.sin(angle),np.sin(angle)**2]])
        
        if excecute:
            print('To do')
            pass

        if pulse_object is not None:
            if type(pulse_object) is str:
                 pulse_object = pg.load_pulse(pulse_object)
            out_pulse_object = pulse_object.copy_pulse()
            out_pulse_object.clear_all()
            field_add_x = np.zeros_like(pulse_object.temporal_representation_x).astype(complex)
            field_add_y = np.zeros_like(pulse_object.temporal_representation_y).astype(complex)
            for i_temp, field_x in enumerate(pulse_object.temporal_representation_x):
                field_y = pulse_object.temporal_representation_y[i_temp]
                field_add_x[i_temp] = np.matmul(lp_matrix,np.array([field_x,field_y]))[0]
                field_add_y[i_temp] = np.matmul(lp_matrix,np.array([field_x,field_y]))[1]

            out_pulse_object._add_time(field_add_x,field_add_y)
            return out_pulse_object
        
        if pol_vector is not None:
            return np.matmul(lp_matrix,pol_vector)
  
class simulator():
    def __init__(self,qd_calibration = None,name = 'Simulator',temp_dir = '', sim_kind = 'ACE', decay = False, phonons = False, dot_size = 5,temperature = 4, hamiltonian = None, photon_pulse_bool = False, dipole_orientation = 0) -> None:
        if qd_calibration is None and hamiltonian is None:
            print('No calibration file provided')
            return None
        
        
        
        self.qd_calibration = qd_calibration
        self.name = name
        self.temp_dir = temp_dir
        self.decay = decay
        self.phonons = phonons
        self.dot_size = dot_size
        self.temperature = temperature
        self.sim_out = None
        #self.photon_emission = [0,0,0,0]
        #self.photon_wavelength = [0,0,0,0]
        self.photon_pulse = None
        self.photon_pulse_bool = photon_pulse_bool
        self.dipole_orientation = dipole_orientation
        if hamiltonian is not None:
            self.hamiltonian = self.read_Hamiltonian(hamiltonian)
        else:
            self.hamiltionan = None
        
        if sim_kind is None:
            # check operating system and default to 'ACE' for linux and 'Qutip' for windows
            if os.name == 'posix':
                sim_kind = 'ACE'
            elif os.name == 'nt':
                sim_kind = 'Qutip'
        self.sim_kind = sim_kind
        self.hwp = generic_wave_plate(name = 'HWP', phase = np.pi, unit = 'deg')
        self.print_info()
        
        self.decay_scale_x = 0
        self.decay_scale_y = 0
        
        self.set_ace_six_level()
        self.refresh_num_states()
        
    
    def refresh_num_states(self):
        if any(self.sim_kind.lower() == x for x in ['ace','ace_ds','qutip']):
            self.num_states = 4 
            
        if any(self.sim_kind.lower() == x for x in ['ace_6ls']):
            self.num_states = 6 
            
            
        
        self.photon_emission = []
        self.photon_wavelength = []
        self.photon_polarisation = []
        for i in range(self.num_states + self.num_states - 4): # stupid fix, but we need 8 transitions for 6 level system
            self.photon_emission.append(0)
            self.photon_wavelength.append(0)
            self.photon_polarisation.append([])
    
    def set_qd_calibration(self,qd_calibration):
        self.qd_calibration = qd_calibration
        self.print_info()
    
    def set_phonons(self,phonons):
        self.phonons = phonons
    
    def set_temperature(self,temperature):
        self.temperature = temperature    
    
    def toggle_phonons(self):
        self.phonons = not self.phonons
    
    def set_dipole_orientation(self,dipole_orientation):
        self.dipole_orientation = dipole_orientation
    
    def get_num_states(self):
        return self.num_states
    
    def set_decay(self,decay = True):
        self.decay = decay

    def print_info(self):
        print('\nSimulator connected!')
        print('Quantum dot calibration file:',self.qd_calibration)
        print('Simulation method:',self.sim_kind)
        print('Temp directory:',self.temp_dir)

    def read_Hamiltonian(self, hamiltonian_file):
        config = configparser.ConfigParser()
        config.read(hamiltonian_file) 

        dimension = int(config['GENERAL']['dimension'])
        rotatin_frame = float(config['GENERAL']['rotating_frame'])
        field_coupling = int(config['GENERAL']['field_coupling'])
        rows = []
        for i in range(dimension): 
            index = i+1 
            rows.append(config['HAMILTONIAN']['row_'+str(index)])

        
        if self.sim_kind.lower() == 'qutip':
            current_row = rows[0].split(',')
            H_energy = float(current_row[0])*qt.basis(dimension,0)*qt.basis(dimension,0).dag()
            for i in range(dimension-1):
                current_row = rows[i+1].split(',')
                H_energy += (float(current_row[i+1])-rotatin_frame)/HBAR*1e-3*qt.basis(dimension,i+1)*qt.basis(dimension,i+1).dag() 
            
            print(H_energy)
            coupling_index = []
            coupling_type = []
            for i in range(dimension):
                current_row = rows[i].split(',')
                for j in range(i+1,dimension): 
                    current_cell = str(current_row[j]).strip()
                    if current_cell[0].lower() == 'c':
                        coupling_index.append([i,j])
                        coupling_type.append(int(current_cell[1:])) 
            
            def H_qutip():
                pass
                        
                    
        return None

    def simulate(self,pulse_object,sim_dt = None, dipole_moment = 1, plot = False):
        pulse_object.set_rotating_frame(self.qd_calibration)
        self.refresh_num_states()
        if self.dipole_orientation != 0:
            pulse_object = self.hwp.rotate(pulse_object,angle=self.dipole_orientation/2,excecute=False)
        if self.sim_kind.lower() == 'ace':
            self.decay_scale_x = 0
            self.decay_scale_y = 0
            self.sim_out = self.ace_four_level(pulse_object,sim_dt,self.decay,self.phonons,plot,dipole_moment)
        
        elif self.sim_kind.lower() == 'ace_ds':
            self.decay_scale_x = 0
            self.decay_scale_y = 0
            self.sim_out = self.ace_four_level_ds(pulse_object,sim_dt,self.decay,self.phonons,plot,dipole_moment)

        elif self.sim_kind.lower() == 'qutip':
            self.decay_scale_x = 0
            self.decay_scale_y = 0
            self.sim_out =  self.qutip_four_level(pulse_object,sim_dt,self.decay,self.phonons,plot,dipole_moment*np.pi)
            
        elif self.sim_kind == 'ace_6ls':
            self.sim_out = self.ace_six_level(pulse_object,sim_dt,self.decay,plot,dipole_moment)
        
        else:
            print('No valid simulation kind provided')
        self.photon_wavelength[0] = pulse_object._Units_inverse(pulse_object.exciton_x_emission, 'nm')
        self.photon_wavelength[1] = pulse_object._Units_inverse(pulse_object.exciton_y_emission, 'nm')
        self.photon_wavelength[2] = pulse_object._Units_inverse(pulse_object.biexciton_x_emission, 'nm')
        self.photon_wavelength[3] = pulse_object._Units_inverse(pulse_object.biexciton_y_emission, 'nm')
        
        self.photon_polarisation[0] = [1,0]
        self.photon_polarisation[1] = [0,1]
        self.photon_polarisation[2] = [1,0]
        self.photon_polarisation[3] = [0,1]
            
        
        if self.decay:
            # photon emission estimated from integral over time
            self.photon_emission[0] = np.trapz(x = np.real(self.sim_out[0]),y = np.abs(self.sim_out[2]))
            # normalize to exciton lifetime
            self.photon_emission[0] *= 1/pulse_object.lifetime_exciton*(1-self.decay_scale_x)
            # plus last value of dynamics + half of last biexcion value
            self.photon_emission[0] += np.abs(self.sim_out[2][-1]) + np.abs(self.sim_out[4][-1])/2*(1-self.decay_scale_x)

            self.photon_emission[1] = np.trapz(x = np.real(self.sim_out[0]),y = np.abs(self.sim_out[3]))
            self.photon_emission[1] *= 1/pulse_object.lifetime_exciton*(1-self.decay_scale_y)
            self.photon_emission[1] += np.abs(self.sim_out[3][-1]) + np.abs(self.sim_out[4][-1])/2*(1-self.decay_scale_y)

            # for biexciton emission same, but split in two 
            self.photon_emission[2] = np.trapz(x = np.real(self.sim_out[0]),y = np.abs(self.sim_out[4]))/2
            self.photon_emission[2] *= 1/pulse_object.lifetime_biexciton*(1-self.decay_scale_x)
            self.photon_emission[2] += np.abs(self.sim_out[4][-1])/2*(1-self.decay_scale_x)

            self.photon_emission[3] = np.trapz(x = np.real(self.sim_out[0]),y = np.abs(self.sim_out[4]))/2
            self.photon_emission[3] *= 1/pulse_object.lifetime_biexciton*(1-self.decay_scale_y)
            self.photon_emission[3] += np.abs(self.sim_out[4][-1])/2*(1-self.decay_scale_y)
        else:
            self.photon_emission[0] = np.abs(self.sim_out[2][-1]) + np.abs(self.sim_out[4][-1])/2*(1-self.decay_scale_x)
            self.photon_emission[1] = np.abs(self.sim_out[3][-1]) + np.abs(self.sim_out[4][-1])/2*(1-self.decay_scale_y)
            self.photon_emission[2] = np.abs(self.sim_out[4][-1])/2*(1-self.decay_scale_x)
            self.photon_emission[3] = np.abs(self.sim_out[4][-1])/2*(1-self.decay_scale_y)
            
        if self.num_states == 6:
            #print(self.six_ls_energy)
            new_photon_energy = np.zeros(8)
            new_photon_energy[0] = self.six_ls_energy[1] - self.six_ls_energy[0]
            new_photon_energy[1] = self.six_ls_energy[2] - self.six_ls_energy[0]
            new_photon_energy[2] = self.six_ls_energy[5] - self.six_ls_energy[1]
            new_photon_energy[3] = self.six_ls_energy[5] - self.six_ls_energy[2]
            new_photon_energy[4] = self.six_ls_energy[3] - self.six_ls_energy[0]
            new_photon_energy[5] = self.six_ls_energy[4] - self.six_ls_energy[0]
            new_photon_energy[6] = self.six_ls_energy[5] - self.six_ls_energy[3]
            new_photon_energy[7] = self.six_ls_energy[5] - self.six_ls_energy[4]
            
            for i in range(8):
                #print(pulse_object._Units_inverse(pulse_object._Units(new_photon_energy[i],'meV'),'nm'))
                self.photon_wavelength[i] = pulse_object._Units_inverse(pulse_object._Units(new_photon_energy[i],'meV'),'nm')
            #for i in range(6):
                #print(pulse_object._Units_inverse(pulse_object._Units(self.six_ls_energy[i],'meV'),'nm'))
            #print(self.photon_wavelength)
            # rethink the emission wavelength.. but looks good
            
            # config = configparser.ConfigParser()
            # config.read(self.qd_calibration)
            # dark_wavelength = float(config['EMISSION']['dark_wavelength'])
            # dark_splitting = pulse_object._Units_inverse(float(config['SPLITTING']['fss_dark'])*1e-3,'meV')
            # dark_splitting = pulse_object._Units(dark_splitting,'nm')
                                                  
            # self.photon_wavelength[4] =  dark_wavelength + dark_splitting
            # self.photon_wavelength[5] =  dark_wavelength - dark_splitting
            # self.photon_wavelength[6] = self.photon_wavelength[2] + self.photon_wavelength[0]-self.photon_wavelength[4]
            # self.photon_wavelength[7] = self.photon_wavelength[3] + self.photon_wavelength[1]-self.photon_wavelength[5]
           
            
            if self.decay:
                self.photon_emission[4] = np.trapz(x = np.real(self.sim_out[0]),y = np.abs(self.sim_out[5]))
                self.photon_emission[4] *= 1/pulse_object.lifetime_exciton*(self.decay_scale_x)
                self.photon_emission[4] += np.abs(self.sim_out[5][-1]) + np.abs(self.sim_out[4][-1])/2*self.decay_scale_x
                
                self.photon_emission[5] = np.trapz(x = np.real(self.sim_out[0]),y = np.abs(self.sim_out[6]))
                self.photon_emission[5] *= 1/pulse_object.lifetime_exciton*(self.decay_scale_y)
                self.photon_emission[5] += np.abs(self.sim_out[6][-1]) + np.abs(self.sim_out[4][-1])/2*self.decay_scale_y
                
                self.photon_emission[6] = np.trapz(x = np.real(self.sim_out[0]),y = np.abs(self.sim_out[4]))/2
                self.photon_emission[6] *= 1/pulse_object.lifetime_biexciton*(self.decay_scale_x)
                self.photon_emission[6] += np.abs(self.sim_out[4][-1])/2*self.decay_scale_x
                
                self.photon_emission[7] = np.trapz(x = np.real(self.sim_out[0]),y = np.abs(self.sim_out[4]))/2
                self.photon_emission[7] *= 1/pulse_object.lifetime_biexciton*(self.decay_scale_y)
                self.photon_emission[7] += np.abs(self.sim_out[4][-1])/2*self.decay_scale_y
            else:
                self.photon_emission[4] = np.abs(self.sim_out[5][-1])  + np.abs(self.sim_out[4][-1])/2*self.decay_scale_x
                self.photon_emission[5] = np.abs(self.sim_out[6][-1]) + np.abs(self.sim_out[4][-1])/2*self.decay_scale_y
                self.photon_emission[6] = np.abs(self.sim_out[4][-1])/2*self.decay_scale_x
                self.photon_emission[7] = np.abs(self.sim_out[4][-1])/2*self.decay_scale_y
            
            self.photon_polarisation[4] = [1,0]
            self.photon_polarisation[5] = [0,1]
            self.photon_polarisation[6] = [1,0]
            self.photon_polarisation[7] = [0,1]

        if self.photon_pulse_bool:
            pass
        self.photon_pulse = pulse_object.copy_pulse()
        self.photon_pulse.clear_all()
        self.photon_pulse.add_time_field(self.sim_out[0],np.sqrt(self.sim_out[2]), polarisation=[1,0],frequency=pulse_object.exciton_x_emission,power = self.photon_emission[0])
        self.photon_pulse.add_time_field(self.sim_out[0],np.sqrt(self.sim_out[3]), polarisation=[0,1],frequency=pulse_object.exciton_y_emission, power = self.photon_emission[1])
        self.photon_pulse.add_time_field(self.sim_out[0],np.sqrt(self.sim_out[4]), polarisation=[1,0],frequency=pulse_object.biexciton_x_emission, power = self.photon_emission[2])
        self.photon_pulse.add_time_field(self.sim_out[0],np.sqrt(self.sim_out[4]), polarisation=[0,1],frequency=pulse_object.biexciton_y_emission, power = self.photon_emission[3])

        
        
        
        # biexciton_x_emission = np.abs(sim_input[4][-1])/2 
        #     biexciton_y_emission = np.abs(sim_input[4][-1])/2 
        #     exciton_x_emission = np.abs(sim_input[2][-1]) + np.abs(sim_input[4][-1])/2
        #     exciton_y_emission = np.abs(sim_input[3][-1]) + np.abs(sim_input[4][-1])/2

        # if unit.lower()[0] == 'n':
        #     plot_domain = pulse_object.wavelengths
        #     position_x_h = pulse_object._Units_inverse(pulse_object.exciton_x_emission, 'nm')
        #     position_x_v = pulse_object._Units_inverse(pulse_object.exciton_y_emission, 'nm')
        #     postion_xx_h = pulse_object._Units_inverse(pulse_object.biexciton_x_emission, 'nm')
        #     postion_xx_v = pulse_object._Units_inverse(pulse_object.biexciton_y_emission, 'nm')
        
    def get_simulation_results(self):
        return self.sim_out

    def ace_four_level(self,pulse_object,sim_dt = None, decay = False, phonons = False,plot = False,dipole_moment = 1):
        if type(pulse_object) is str:
            pulse_object = pg.load_pulse(pulse_object)
        sim_pulse_object = pulse_object.copy_pulse()
        sim_pulse_object.clear_filter()
        sim_pulse_object.add_filter_rectangle(transmission=dipole_moment,cap_transmission=False)
        sim_pulse_object.apply_frequency_filter()
        
        
        if sim_dt is None:
            sim_dt = sim_pulse_object.dt

        # find next power of 2 dt*2^n thats just greater than 10 ps 
        n = np.ceil(np.log2(10/sim_dt))     
        t_mem = sim_dt*2**n
        pulse_x, pulse_y = sim_pulse_object.generate_pulsefiles(temp_dir=self.temp_dir,precision=8)
        t,g,x,y,b = biexciton(sim_pulse_object.t0,sim_pulse_object.tend,dt=sim_dt,delta_xy=0, delta_b=4, temp_dir=self.temp_dir,
                                lindblad=decay,pulse_file_x=pulse_x,pulse_file_y=pulse_y,
                                output_ops=['|0><0|_4','|1><1|_4','|2><2|_4','|3><3|_4'],phonons=self.phonons, ae=self.dot_size,
                                temperature=self.temperature,suffix='pulse',calibration_file=self.qd_calibration,t_mem = t_mem)
        
        if plot:
            sim_pulse_object.plot_pulses(domain='nm',plot_frequ_intensity=True,sim_input = [t,g,x,y,b],sim_label=['g','x','y','b'])
        
        return [np.real(t),np.abs(g),np.abs(x),np.abs(y),np.abs(b)]
    
    def set_ace_six_level(self, b_x = None, b_z = None, b_field_frame = True):
        if b_x is not None:
            self.b_x = b_x
        else: 
            self.b_x = 0
        if b_z is not None:
            self.b_z = b_z
        else: 
            self.b_z = 0
        self.b_field_frame = b_field_frame
    
    def set_mag_field(self,b_x,b_z):
        self.b_x = b_x
        self.b_z = b_z
    
    def set_b_field_frame(self,b_field_frame):
        self.b_field_frame = b_field_frame
        
    def toggle_b_field_frame(self):
        self.b_field_frame = not self.b_field_frame
    
    def ace_six_level(self,pulse_object,sim_dt = None, decay = False, plot = False,dipole_moment = 1):
        
        if type(pulse_object) is str:
            pulse_object = pg.load_pulse(pulse_object)
        sim_pulse_object = pulse_object.copy_pulse()
        sim_pulse_object.clear_filter()
        sim_pulse_object.add_filter_rectangle(transmission = dipole_moment,cap_transmission=False)
        sim_pulse_object.apply_frequency_filter()
        
        if sim_dt is None:
            sim_dt = sim_pulse_object.dt
            
        n = np.ceil(np.log2(10/sim_dt))
        t_mem = sim_dt*2**n
        pulse_x, pulse_y = sim_pulse_object.generate_pulsefiles(temp_dir=self.temp_dir,precision=8) 
        ds_t, _, ds_occ, _, rho = sixls_linear_dressed_states(sim_pulse_object.t0,sim_pulse_object.tend,dt=sim_dt,pulse_file_x=pulse_x,pulse_file_y=pulse_y,temp_dir=self.temp_dir, suffix='pulse', initial = '|0><0|_6', rf = False, calibration_file = self.qd_calibration, bx = self.b_x, bz = self.b_z, lindblad = decay,phonons = self.phonons, ae=self.dot_size, temperature=self.temperature) 
        
        return_vector = [np.real(ds_t)]
        if self.b_field_frame:
            rho, index = self.bx_field_basis_transformation(rho,self.b_x,self.qd_calibration, bz=self.b_z)
            new_index = [index[0],index[1],index[2],index[5],index[3],index[4]]
            for i in new_index: #range(self.num_states)
                cur_return = []
                for j in range(len(ds_t)):
                    cur_return.append(np.abs(rho[j][i,i])) 
                return_vector.append(cur_return)
            
        else:
            _,_ = self.bx_field_basis_transformation(rho,self.b_x,self.qd_calibration, bz=self.b_z)
            for i in [0,1,2,5,3,4]: #range(self.num_states)
                cur_return = []
                for j in range(len(ds_t)):
                    cur_return.append(np.abs(rho[j][i,i])) 
                return_vector.append(cur_return)
        
        return return_vector
        
        #ds_t, ds_e, ds_occ, ds_color, rho = sixls_linear_dressed_states(t0,t_end,dt=dt,pulse_file_x=pulse_x,pulse_file_y=pulse_y,temp_dir='sim_dump/', suffix='_ds2',
        #                       verbose=False,initial='|5><5|_6',
                            # lindblad=False,rf = True, rf_file = ph_x,bx=bx,temperature = 1.5, ae = 5, phonons = False,
                            #  calibration_file=calib_file, make_transparent=[0,2,4])
    
    def light_dressed_states(self,pulse_object):
        if type(pulse_object) is str:
            pulse_object = pg.load_pulse(pulse_object)
        sim_pulse_object = pulse_object.copy_pulse()
        E_X, E_Y, E_S, E_F, E_B, _, _, g_ex, g_hx, g_ez, g_hz = read_calibration_file(self.qd_calibration)
        mu_b = 5.7882818012e-2   # meV/T
        hbar = 0.6582173  # meV*ps
        #field_x, field_y = sim_pulse_object.generate_field_functions() # _lab_frame
        
        
        if self.num_states == 4:
            energy_mat = [[],[],[],[]]
            def H_func(field_x,field_y):
                H = np.array([[0,field_x,field_y,0],
                            [np.conj(field_x),E_X,0,field_x],
                            [np.conj(field_y),0,E_Y,field_y],
                            [0,np.conj(field_x),np.conj(field_y),E_B]],dtype=complex)
                
                
                return H
        elif self.num_states == 6:
            energy_mat = [[],[],[],[],[],[]]
            A = -0.5*mu_b*self.b_x*(g_ex+g_hx)
            B = -0.5*mu_b*self.b_x*(g_ex-g_hx) 
            
            C = -1j*0.5*mu_b*self.b_z*(g_ez-3*g_hz)
            D = 1j*0.5*mu_b*self.b_z*(g_ez+3*g_hz)
            def H_func(field_x,field_y):
                H = np.array([[0,field_x,field_y,0,0,0],
                    [np.conj(field_x),E_X,C,A,0,field_x],
                    [np.conj(field_y),-C,E_Y,0,B,field_y],
                    [0,A,0,E_S,D,0],
                    [0,0,B,-D,E_F,0],
                    [0,np.conj(field_x),np.conj(field_y),0,0,E_B]],dtype=complex)
                return H
        
        #self.system_hamiltonian = H_func(0,0)
        
        for i, t in enumerate(sim_pulse_object.time):
            eigenvalue, eigenvector = np.linalg.eig(H_func(sim_pulse_object.temporal_representation_x[i],sim_pulse_object.temporal_representation_y[i]))
            index = []
            for vec in eigenvector:
                index.append(np.argmax(np.abs(vec)))
            for i in range(self.num_states):
                energy_mat[i].append(np.real(eigenvalue[index[i]]))
                
        return [sim_pulse_object.time,energy_mat]
    
    # def save_system_hamiltonian_numpy(self,filename):
    #     np.save(filename,self.system_hamiltonian)
    
    # def save_system_hamiltonian_txt(self,filename):
    #     if self.sim_kind == 'ace':
    #         system_ham_save_str = np.array([])
            
    #         for i in range(self.num_states):
    #             for j in range(self.num_states):
    #                  #"{}*(|1><3|_6 + |3><1|_6 )".format(-0.5*mu_b*bx*(g_ex+g_hx))
    #                 system_ham_save_str.append(str(self.system_hamiltonian[i,j])+' ')
            
    #     elif self.sim_kind == 'qutip':
    #         pass
    #     else: 
    #         self.save_system_hamiltonian_numpy(filename)
            
        
        
        
        
    
    def bx_field_basis_transformation(self,rho,bx,calibration_file, bz = 0): 
        E_X, E_Y, E_S, E_F, E_B, _, _, g_ex, g_hx, g_ez, g_hz = read_calibration_file(calibration_file)
        mu_b = 5.7882818012e-2   # meV/T
        hbar = 0.6582173  # meV*ps
        A = -0.5*mu_b*bx*(g_ex+g_hx)
        B = -0.5*mu_b*bx*(g_ex-g_hx) 
        
        #print('+ mixing'+str(A))
        #print('- mixing'+str(B))
        
        C = -1j*0.5*mu_b*bz*(g_ez-3*g_hz)
        D = 1j*0.5*mu_b*bz*(g_ez+3*g_hz) # -
        # system_op.append("i*{}*(|1><2|_6 -|2><1|_6)".format(0.5*mu_b*bz*(g_ez-3*g_hz)))
        # system_op.append("i*{}*(|4><3|_6 - |3><4|_6 )".format(-0.5*mu_b*bz*(g_ez+3*g_hz)))
        bare_state_index = np.argsort([0,E_X,E_Y,E_S,E_F,E_B])

        H = np.array([[0,0,0,0,0,0],
                    [0,E_X,C,A,0,0],
                    [0,-C,E_Y,0,B,0],
                    [0,A,0,E_S,D,0],
                    [0,0,B,-D,E_F,0],
                    [0,0,0,0,0,E_B]])

        sub_H_x_p = np.array([[E_X,A],
                              [A,E_S]])
        sub_H_x_m = np.array([[E_Y,B],
                              [B,E_F]])
        
        _, U_x_p = np.linalg.eig(sub_H_x_p)
        _, U_x_m = np.linalg.eig(sub_H_x_m) 
        
        
        self.decay_scale_x = min(abs(U_x_p[0]))**2
        self.decay_scale_y = min(abs(U_x_m[0]))**2
        
        eigenvalue, eigenvector = np.linalg.eig(H)
        
        index = []
        for vec in eigenvector:
            index.append(np.argmax(np.abs(vec)))
        
        self.six_ls_energy = np.real(eigenvalue[index])
        #print(np.real(eigenvalue[index]))
        # energy shift -> Monday!!!
        
        #print(eigenvector.transpose())
        #dress_state_index = np.argsort(eigenvalue)
        #print([0,E_X,E_Y,E_S,E_F,E_B])
        #print(bare_state_index)
        #print(dress_state_index)
        new_rho = []
        s0 = []
        s1 = []
        s2 = []
        s3 = []
        s4 = []
        s5 = []
        for r in rho: 
            new_rho.append(eigenvector.transpose()@r@eigenvector)

            s0.append(new_rho[-1][0,0])
            s1.append(new_rho[-1][1,1])
            s2.append(new_rho[-1][2,2])
            s3.append(new_rho[-1][3,3])
            s4.append(new_rho[-1][4,4])
            s5.append(new_rho[-1][5,5])
        # print(H)
        # print('bare_states')
        # print(bare_state_index)
        # print('dressed_states')
        # print(dress_state_index)

        # print(eigenvalue)
        return np.array(new_rho), index 
    
    def ace_four_level_ds(self,pulse_object,sim_dt = None, decay = False, phonons = False,plot = False,dipole_moment = 1): 
        if type(pulse_object) is str:
            pulse_object = pg.load_pulse(pulse_object)
        sim_pulse_object = pulse_object.copy_pulse()
        sim_pulse_object.clear_filter()
        sim_pulse_object.add_filter_rectangle(transmission=dipole_moment,cap_transmission=False)
        sim_pulse_object.apply_frequency_filter()
        
        if sim_dt is None:
            sim_dt = sim_pulse_object.dt

        # find next power of 2 dt*2^n thats just greater than 10 ps 
        n = np.ceil(np.log2(10/sim_dt))     
        t_mem = sim_dt*2**n
        pulse_x, pulse_y = sim_pulse_object.generate_pulsefiles(temp_dir=self.temp_dir,precision=8)
        t,ev,ds_occ,s_colors,rho = biexciton_dressed_states(sim_pulse_object.t0,sim_pulse_object.tend,dt=sim_dt,delta_xy=0, delta_b=4, temp_dir=self.temp_dir,
                                lindblad=decay,pulse_file_x=pulse_x,pulse_file_y=pulse_y,
                                output_ops=['|0><0|_4','|1><1|_4','|2><2|_4','|3><3|_4'],phonons=self.phonons, ae=self.dot_size,
                                temperature=self.temperature,suffix='pulse',calibration_file=self.qd_calibration,t_mem = t_mem)
        t = np.real(t)
        g = rho[:,0,0]
        x = rho[:,1,1]
        y = rho[:,2,2]
        b = rho[:,3,3]
    

        s0,s1,s2,s3,e0,e1,e2,e3,new_rho = self.four_level_pulse_basis_transformation(rho,sim_pulse_object)

        # e0 = ev[:,0]
        # e1 = ev[:,1]
        # e2 = ev[:,2]
        # e3 = ev[:,3]

        hbar = 0.6582173  # meV*ps

        phot_pulse = pg.PulseGenerator(t[0],t[-1],np.abs(t[1]-t[0]),calibration_file=self.qd_calibration)
        phot_x = np.array(x,dtype=complex)*np.exp(-1j*np.ones_like(x)*phot_pulse.exciton_x_emission/hbar*2*np.pi*t)
        phot_pulse._add_time(phot_x,np.zeros_like(phot_x))

        self.phot_pulse = phot_pulse
        #phot_pulse._add_time(np.array(s0,dtype=complex)*np.exp(-1j*np.array(e0,dtype=complex)/hbar*2*np.pi*phot_pulse.time),np.zeros_like(s0))
        if plot:
            sim_pulse_object.plot_pulses(domain='nm',plot_frequ_intensity=True,sim_input = [t,g,x,y,b,s0,s1,s2,s3],sim_label=['g','x','y','b','ds0','ds1','ds2','ds3'])
            # plt.figure()
            # plt.plot(t,e0,label='s0')
            # plt.plot(t,e1,label='s1')
            # plt.plot(t,e2,label='s2')
            # plt.plot(t,e3,label='s3')
            # plt.xlabel('time / ps')
            # plt.ylabel('DS energy / meV')
            # plt.title('Dressed state energy')
            # #plt.ylim([-1,1])
            
            # plt.figure()
            # plt.plot(t,np.abs(g),label= 'g')
            # plt.plot(t,np.abs(x),label= 'x')
            # plt.plot(t,np.real(rho[:,0,1]),label='gx')
            # plt.plot(t,np.real(rho[:,1,0]), label = 'xg')

            # phot_pulse.plot_pulses(domain='nm')
            # plt.show()

        return [np.real(t),np.abs(s0),np.abs(s1),np.abs(s2),np.abs(s3)]
        

    def qutip_four_level(self,pulse_object,sim_dt = None, decay = False, phonons = False,plot = False,dipole_moment = np.pi):
        if type(pulse_object) is str:
            pulse_object = pg.load_pulse(pulse_object)

        sim_pulse_object = pulse_object.copy_pulse()
        sim_pulse_object.clear_filter()
        sim_pulse_object.add_filter_rectangle(transmission=dipole_moment,cap_transmission=False)
        sim_pulse_object.apply_frequency_filter()

        if sim_dt is None:
            sim_dt = sim_pulse_object.dt

        pulse_x, pulse_y = sim_pulse_object.generate_field_functions() 

        t,g,x,y,b, _, _, _, _, _, _ = fourlevel_system(t0=sim_pulse_object.t0,tend=sim_pulse_object.tend,dt=sim_dt,calibration_file=self.qd_calibration,pulse_x=pulse_x, pulse_y=pulse_y,timeAxis='PULSE', collapse=decay,timeAxis_smart=False) 

        if plot:
            sim_pulse_object.plot_pulses(domain='nm',plot_frequ_intensity=True,sim_input = [t,g,x,y,b],sim_label=['g','x','y','b'])
        
        return [np.real(t),np.abs(g),np.abs(x),np.abs(y),np.abs(b)]
        # t_axis, g_occ, x_occ, y_occ, b_occ, polar_gx, polar_xb, polar_gb, t_axis, pulse1, pulse2

        #tbd
    

    def four_level_pulse_basis_transformation(self,rho,pulse_object): 
        E_X, E_Y, _, _, E_B, _, _, _, _, _, _ = read_calibration_file(self.qd_calibration)
        mu_b = 5.7882818012e-2   # meV/T
        hbar = 0.6582173  # meV*ps


        bare_state_index = np.argsort([0,E_X,E_Y,E_B])

        # x_upper = np.conj(pulse_object.temporal_representation_x)
        # y_upper = np.conj(pulse_object.temporal_representation_x)
        # x_lower = pulse_object.temporal_representation_x
        # y_lower = pulse_object.temporal_representation_y

        H = np.array([[0,0,0,0],
                    [0,E_X,0,0],
                    [0,0,E_Y,0],
                    [0,0,0,E_B]],dtype=complex)
        #print(eigenvector)
        #print(eigenvector.transpose())
        #dress_state_index = np.argsort(eigenvalue)
        #print([0,E_X,E_Y,E_S,E_F,E_B])
        #print(bare_state_index)
        #print(dress_state_index)
        new_rho = []
        s0 = []
        s1 = []
        s2 = []
        s3 = []

        e0 = []
        e1 = []
        e2 = []
        e3 = []
        for i, r in enumerate(rho): 
            H[0,1] = np.conj(pulse_object.temporal_representation_x[i])
            H[0,2] = np.conj(pulse_object.temporal_representation_y[i])
            H[1,3] = np.conj(pulse_object.temporal_representation_x[i])
            H[2,3] = np.conj(pulse_object.temporal_representation_y[i])

            H[1,0] = pulse_object.temporal_representation_x[i]
            H[2,0] = pulse_object.temporal_representation_y[i]
            H[3,1] = pulse_object.temporal_representation_x[i]
            H[3,2] = pulse_object.temporal_representation_y[i]

            
            eigenvalue, eigenvector = np.linalg.eig(H)
            new_rho.append(eigenvector.transpose()@r@eigenvector)

            s0.append(new_rho[-1][0,0])
            s1.append(new_rho[-1][1,1])
            s2.append(new_rho[-1][2,2])
            s3.append(new_rho[-1][3,3])

            e0.append(np.real(eigenvalue[0]))
            e1.append(np.real(eigenvalue[1]))
            e2.append(np.real(eigenvalue[2]))
            e3.append(np.real(eigenvalue[3]))
        # print(H)
        # print('bare_states')
        # print(bare_state_index)
        # print('dressed_states')
        # print(dress_state_index)

        # print(eigenvalue)
        return s0,s1,s2,s3,e0,e1,e2,e3,np.array(new_rho)


    def simulate_folder(self,folder,decay = False, phonons = False,plot = False,dipole_moment = 1, qutip = False):
        files = os.listdir(folder)
        files = [os.path.join(folder, f) for f in files]
        files.sort(key=lambda x: os.path.getmtime(x))
        c = 0
        sim_output = []
        for pulse_file in files:
            print('Simulation progress: ',str(c+1),'/',str(len(files)))
            cur_pulse = pg.load_pulse(pulse_file)
            if qutip:
                sim_output.append(self.qutip_four_level(cur_pulse,decay,phonons,plot,dipole_moment*np.pi))
            else:
                sim_output.append(self.ace_four_level(cur_pulse,decay,phonons,plot,dipole_moment))
            c += 1
        return sim_output
        
    def spectral_plot_four_level(self,pulse_object,sim_input, model = 'static', unit = 'nm',pulse_intesity = True, plot = True):
        if type(pulse_object) is str:
            pulse_object = pg.load_pulse(pulse_object)
        if model == 'static':
            biexciton_x_emission = np.abs(sim_input[4][-1])/2 
            biexciton_y_emission = np.abs(sim_input[4][-1])/2 
            exciton_x_emission = np.abs(sim_input[2][-1]) + np.abs(sim_input[4][-1])/2
            exciton_y_emission = np.abs(sim_input[3][-1]) + np.abs(sim_input[4][-1])/2

        if unit.lower()[0] == 'n':
            plot_domain = pulse_object.wavelengths
            position_x_h = pulse_object._Units_inverse(pulse_object.exciton_x_emission, 'nm')
            position_x_v = pulse_object._Units_inverse(pulse_object.exciton_y_emission, 'nm')
            postion_xx_h = pulse_object._Units_inverse(pulse_object.biexciton_x_emission, 'nm')
            postion_xx_v = pulse_object._Units_inverse(pulse_object.biexciton_y_emission, 'nm')
        elif unit.lower()[0] == 'f':
            plot_domain = pulse_object.frequencies + pulse_object.central_frequency
            position_x_h = pulse_object.exciton_x_emission
            position_x_v = pulse_object.exciton_y_emission
            postion_xx_h = pulse_object.biexciton_x_emission
            postion_xx_v = pulse_object.biexciton_y_emission
        elif unit.lower()[0] == 'm':
            plot_domain = pulse_object.energies + pulse_object.central_energy
            position_x_h = pulse_object._Units_inverse(pulse_object.exciton_x_emission, 'meV')
            position_x_v = pulse_object._Units_inverse(pulse_object.exciton_y_emission, 'meV')
            postion_xx_h = pulse_object._Units_inverse(pulse_object.biexciton_x_emission, 'meV')
            postion_xx_v = pulse_object._Units_inverse(pulse_object.biexciton_y_emission, 'meV')

        if plot:
            if pulse_intesity: 
                plot_pulse_h = np.abs(pulse_object.frequency_representation_x)**2
                plot_pulse_v = np.abs(pulse_object.frequency_representation_y)**2
            else:
                plot_pulse_h = np.abs(pulse_object.frequency_representation_x)
                plot_pulse_v = np.abs(pulse_object.frequency_representation_y)

            plt.figure()
            fig,ax_pulse = plt.subplots()
            ax_emission=ax_pulse.twinx() 
            ax_pulse.plot(plot_domain,plot_pulse_h,'b-',label='H pol')
            ax_pulse.plot(plot_domain,plot_pulse_v,'r-',label='V pol') 
            ax_pulse.set_xlabel(unit) 
            ax_pulse.set_ylim([0,np.max([plot_pulse_h,plot_pulse_v])*1.1])
            if pulse_intesity:
                ax_pulse.set_ylabel('Pulse intensity (a.u.)')
            else:
                ax_pulse.set_ylabel('Pulse amplitude (a.u.)')
            
            ax_emission.plot(np.array([1,1])*position_x_h,np.array([0,1])*exciton_x_emission,'b-',label='H pol')
            ax_emission.plot(np.array([1,1])*position_x_v,np.array([0,1])*exciton_y_emission,'r-',label='V pol')
            ax_emission.plot(np.array([1,1])*postion_xx_h,np.array([0,1])*biexciton_x_emission,'b-')
            ax_emission.plot(np.array([1,1])*postion_xx_v,np.array([0,1])*biexciton_y_emission,'r-')
            ax_emission.set_ylabel('Emission (a.u.)') 
            ax_emission.set_ylim([0,1])
        
    def open_control(self,pulse_object = None,open_gui = True,parent_window = None,previous_control = None):
        control_object = simc.simulator_control(self,pulse_object = pulse_object, open_gui = open_gui,previous_control= previous_control)
        return control_object
        
        # def open_control(self,pulse_object = None,simulation_object = None,open_gui = True,parent_window = None,previous_control = None):
        # control_object = smc.spectrometer_control(self,pulse_object = pulse_object,simulation_object = simulation_object, open_gui = open_gui,parent_window=parent_window, previous_control = previous_control)
        # return control_object



            
        pass

# fake devices for testing 

class fake_spectrometer:
    def __init__(self,start_wl= 785, end_wl = 803, n_wl = 1340) -> None:
        self.start_wl = start_wl
        self.end_wl = end_wl
        self.n_wl = n_wl
        pass

    def get_spectrum(self,data_file = None, pulse_object = None, draw = False, pandas = False):

        if data_file is not None:
            f = open(data_file, "r")
            lines = f.readlines()
            wavelength = []
            counts = []
            for line in lines:
                if pandas:
                    line = line.split(',')
                else:
                    line = line.split()
                #print(line)
                if line == ['\n']:
                    pass
                elif line[0] != 'wl':
                    wavelength.append(float(line[0]))
                    counts.append(float(line[1]))
            
        elif pulse_object is not None:
            wavelength = np.linspace(self.start_wl,self.end_wl,self.n_wl)
            pulse_wavelength = 299792.458/pulse_object.central_wavelength + pulse_object.frequencies
            
            pulse_wavelength = 299792.458/pulse_wavelength

            counts = np.interp(wavelength,pulse_wavelength,np.abs(pulse_object.frequency_representation_x)**2+
                            np.abs(pulse_object.frequency_representation_y)**2) 
            #counts = np.abs(60 + counts + np.random.normal(0, 10, len(wavelength)))
            counts = np.abs(counts +60)
        else:
            wavelength = np.linspace(self.start_wl,self.end_wl,self.n_wl)
            counts = np.abs(np.ones_like(wavelength)*60 + np.random.normal(0, 20, len(wavelength)))

        if draw:
            plt.figure()
            plt.plot(wavelength,counts)
            plt.xlabel('Wavelength (nm)')
            plt.ylabel('Counts')
            plt.show()

        return wavelength, counts 
    
    def close(self):
        print('Spectrometer closed')
        pass
    
class fake_motor:
    def __init__(self) -> None:
        pass

    def set_position(self,position):
        print('Moving to position: ',position)
        pass

    def close(self):
        print('Motor closed')
        pass

class fake_motor_alt:
    def __init__(self) -> None:
        pass

    def move_to(self,position):
        print('Moving to position: ',position)
        pass

    def close(self):
        print('Motor closed')
        pass

class fake_attenuator:
    def __init__(self) -> None:
        pass

    def set_attenuation(self,attenuation):
        print('Setting attenuation to: ',attenuation)
        pass

    def close(self):
        print('Attenuator closed')
        pass

class fake_power_meter:
    def __init__(self) -> None:
        pass

    def get_power(self,pulse_object = None):
        if pulse_object is not None:
            power = pulse_object.pulse_power
        else:
            power = 0
        print('got power:',power)
        return power
    
    def close(self):
        print('Power meter closed')
        pass

# extra functions 
    
def propagate_medium(pulse_object,length = 0, medium = 'air', correct_time_shift = False):
    # load csv file 'Sellmeier_coefficients.csv'
    # check if medium is in the file
    # calculate refractive index
    # calculate phase shift
    # add phase shift to pulse_object
    # return pulse_object
    # length in mm
    if type(pulse_object) is str:
        pulse_object = pg.load_pulse(pulse_object)
    pulse_object.clear_filter()
    pulse_object.add_filter_rectangle()

    K = [0,0,0]
    L = [0,0,0]
    with open('Sellmeier_coefficients.csv', newline='') as csvfile:
        data = csv.reader(csvfile, delimiter = ',')
        for row in data:
            if row[0] == medium.lower(): 
                K[0] = float(row[1])
                K[1] = float(row[2])
                K[2] = float(row[3])
                L[0] = float(row[4])
                L[1] = float(row[5])
                L[2] = float(row[6])
    
    length = length*1e6 # mm to nm
    refractive_index = _sellmeier(pulse_object.wavelengths,K,L)

    phase_shift = 2*np.pi*refractive_index*length/pulse_object.wavelengths

    if correct_time_shift: #not working correctly yet
        # shift_wl_ind = np.argmax(np.abs(pulse_object.frequency_representation_x) + np.abs(pulse_object.frequency_representation_y))
        # #shift_wl_ind = np.average(pulse_object.wavelengths,weights = np.abs(pulse_object.frequency_representation_x) + np.abs(pulse_object.frequency_representation_y))
        # phase_shift -= 2*np.pi*_sellmeier_d(pulse_object.wavelengths[shift_wl_ind],K,L)*length/pulse_object.wavelengths[shift_wl_ind]*pulse_object.frequencies
        old_max_time_ind = np.argmax(np.abs(pulse_object.temporal_representation_x) + np.abs(pulse_object.temporal_representation_y))
        old_max_time = pulse_object.time[old_max_time_ind]

    pulse_object._add_filter(np.exp(1j*phase_shift),pol='b',merging='*')

    pulse_object.apply_frequency_filter()

    if correct_time_shift:
        new_max_time_ind = np.argmax(np.abs(pulse_object.temporal_representation_x) + np.abs(pulse_object.temporal_representation_y))
        new_max_time = pulse_object.time[new_max_time_ind]
        time_shift = old_max_time - new_max_time
        pulse_object = pulse_time_delay(pulse_object,time_shift)

    return pulse_object


    


def polarising_beam_splitter(pulse_object_1,pulse_object_2 = None):
    #horizontal (x) gets transmitted, vertical (y) gets reflected
    if type(pulse_object_1) is str:
        pulse_object_1 = pg.load_pulse(pulse_object_1)
    if pulse_object_2 is not None:
        if type(pulse_object_2) is str:
            pulse_object_2 = pg.load_pulse(pulse_object_2)
    else:
        pulse_object_2 = pulse_object_1.copy_pulse()
        pulse_object_2.clear_all()

    pulse_object_1.clear_filter()
    pulse_object_2.clear_filter()

    pulse_out_1_1 = pulse_object_1.copy_pulse()
    pulse_out_1_1.add_filter_rectangle(transmission = 1,cap_transmission=False, polarisation = 'x')

    pulse_out_2_2 = pulse_object_2.copy_pulse()
    pulse_out_2_2.add_filter_rectangle(transmission = 1,cap_transmission=False, polarisation = 'x')

    pulse_out_1_2 = pulse_object_1.copy_pulse()
    pulse_out_1_2.add_filter_rectangle(transmission = 1,cap_transmission=False, polarisation = 'y')
    pulse_out_1_2.add_phase_filter(phase_taylor=[np.pi/2])

    pulse_out_2_1 = pulse_object_2.copy_pulse()
    pulse_out_2_1.add_filter_rectangle(transmission = 1,cap_transmission=False, polarisation = 'y')
    pulse_out_2_1.add_phase_filter(phase_taylor=[np.pi/2])

    pulse_out_1_1.apply_frequency_filter()
    pulse_out_2_2.apply_frequency_filter()
    pulse_out_1_2.apply_frequency_filter()
    pulse_out_2_1.apply_frequency_filter()

    pulse_out_1_1.merge_pulses(pulse_out_2_1)
    pulse_out_2_2.merge_pulses(pulse_out_1_2)

    pulse_out_1_1.clear_filter()
    pulse_out_2_2.clear_filter()
    return pulse_out_1_1, pulse_out_2_2
  

    

def beam_splitter(pulse_object_1, pulse_object_2 = None, transmission = 0.5, reflection = 0.5):
    if type(pulse_object_1) is str:
        pulse_object_1 = pg.load_pulse(pulse_object_1)
    if pulse_object_2 is not None:
        if type(pulse_object_2) is str:
            pulse_object_2 = pg.load_pulse(pulse_object_2)
    else:
        pulse_object_2 = pulse_object_1.copy_pulse()
        pulse_object_2.clear_all()
        
    pulse_object_1.clear_filter()
    pulse_object_2.clear_filter()
    
    pulse_out_1_1 = pulse_object_1.copy_pulse()
    pulse_out_2_2 = pulse_object_2.copy_pulse()
    pulse_out_1_2 = pulse_object_1.copy_pulse()
    pulse_out_2_1 = pulse_object_2.copy_pulse()



    pulse_out_1_1.add_filter_rectangle(transmission = np.sqrt(transmission),cap_transmission=True)
    pulse_out_2_2.add_filter_rectangle(transmission = np.sqrt(transmission),cap_transmission=True)
    pulse_out_1_2.add_filter_rectangle(transmission = np.sqrt(reflection),cap_transmission=True)
    pulse_out_2_1.add_filter_rectangle(transmission = np.sqrt(reflection),cap_transmission=True)
    
    pulse_out_1_2.add_phase_filter(phase_taylor=[np.pi/2])
    pulse_out_2_1.add_phase_filter(phase_taylor=[np.pi/2])

    pulse_out_1_1.apply_frequency_filter()
    pulse_out_2_2.apply_frequency_filter()
    pulse_out_1_2.apply_frequency_filter()
    pulse_out_2_1.apply_frequency_filter()

    pulse_out_1_1.merge_pulses(pulse_out_2_1)
    pulse_out_2_2.merge_pulses(pulse_out_1_2)

    pulse_out_1_1.clear_filter()
    pulse_out_2_2.clear_filter()
    
    return pulse_out_1_1, pulse_out_2_2
    
    
def phase_plate(pulse_object,phase):
    if type(pulse_object) is str:
        pulse_object = pg.load_pulse(pulse_object)
    pulse_object.clear_filter()
    pulse_object.add_filter_rectangle()
    pulse_object.add_phase_filter(phase_taylor=[phase])
    pulse_object.apply_frequency_filter()
    pulse_object.clear_filter()
    return pulse_object

def chirp_filter(pulse_object,chirp,central_wavelength = None):
    if type(pulse_object) is str:
        pulse_object = pg.load_pulse(pulse_object)
    if central_wavelength is None: 
        central_wavelength = pulse_object.central_wavelength
    pulse_object.clear_filter()
    pulse_object.add_filter_rectangle()
    pulse_object.add_phase_filter(unit='nm',central_f = central_wavelength, phase_taylor = [0,0,chirp])
    pulse_object.apply_frequency_filter()
    pulse_object.clear_filter()
    return pulse_object
        

def pulse_time_delay(pulse_object, time_delay, central_wavelength = None):
    if type(pulse_object) is str:
        pulse_object = pg.load_pulse(pulse_object)

    if central_wavelength is None: 
        central_wavelength = pulse_object.central_wavelength
    pulse_object.clear_filter()
    pulse_object.shift_in_time(time_delay)
    # pulse_object.add_filter_rectangle()
    # lin_phase = time_delay 
    # pulse_object.add_phase_filter(central_f = central_wavelength, phase_taylor=[0,lin_phase], polarisation = 'b',unit = 'nm',f_start = None, f_end = None)
    # pulse_object.apply_frequency_filter()
    # pulse_object.clear_filter()
    return pulse_object

def _interpolate_pulse_to_spectrum(pulse, spectrum_wavelength,polarisation = 'both',intensity = True):
    if intensity:
        power = 2
    else:
        power = 1
    if polarisation.lower()[0] == 'x':
        pulse_spectrum = np.abs(pulse.frequency_representation_x)**power
    elif polarisation.lower()[0] == 'y':
        pulse_spectrum = np.abs(pulse.frequency_representation_y)**power
    else:
        pulse_spectrum = np.abs(pulse.frequency_representation_x)**power + np.abs(pulse.frequency_representation_y)**power
    
    pulse_spectrum_interp = np.interp(spectrum_wavelength,pulse.wavelengths,pulse_spectrum)
    return pulse_spectrum_interp

def _check_method(obj,method):
    if obj is None:
        return
    if not hasattr(obj, method):
        print('ERROR: Method "'+method+'" not found!')
        print('Add method "'+method+'" to '+obj.__class__.__name__+' class.')
        return False
    else:
        return True

def save_device(obj, name,suffix = 'pu',experiment = None):
    if experiment is not None:
        name = experiment.device_path+'/'+name
    with open(name+'.pul_'+suffix, 'wb') as f:
        pickle.dump(obj, f)

def save_pulse(pulse_object, name='/pulse',experiment = None):
    pulse_object.save_pulse(save_dir = experiment.pulse_path, save_name = name)

def save_temp(pulse_object, name='/temp_pulse',experiment = None):
    pulse_object.save_pulse(save_dir = experiment.temp_path, save_name = name)

def load_pulse_device(file_name):
    if file_name[-7:] not in  ['.pul_mo',
     '.pul_ps',
     '.pul_at',
     '.pul_sp']:
        print('ERROR: File is not supported!')
        return
    with open(file_name, "rb") as f:
        device = pickle.load(f)
        device.print_info()
    return device

def excecute_folder(folder):
    
    files = os.listdir(folder)
    files = [os.path.join(folder, f) for f in files]
    files.sort(key=lambda x: os.path.getmtime(x))
    for file in files:
        print(file)
        if file[-7:] in  ['.pul_mo',
                         '.pul_ps',
                         '.pul_at',
                         '.pul_sp']:
            device = load_pulse_device(file)
            device.excecute()

def _sellmeier(lam, K, L):
    #lam in nm
    lam = lam* 1e-3
    n2 = 1
    for i in range(len(K)): 
        n2 += K[i]*lam**2/(lam**2 - L[i])
    return np.sqrt(n2)


def _sellmeier_d(lam, K, L):
    nom = 2*_sellmeier(lam, K, L)
    lam = lam * 1e-3
    term = 0
    for i in range(len(K)):
        term = term + 2*K[i]*lam/(lam**2-L[i])-2*K[i]*lam**3/(lam**2-L[i])**2
    dn = term/nom
    return dn #*1e-3


### control functions 

def update_previous_control(control_object):
    if control_object.previous_control is not None:
        control_object.previous_control.update_previous_control()
        
    
    

if __name__ == "__main__":
    print('I am main!')
    # test code
    #pulse = pg.PulseGenerator(0,200,0.02,central_wavelength=800)
    #pulse.add_filter_rectangle(unit = 'nm', central_f = 800, width_f=20, merging='+',transmission=1,cap_transmission=False)
    #pulse.add_filter_rectangle(unit = 'nm', central_f = 800, width