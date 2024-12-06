import numpy as np
from matplotlib import pyplot as plt
from scipy.io import savemat
import pyaceqd.pulsegenerator as pg
from pyaceqd.six_level_system.linear import sixls_linear
from pyaceqd.six_level_system.linear import energies_linear
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
import tkinter as tk
from tkinter import ttk
from matplotlib.backend_bases import key_press_handler
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,
                                               NavigationToolbar2Tk)
from matplotlib.figure import Figure

class attenuator_control:
    def __init__(self,attenuator_object, pulse_object = None, step_size = 0.02, large_step = 10, initial_attenuation = None, open_gui = True, parent_window = None, previous_control = None):
        
        self.attenuator_object = attenuator_object
        self.name = attenuator_object.name
        self.step_size = step_size
        self.large_step = large_step
        self.initial_attenuation = initial_attenuation
        self.open_gui = open_gui
        self.pulse_object_out = None
        self.parent_window = parent_window
        self.previous_control = previous_control
        
        if pulse_object is not None:
            if type(pulse_object) is str:
                 pulse_object = pg.load_pulse(pulse_object)
            self.pulse_object = pulse_object.copy_pulse()
            self.pulse_object.clear_filter()
            
        else:
            self.pulse_object = None
        pass
    
        if initial_attenuation is None:
            self.current_attenuation = 0.5
            
        if self.attenuator_object.device is None:
            self.excecute = False
        else:
            self.excecute = True
            
        # init# initialize the buttons
        self.button_step_increase = 0
        self.button_step_decrease = 0
        
        self.button_step_small_up = 0
        self.button_step_small_down = 0
        
        self.button_step_large_up = 0
        self.button_step_large_down = 0
        
        self.update_previous_control()
        self.update_attenuation()
        self.print_info()
        
        
        if self.open_gui:
            #print('hier könnte ihre GUI stehen!')
            self.gui()
        
    def reset_buttons(self):
        self.button_step_increase = 0
        self.button_step_decrease = 0
        
        self.button_step_small_up = 0
        self.button_step_small_down = 0
        
        self.button_step_large_up = 0
        self.button_step_large_down = 0
        
        
    def rescale_step_size(self):
        self.step_size = self.step_size + self.step_size*self.button_step_increase - 0.5*self.step_size*self.button_step_decrease
        
        if self.step_size < 0.001:
            self.step_size = 0.001
    
    def update_attenuation(self):
        #self.update_previous_control()
        self.rescale_step_size()
        self.current_attenuation = self.current_attenuation + self.step_size*self.button_step_small_up - self.step_size*self.button_step_small_down + self.step_size*self.button_step_large_up*self.large_step - self.step_size*self.button_step_large_down*self.large_step
        
        if self.current_attenuation < 0:
            self.current_attenuation = 0
        if self.current_attenuation > 1:
            self.current_attenuation = 1
        
        
        self.pulse_object_out = self.attenuator_object.set_attenuation(attenuation = self.current_attenuation, pulse_object = self.pulse_object,excecute = self.excecute)
    
    def update_control(self):
        self.update_attenuation()
        self.reset_buttons()
        
    def set_button_step_increase(self,botton_state = 1):
        self.button_step_increase = botton_state
        
    def set_button_step_decrease(self,botton_state = 1):
        self.button_step_decrease = botton_state
        
    def set_button_step_small_up(self,botton_state = 1):
        self.button_step_small_up = botton_state
    
    def set_button_step_small_down(self,botton_state = 1):
        self.button_step_small_down = botton_state
    
    def set_button_step_large_up(self,botton_state = 1):
        self.button_step_large_up = botton_state
    
    def set_button_step_large_down(self,botton_state = 1):
        self.button_step_large_down = botton_state
    
    def set_step_size(self,step_size):
        self.step_size = step_size
        
    def set_current_attenuation(self,attenuation):
        self.current_attenuation = attenuation
        
    def set_control_value(self,control_value):
        self.current_attenuation = control_value
    
    def get_control_value(self):
        return self.current_attenuation
    
    def get_current_attenuation(self):
        return self.current_attenuation
    
    def get_pulse_object(self):
        if self.pulse_object_out is None:
            return None
        return self.pulse_object_out.copy_pulse()
    
    # def update_previous_control(self):
    #     if self.previous_control is not None:
    #         if self.open_gui:
    #             #self.previous_control.update_gui()
    #             pass
    #         self.previous_control.update_control()
    #         self.pulse_object = self.previous_control.get_pulse_object()
    
    def update_previous_control(self):
        if self.previous_control is not None:
            if self.open_gui:
                #self.previous_control.update_gui()
                pass
            self.previous_control.update_previous_control()
            self.previous_control.update_control()
            self.pulse_object = self.previous_control.get_pulse_object()
            self.update_control()
        
    def print_info(self):
        print('Current attenuation: ', self.current_attenuation)
        print('Step size: ', self.step_size)
        print('Large step size: ', self.large_step)
        print('Initial attenuation: ', self.initial_attenuation)
        print('Open GUI: ', self.open_gui)
        print('Parent window: ', self.parent_window)
        print('Previous control: ', self.previous_control)
        
    def get_button_state(self):
        print('Attenuator buttons')
        print('Step increase: ', self.button_step_increase)
        print('Step decrease: ', self.button_step_decrease)
        print('Small step up: ', self.button_step_small_up)
        print('Small step down: ', self.button_step_small_down)
        print('Large step up: ', self.button_step_large_up)
    
    def close(self):
        self.attenuator_object.close()
    
    def update_gui(self):
        #self.button_press(button)
        self.set_step_size(float(self.step_size_entry.get()))
        self.update_attenuation()
        #self.print_info()
        #self.get_button_state()
        self.reset_buttons()
        self.current_attenuation_entry.delete(0, 'end')
        self.current_attenuation_entry.insert(0,str(np.round(self.current_attenuation,decimals=3)))
        
        self.step_size_entry.delete(0, 'end')
        self.step_size_entry.insert(0,str(np.round(self.step_size,decimals=3)))
        
        self.attenuation_bar['value'] = self.current_attenuation
        
        
    def gui(self):
        self.open_gui = True
        if self.previous_control is None:
            self.gui_window = tk.Tk()
        else:
            self.gui_window = tk.Toplevel(self.previous_control.gui_window)
            
        def close_gui():
            self.gui_window.destroy()
            if self.attenuator_object.device is not None:
                self.attenuator_object.close()
            
            print('Window closed')
        
        def apply_attenuation():
            self.set_current_attenuation(float(self.current_attenuation_entry.get()))
            self.update_gui()

        # # update the gui after 100 ms to catch the
        # button pres
            
        self.gui_window.protocol("WM_DELETE_WINDOW", close_gui)
        self.gui_window.geometry('500x100')
        self.gui_window.title('Attenuator control: '+ self.attenuator_object.name)
        
        tk.Button(self.gui_window, text = '<<', width= 10, command=lambda: [self.set_button_step_large_down(), self.update_gui()]).grid(row=1,column=0)
        
        tk.Button(self.gui_window, text = '<', width= 10, command=lambda: [self.set_button_step_small_down(), self.update_gui()]).grid(row=1,column=1)
        
        tk.Button(self.gui_window, text = '>', width= 10, command=lambda: [self.set_button_step_small_up(), self.update_gui()]).grid(row=1,column=2)
        
        tk.Button(self.gui_window, text = '>>', width= 10, command=lambda: [self.set_button_step_large_up(), self.update_gui()]).grid(row=1,column=3)
        
        tk.Button(self.gui_window, text = 'Set attenuation:', command= apply_attenuation).grid(row=2,column=0)
        
        self.current_attenuation_entry = tk.Entry(master= self.gui_window, width=10)
        self.current_attenuation_entry.insert(0,str(np.round(self.current_attenuation,decimals=3)))
        self.current_attenuation_entry.grid(row=2,column=1)
        
        tk.Button(self.gui_window, text='step +', width= 10, command=lambda: [self.set_button_step_increase(), self.update_gui()]).grid(row=2,column=2)
        
        tk.Button(self.gui_window, text='step -', width= 10, command=lambda: [self.set_button_step_decrease(), self.update_gui()]).grid(row=3,column=2)
        
        self.step_size_entry = tk.Entry(master= self.gui_window, width=10)
        self.step_size_entry.insert(0,str(np.round(self.step_size,decimals=3)))
        self.step_size_entry.grid(row=2,column=3)
        
        self.attenuation_bar = ttk.Progressbar(master=self.gui_window, orient='horizontal', length=200, mode='determinate', maximum=1, value=self.current_attenuation)
        self.attenuation_bar.grid(row=3,column=0)
        
    def start_gui(self):
        self.gui_window.mainloop()
        
if __name__ == '__main__':
    
    
    import os as os
    from Pulse_v2 import fake_spectrometer, fake_motor, pulse_shaper_obj, simulator, attenuator, half_wave_plate, fake_attenuator, motor, spectrometer, load_pulse_device, create_experiment, save_pulse, save_temp, save_device, excecute_folder, power_meter
    
    cur_dir = os.getcwd()
    if 'PulseGenerationACE' in cur_dir:
        pass
    else:
    # change the directory to the folder where the calibration files are stored
        os.chdir(cur_dir+'/PulseGenerationACE/PULSE')
    print(cur_dir) 
    
    
    qd_calibration = 'QD_Iker_April.txt'
    ps_calibration2 = 'calibration_slit_38.txt' 
    ps_calibration = 'calibration_Cam_230424.txt' 
    
    t_0 = 0 
    t_end = 100
    
    initial_pulse = pg.PulseGenerator(t_0,t_end,0.2,calibration_file=qd_calibration)
    initial_pulse.add_gaussian_time(width_t=0.5,t0=t_end/2)
    #initial_pulse.plot_pulses(domain='nm')
    
    lab_motor = motor(fake_motor(),name='motor1')
    lab_motor2 = motor(fake_motor(),name='motor2')
    
    pulse_shaper = pulse_shaper_obj(device=lab_motor, calibration_file=ps_calibration, name='heisl')
    ps_controller = pulse_shaper.open_control(pulse_object=initial_pulse,previous_control = None)
    
    pulse_shaper2 = pulse_shaper_obj(device=lab_motor2, calibration_file=ps_calibration, name='bert')
    ps_controller2 = pulse_shaper2.open_control(pulse_object = initial_pulse,previous_control = None)
    
    lab_att = fake_attenuator()
    lab_att2 = fake_attenuator()
    
    att_object = attenuator(lab_att, name = 'F Att')
    att_controller = att_object.open_control(previous_control=ps_controller)
    
    att_object2 = attenuator(lab_att2, name = 'U Att')
    att_controller2 = att_object2.open_control(previous_control=ps_controller2)
    
    
    #pulse_shaper3 = pulse_shaper_obj(device=lab_motor, calibration_file=ps_calibration, name='gregor')
    #ps_controller3 = pulse_shaper3.open_control(pulse_object = initial_pulse,previous_control = None) 
    
    lab_spectromter = fake_spectrometer(start_wl=774, end_wl = 790, n_wl = 1024)

    spec = spectrometer(device=lab_spectromter, name='MF_spec') 
    sm_controller = spec.open_control(pulse_object=None,previous_control = [att_controller, att_controller2], open_gui= True) #, ps_controller3#ps_controller2.get_pulse_object()
    #sm_controller.toggle_display_pulse()
    #sm_controller.toggle_display_experiment()
    #sm_controller.gui()
    
    sm_controller.start_gui()
    
    
    
    

    pass
        