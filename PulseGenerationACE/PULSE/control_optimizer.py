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

from scipy.optimize import minimize
from skopt import gp_minimize


class control_optimizer():
    def __init__(self, device_control: list = [], measururement_control: list = [], scan_limits = [], measurement_kind = 'spectrometer', num_calls = 20, num_random_starts = 10, open_gui = True, active_measurement_index = 0, name = 'control_optimizer') -> None: 
        
        if type(device_control) != list:
            device_control = [device_control]
            
        if type(measururement_control) != list:
            measururement_control = [measururement_control]
            
        self.device_control_initial = device_control
        self.measururement_control = measururement_control
        
        
        
        # for i in range(len(device_control)):
        #     if len(scan_limits) < i: 
        #         scan_limits.append((float(0.0),float(1.0)))
        #self.scan_limits_initial = scan_limits
        self.name = name
        
        self.device_participating = []
        
        self.snipping_window = []
        self.snipping_bool = []
        for i in range(len(device_control)):
            self.device_participating.append(True)
            self.snipping_window.append(float(0.5))
            self.snipping_bool.append(True)
        
        self.set_scan_limits(scan_limits)
        self.set_participating()
        self.measurement_kind = measurement_kind
        
        self.num_calls = num_calls
        self.num_random_starts = num_random_starts
        #print(scan_limits)
        
        self.open_gui = open_gui
        self.run_counter = 0
        self.active_measurement_index = active_measurement_index
        
        if self.measurement_kind == 'spectrometer':
            #self.set_spectrometer_measurement()
            pass
        
        if self.open_gui:
            self.gui()
    
    def set_device_participating(self,device_index,participating = True):
        self.device_participating[device_index] = participating
    
    def set_participating(self):
        self.device_control = []
        self.scan_limits = []
        for i, device in enumerate(self.device_control_initial):
            if self.device_participating[i]:
                self.device_control.append(device)  
                self.scan_limits.append(self.scan_limits_initial[i])  
                
    def set_active_measurement_index(self,active_measurement_index):
        self.active_measurement_index = active_measurement_index
                
    def set_snipping(self):
        for i, device in enumerate(self.device_control_initial):
            if self.snipping_bool[i]:
                control_value = device.get_control_value()
                self.scan_limits_initial[i] = (control_value - self.snipping_window[i],control_value + self.snipping_window[i])
                #print('snipped')
            elif self.open_gui:
                self.scan_limits_initial[i] = (float(self.scan_limit_entries_min[i].get()),float(self.scan_limit_entries_max[i].get()))
        #print(self.scan_limits_initial)
        
    def set_scan_limits(self,scan_limits=[]):
        self.scan_limits_initial = scan_limits
        if len(scan_limits) < len(self.device_control_initial):
            for i in range(len(self.device_control_initial) - len(scan_limits)):
                self.scan_limits_initial.append((float(0.0),float(1)))
                
        self.scan_limits = self.scan_limits_initial
    def update_device_control(self):
        for device in self.device_control:
            device.update_control()
            #print(device.name)
            if device.open_gui:
                device.update_gui()
                #print('gui updated')
                #device.gui_window.update_idletasks()
                #device.gui_window.update()
    
    def set_control_value(self,control_values):
        for i, device in enumerate(self.device_control):
            device.set_control_value(control_values[i])
            #print('Control Value: ',control_values[i])
            #print('Device Value: ',device.get_control_value())
    
    def set_num_calls(self,num_calls):
        self.num_calls = num_calls
        
    def set_num_random_starts(self,num_random_starts):
        self.num_random_starts = num_random_starts
    
    def get_control_value(self):
        control_values = []
        for device in self.device_control:
            control_values.append(device.get_control_value())
        return control_values
            
    def set_measurement_kind(self,measurement_kind):
        self.measurement_kind = measurement_kind
    
    def update_measurement_control(self):
        for device in self.measururement_control:
            if device.open_gui:
                device.gui_window.update_idletasks()
                device.update_gui()
            
            #device.update_measurement()
            
    def set_spectrometer_measurement(self, center_wl = 0, width_wl = 0, method = 'max'):
        #self.spectrometer_center_wl = center_wl
        #self.spectrometer_width_wl = width_wl
        #self.spectrometer_weight = weight
        self.arguments = [center_wl,width_wl,method]
        
    def run_optimization(self):
        self.run_counter = 0
        #measurement_gui_state = self.measururement_control[0].force_gui_update
        self.measururement_control[0].force_gui_update = False
        self.stop_flag = False
        self.optimized_control_values = gp_minimize(self.optimize_function, n_calls = self.num_calls, n_random_starts = self.num_random_starts, dimensions=self.scan_limits, x0=self.get_control_value(),y0=self.optimize_function(self.get_control_value()),initial_point_generator='halton')
        
        self.optimize_function(self.optimized_control_values.x)
        self.measururement_control[0].force_gui_update = True #measurement_gui_state
        self.measururement_control[0].update_gui()
        pass
    
    def optimize_function(self, optimize_values:list): 
        if self.stop_flag:
            return 0
        
        self.set_control_value(optimize_values)
        self.update_device_control()
        self.update_measurement_control()
        
        self.run_counter += 1
        
        if self.open_gui:
            self.progress_bar['value'] = self.run_counter/self.num_calls
            self.progress_bar.update()
            
            self.run_counter_label['text'] = 'Run Counter: '+str(self.run_counter-1)
            self.run_counter_label.update()
        
        if self.measurement_kind == 'spectrometer':
            #output = self.spectrometer_measurement()
            output = self.measururement_control[0].measurement_method(self.arguments) # only works with gui open ... hmmm
            
            if type(output) == list:
                print('Output: ',sum(output))
                return -sum(output)
            else:
                print('Output: ',output)
                return -output
        pass
    
    
    # def spectrometer_measurement(self):
    #     #intensity_vector = []
    #     #for device in self.measururement_control:
    #     device = self.measururement_control[self.active_measurement_index]
    #     if device.display_experiment:
    #         wavelength = device.get_wavelength_experiment()
    #         intensity = device.get_intensity_experiment()
    #     elif any([device.display_pulse,device.display_simulation]):
    #         wavelength = device.get_wavelength_pulse()
    #         intensity = device.get_combined_pulse_simulation()
            
    #     if self.spectrometer_width_wl == 0:
    #         wavelength_window = np.amin(np.abs(wavelength - self.spectrometer_center_wl))  
    #     else:  
    #         wavelength_window = np.where(np.abs(wavelength - self.spectrometer_center_wl) < self.spectrometer_width_wl/2)
            
    #     #print(np.mean(wavelength[wavelength_window]))
        
    #     if self.spectrometer_weight == 'max':
    #         intensity = np.max(intensity[wavelength_window])
    #     if self.spectrometer_weight == 'mean':
    #         intensity = np.mean(intensity[wavelength_window])
    #     if self.spectrometer_weight == 'sum':
    #         intensity = np.sum(intensity[wavelength_window])
    #     #    intensity_vector.append(intensity)
    #     #self.intensity_vector = intensity_vector
    #     return intensity #np.sum(intensity_vector)
    
    def run(self):
        self.run_button.config(state=tk.DISABLED)
        for i, device in enumerate(self.device_control_initial):
            self.device_participating[i] = self.scan_participating_var[i].get()
        
            
        for i, device in enumerate(self.device_control_initial):
            self.scan_limits_initial[i] = (float(self.scan_limit_entries_min[i].get()),float(self.scan_limit_entries_max[i].get()))
        
        self.set_participating()
        self.num_calls = int(self.num_calls_entry.get())
        self.num_random_starts = int(self.num_random_starts_entry.get())
        
        if self.measurement_kind == 'spectrometer':
            #self.spectrometer_center_wl = float(self.center_wl_entry.get())
            #self.spectrometer_width_wl = float(self.width_wl_entry.get())
            
            self.arguments = self.measururement_control[0].get_measurement_args(control = self)
            # for i in range(len(self.measurement_args)):
            #     if self.measurement_args[i] is None:
            #         self.arguments.append(None)
            #     else:
            #         self.arguments.append(self.measurement_args[i].get())

        self.stop_button.config(state=tk.NORMAL)
        self.run_optimization()
        self.stop_button.config(state=tk.DISABLED)
        self.run_button.config(state=tk.NORMAL)
        self.button_color(self.run_button,False)
    
    def stop_run(self):
        self.stop_flag = True
        self.stop_button.config(state=tk.DISABLED)
        self.run_button.config(state=tk.NORMAL)
        self.button_color(self.run_button,False)
        
    def button_color(self, button_object, state):
            if state:
                button_object.config(bg='green')
            else:
                button_object.config(bg='white')
    
    def gui(self):
        
        self.open_gui = True
        if self.measururement_control[0].open_gui:
            self.gui_window = tk.Toplevel(self.measururement_control[0].gui_window)
        #self.gui_window = tk.Toplevel(self.measururement_control[0].gui_window)
        else:
            self.gui_window = tk.Tk()
            
        def close_gui():
            self.gui_window.destroy()
            print('Optimiser closed')
        
        self.gui_window.protocol("WM_DELETE_WINDOW", close_gui)
        self.gui_window.title('Control Optimization')
        
        self.set_scan_limits()
        
        self.run_button = tk.Button(self.gui_window, text='Run Optimization', command=self.run)
        self.run_button.grid(row=0,column=0)
        self.run_button.config(command= lambda: [self.button_color(self.run_button,True),self.run()])
        self.button_color(self.run_button,False)
        
        self.progress_bar = ttk.Progressbar(self.gui_window, orient='horizontal', length=150, mode='determinate', value=self.run_counter/self.num_calls, maximum=1)
        self.progress_bar.grid(row=0,column=1,columnspan=1)
        
        self.run_counter_label = tk.Label(self.gui_window, text='Run Counter: '+str(self.run_counter))
        self.run_counter_label.grid(row=0,column=2)
        
        self.stop_button = tk.Button(self.gui_window, text='Stop', command=self.stop_run, state=tk.DISABLED)
        self.stop_button.grid(row=0,column=3)
        
        tk.Label(self.gui_window, text='Number of Calls').grid(row=1,column=0)
        self.num_calls_entry = tk.Entry(self.gui_window)
        self.num_calls_entry.grid(row=1,column=1)
        self.num_calls_entry.insert(0,str(self.num_calls))
        
        tk.Label(self.gui_window, text='Number of Random Starts').grid(row=2,column=0)
        self.num_random_starts_entry = tk.Entry(self.gui_window)
        self.num_random_starts_entry.grid(row=2,column=1)
        self.num_random_starts_entry.insert(0,str(self.num_random_starts))
        
        tk.Button(self.gui_window, text='Snipp limits', command=self.update_gui).grid(row=2,column=4)
        
        tk.Label(self.gui_window, text='Snipping Window').grid(row=2,column=5)
        
        
        row_offset = 3
        self.scan_limit_entries_min = []
        self.scan_limit_entries_max = []
        self.scan_participating_checkbox = []
        self.scan_participating_var = []
        self.snipping_checkbox = []
        self.snipping_window_entries = []
        self.snipping_bool_var = []
        
        for i, device in enumerate(self.device_control_initial):
            tk.Label(self.gui_window, text=device.name).grid(row=i+row_offset,column=0)
            
            self.scan_limit_entries_min.append(tk.Entry(self.gui_window))
            self.scan_limit_entries_min[i].grid(row=i+row_offset,column=1)
            self.scan_limit_entries_min[i].insert(0,str(self.scan_limits_initial[i][0]))
            
            self.scan_limit_entries_max.append(tk.Entry(self.gui_window))
            self.scan_limit_entries_max[i].grid(row=i+row_offset,column=2)
            self.scan_limit_entries_max[i].insert(0,str(self.scan_limits_initial[i][1]))
            
            self.scan_participating_var.append(tk.BooleanVar(value=self.device_participating[i]))
            self.scan_participating_checkbox.append(tk.Checkbutton(self.gui_window, text='Participating', onvalue=True, offvalue=False, variable= self.scan_participating_var[i]))
            self.scan_participating_checkbox[i].grid(row=i+row_offset,column=3)
            
            self.snipping_bool_var.append(tk.BooleanVar(value=self.snipping_bool[i]))
            self.snipping_checkbox.append(tk.Checkbutton(self.gui_window, text='Snipping', onvalue=True, offvalue=False, variable= self.snipping_bool_var[i]))
            self.snipping_checkbox[i].grid(row=i+row_offset,column=4)
            
            self.snipping_window_entries.append(tk.Entry(self.gui_window))
            self.snipping_window_entries[i].grid(row=i+row_offset,column=5)
            self.snipping_window_entries[i].insert(0,str(self.snipping_window[i]))
        
        
        
        if self.measurement_kind == 'spectrometer': 
            # tk.Label(self.gui_window, text='Spectrometer Measurement').grid(row=i+row_offset+1,column=0)
            # tk.Label(self.gui_window, text='Center Wavelength').grid(row=i+row_offset+2,column=0)
            # self.center_wl_entry = tk.Entry(self.gui_window)
            # self.center_wl_entry.grid(row=i+row_offset+2,column=1)
            # self.center_wl_entry.insert(0,str(self.spectrometer_center_wl))
            
            # tk.Label(self.gui_window, text='Width Wavelength').grid(row=i+row_offset+3,column=0)
            # self.width_wl_entry = tk.Entry(self.gui_window)
            # self.width_wl_entry.grid(row=i+row_offset+3,column=1)
            # self.width_wl_entry.insert(0,str(self.spectrometer_width_wl))
            
            #self.measururement_control[0].single_line_measurement_gui(control = self, row_offset = #row_offset +1 +i)
            
            self.measururement_control[0].measurement_method_gui(control = self, row_offset = row_offset + row_offset +1 +i)
        
        self.update_gui()
        
    def update_gui(self):
        
        for i in range(len(self.device_control_initial)):
            self.snipping_window[i] = float(self.snipping_window_entries[i].get())
            self.snipping_bool[i] = self.snipping_bool_var[i].get()
        
        self.set_snipping()
        self.set_scan_limits(self.scan_limits_initial)
        for i, device in enumerate(self.device_control_initial):
            self.scan_limit_entries_min[i].delete(0,tk.END)
            self.scan_limit_entries_min[i].insert(0,str(self.scan_limits_initial[i][0])) 
            
            self.scan_limit_entries_max[i].delete(0,tk.END)
            self.scan_limit_entries_max[i].insert(0,str(self.scan_limits_initial[i][1]))
        
        #self.gui_window.update_idletasks()
        self.gui_window.update()
        
    def start_gui(self):
        self.gui_window.mainloop()
            
            
            
        

if __name__ == '__main__':
    import numpy as np
    from matplotlib import pyplot as plt
    from scipy.io import savemat
    import pyaceqd.pulsegenerator as pg 
    import os as os
    from Pulse_v2 import fake_spectrometer, fake_motor, pulse_shaper_obj, simulator, attenuator, half_wave_plate, fake_attenuator, motor, spectrometer, load_pulse_device, create_experiment, save_pulse, save_temp, save_device, excecute_folder, power_meter, time_delay, fake_power_meter  
    
    
    cur_dir = os.getcwd()
    if 'PulseGenerationACE' in cur_dir:
        pass
    else:
    # change the directory to the folder where the calibration files are stored
        os.chdir(cur_dir+'/PulseGenerationACE/PULSE')
    print(cur_dir) 
    
    ps_calibration = 'calibration_slit_13.txt' 
    #ps_calibration2 = 'calibration_slit_38.txt' 
    qd_calibration = 'QD_Iker_April.txt'

    lab_motor = motor(fake_motor(),name='motor1')
    lab_att = fake_attenuator()
    lab_spectromter = fake_spectrometer(start_wl=774.15, end_wl = 793.422, n_wl = 1340)
    
    lab_motor2 = motor(fake_motor(),name='motor2')
    lab_motor3 = motor(fake_motor(),name='motor3')
    
    lab_att2 = fake_attenuator()

    t_0 = 0 
    t_end = 100

    initial_pulse = pg.PulseGenerator(t_0,t_end,0.2,calibration_file=qd_calibration)
    initial_pulse.add_gaussian_time(unit = 'nm', central_f= 779.89, width_t=0.1,t0=30,area_time=20, polarisation=[1,0])

    pulse_shaper = pulse_shaper_obj(device=lab_motor, calibration_file=ps_calibration, name='pulse_shaper_A')
    att_object = attenuator(lab_att, name = 'attenuator_A')
    
    pulse_shaper2 = pulse_shaper_obj(device=lab_motor2, calibration_file=ps_calibration, name='heisl2')
    att_object2 = attenuator(lab_att, name = 'F Att2')
    
    delay_stage = time_delay(device=lab_motor3, name='delay_stage')
    
    
    
    ps_controller = pulse_shaper.open_control(pulse_object=initial_pulse)
    att_controller = att_object.open_control(previous_control=ps_controller)
    
    ps_controller2 = pulse_shaper2.open_control(pulse_object=initial_pulse)
    att_controller2 = att_object2.open_control(previous_control=ps_controller2)
    delay_controller = delay_stage.open_control(previous_control=att_controller2)
    
    simulator_object = simulator(qd_calibration=qd_calibration, sim_kind='ace', temp_dir='sim_dump/')
    sim_controller = simulator_object.open_control(previous_control = [att_controller], open_gui=False) 
    
    spec = spectrometer(device=lab_spectromter, name='MF_spec') 
    sm_controller = spec.open_control(pulse_object=None,previous_control = [att_controller], simulation_control= sim_controller,open_gui=False) #, ps_controller3#ps_controller2.get_pulse_object()
    sm_controller.set_simulation_background(45)
    sm_controller.set_simulation_gaussian_noise(10)#4.35
    sm_controller.set_pulse_scale(1)
    sm_controller.set_simulation_counts(1000)
    sm_controller.toggle_running()
    sm_controller.set_measurement_arguments(arguments=[779.94, 0.05]) 
    
    sim_controller.toggle_running()
    sim_controller.gui()
    
    
    sm_controller.change_view()
    sm_controller.gui()
    
    
    lab_power_meter = fake_power_meter()
    
    power_meter_A = power_meter(device=lab_power_meter,name='power_meter_A')
    pm_controller = power_meter_A.open_control(previous_control=[att_controller])
    
    
    co = control_optimizer(device_control=[ps_controller,att_controller],measururement_control=sm_controller,measurement_kind='spectrometer',open_gui=False) # sm_controller
    #co.set_scan_limits([(780.0,781.0),(0.0,0.1)])
    #co.set_spectrometer_measurement(779.89,0.2,method='max')
    co.gui()
    #co.run_optimization()
    
    
    
    sm_controller.start_gui()
    # 