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
from matplotlib.backend_bases import key_press_handler
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,
                                               NavigationToolbar2Tk)
from matplotlib.figure import Figure 

from datetime import datetime

class powermeter_control:
    def __init__(self, powermeter_object, pulse_object = None, previous_control = [None], open_gui = True, record_time = 120, refresh_rate = 500):
        
        self.powermeter_object = powermeter_object
        self.name = powermeter_object.name
        self.pulse_object = pulse_object
        
        self.previous_control = previous_control
        
        self.display_experiment = False
        self.display_pulse = False
        
        self.force_gui_update = True
        
        if pulse_object is not None:
            if type(pulse_object) is str:
                 pulse_object = pg.load_pulse(pulse_object)
            out_pulse_object = pulse_object.copy_pulse()
            out_pulse_object.clear_filter()
            self.pulse_object = out_pulse_object
        else:
            self.pulse_object = None
        pass
        
        self.open_gui = open_gui
        
        if self.powermeter_object.device is not None:
            self.excecute = True
        else:
            self.excecute = False
        
        self.record_time = record_time
        self.refresh_rate = refresh_rate
        
        self.running_experiment = False
        self.running_pulse = False
        
        self.reset_recording()
        self.set_measurement_method()
        
        if self.open_gui:
            self.gui()
    
    
    def reset_recording(self):
        self.power_experiment_vector = []
        self.power_pulse_vector = []
        self.time_vector = []
        self.time_vector_seconds = []
        
        self.max_power_experiment = 0
        self.max_power_pulse = 0
        
        self.current_power_experiment = 0
        self.current_power_pulse = 0
    
    def update_previous_control(self):
        # if self.previous_control is not None:
        #     self.pulse_object = self.previous_control.get_pulse_object()
        
        if self.previous_control is not None and type(self.previous_control) is list:
            pulse_list = []
            for j in range(len(self.previous_control)):
                if self.previous_control[j].open_gui:
                    #self.previous_control[j].update_gui()
                    pass
                self.previous_control[j].update_previous_control()
                pulse_list.append(self.previous_control[j].get_pulse_object())
                #print('shaper: '+ self.previous_control[j].pulse_shaper_object.name)
                #print(self.previous_control[j].get_pulse_object().pulse_power)
            
            if pulse_list[0] is not None:
                self.pulse_object = pulse_list[0].copy_pulse()
            
                for i in range(len(pulse_list)-1): 
                    #print(i+1)
                    self.pulse_object.merge_pulses(pulse_list[i+1])
                
        elif self.previous_control is not None:
            if self.open_gui:
                #self.previous_control.update_gui()
                pass
            self.previous_control.update_previous_control()
            self.pulse_object = self.previous_control.get_pulse_object()
            
    def toggle_display_experiment(self):
        self.display_experiment = not self.display_experiment
        
    def toggle_display_pulse(self):
        self.display_pulse = not self.display_pulse
        
    def toggle_running_experiment(self):  
        self.running_experiment = not self.running_experiment
        
    def toggle_running_pulse(self):
        self.running_pulse = not self.running_pulse
         
    def set_record_time(self, record_time):
        self.record_time = record_time
        
    def set_refresh_rate(self, refresh_rate):
        self.refresh_rate = refresh_rate
    
    def update_measurement(self):
        self.update_previous_control()
        current_power = self.powermeter_object.get_power(excecute = self.excecute, pulse_object = self.pulse_object)
        self.current_power_experiment = current_power[0]
        self.current_power_pulse = current_power[1]
        
        if self.current_power_experiment > self.max_power_experiment:
            self.max_power_experiment = self.current_power_experiment
        
        if self.pulse_object is not None:
            if self.current_power_pulse > self.max_power_pulse:
                self.max_power_pulse = self.current_power_pulse
        
        #self.power_experiment_vector.append(self.current_power_experiment)
        #self.power_pulse_vector.append(self.current_power_pulse)
        #self.time_vector.append(datetime.now())
        #self.update_recording()
        
    def simple_power_measurement(self, arguments = []):
        self.update_measurement()
        
        if self.running_experiment: 
            self.measurement = self.current_power_experiment
        elif self.running_pulse:
            self.measurement = self.current_power_pulse
        else:
            self.measurement = 0
            
        return self.measurement
    
    def get_measurement_result(self):
        return self.measurement
    
    def simple_power_measurement_gui(self,control = None, row_offset = 0):
        control.measurement_args_pm = [None]
        num_rows = 2
        tk.Label(control.gui_window, text=self.name + ': Simply power!').grid(row=1+row_offset, column=0)
        tk.Label(control.gui_window, text='Nothing to input, just power!').grid(row=2+row_offset, column=0)
        return num_rows
        
    def set_measurement_method(self, method= 'simple'):
        if method == 'simple':
            self.measurement_method = self.simple_power_measurement
            self.measurement_method_gui = self.simple_power_measurement_gui
            self.number_measurement_outputs = 1
    
    def get_number_measurement_outputs(self):
        return self.number_measurement_outputs
    
    def get_measurement_args(self,control):
        return control.measurement_args_pm
    
    def get_power_experiment(self):
        return self.power_experiment_vector
    
    def get_power_pulse(self):
        return self.power_pulse_vector
    
    def get_time_vector_seconds(self):
        return self.time_vector_seconds
    
    def get_time_vector(self):
        return self.time_vector
    
    def update_recording(self):
        current_time = datetime.now()
        
        self.time_vector_seconds = []
        for i in reversed(range(len(self.time_vector))):
            self.time_vector_seconds.append((current_time - self.time_vector[i]).total_seconds())
            if (current_time - self.time_vector[i]).total_seconds() > self.record_time:
                self.power_experiment_vector.pop(i)
                self.power_pulse_vector.pop(i)
                self.time_vector.pop(i)
                self.time_vector_seconds.pop(i)
        self.time_vector_seconds = list(reversed(self.time_vector_seconds))
    def close(self):
        self.powermeter_object.close()
        
    def button_color(self, button_object, state):
            if state:
                button_object.config(bg='green')
            else:
                button_object.config(bg='white')
    
    def gui(self):
        self.update_measurement()
        self.open_gui = True
        
        if self.previous_control is None:
            self.gui_window = tk.Tk()
        
        else:
            if type(self.previous_control) is not list:
                self.previous_control = [self.previous_control] 
                
            if self.previous_control[0].open_gui:
                self.gui_window = tk.Toplevel(self.previous_control[0].gui_window)
            else:
                self.gui_window = tk.Tk()    
        pass
        
        def close_gui():
            if self.powermeter_object.device is not None:
                self.close()
            self.gui_window.destroy()
            print(self.name + ' closed')
        
        self.gui_window.protocol("WM_DELETE_WINDOW", close_gui)
        self.gui_window.title('Powermeter Control: '+self.name)
        
        # fig_powermeter = Figure(figsize=(5, 4), dpi=100)
        # self.ax_experiment = fig_powermeter.add_subplot()
        # self.ax_pulse = self.ax_experiment.twinx()

        # self.ax_experiment.set_xlabel('Time / s')
        # self.ax_experiment.set_ylabel('Power (mW)')
        # self.ax_pulse.set_ylabel('Power (a.u.)')
        
        # self.ax_experiment.set_xlim([0, self.record_time])
        
        # self.ax_experiment_plot, = self.ax_experiment.plot(self.time_vector_seconds,self.power_experiment_vector, 'k-')
        # self.ax_pulse_plot, = self.ax_pulse.plot(self.time_vector_seconds,self.power_pulse_vector, 'r-')
        
        # self.canvas = FigureCanvasTkAgg(fig_powermeter, master=self.gui_window)  # A tk.DrawingArea.
        # self.canvas.draw()
        
        
        # self.canvas.get_tk_widget().grid(row=0, column=0, columnspan=3)
        
        self.experiment_label = tk.Label(self.gui_window, text='Experiment: ')
        self.experiment_label.grid(row=0, column=1)
        
        self.pulse_label = tk.Label(self.gui_window, text='Simulation: ')
        self.pulse_label.grid(row=0, column=2)
        
        current_pow_exp_label = tk.Label(self.gui_window, text='Current Power: ')
        current_pow_exp_label.grid(row=1, column=0)
        
        
        self.current_pow_exp_var = tk.Label(self.gui_window, text=str(self.current_power_experiment))
        self.current_pow_exp_var.grid(row=1, column=1)
        
        self.current_pow_pul_var = tk.Label(self.gui_window, text=str(self.current_power_pulse))
        self.current_pow_pul_var.grid(row=1, column=2)
        
        maximum_pow_exp_label = tk.Label(self.gui_window, text='Max Power: ')
        maximum_pow_exp_label.grid(row=2, column=0)
        
        self.maximum_pow_exp_var = tk.Label(self.gui_window, text=str(self.max_power_experiment))
        self.maximum_pow_exp_var.grid(row=2, column=1)
        
        self.maximum_pow_pul_var = tk.Label(self.gui_window, text=str(self.max_power_pulse))
        self.maximum_pow_pul_var.grid(row=2, column=2) 
        
        self.reset_button = tk.Button(self.gui_window, text='Reset', command=self.reset_recording)
        self.reset_button.grid(row=3, column=0)
        
        self.running_experiment_button = tk.Button(self.gui_window, text='Run Experiment', command= lambda: [self.toggle_running_experiment(), self.button_color(self.running_experiment_button, self.running_experiment), self.update_gui()])
        self.running_experiment_button.grid(row=4, column=1)
        self.button_color(self.running_experiment_button, self.running_experiment)
        
        self.running_pulse_button = tk.Button(self.gui_window, text='Run Pulse', command= lambda: [self.toggle_running_pulse(), self.button_color(self.running_pulse_button, self.running_pulse), self.update_gui()])
        self.running_pulse_button.grid(row=4, column=2)
        self.button_color(self.running_pulse_button, self.running_pulse)
        
        if self.pulse_object is None:
            self.running_pulse_button.config(state='disabled')
        
        self.update_gui()
    
    def update_gui(self):
        self.gui_window.update_idletasks()
        if any([self.running_experiment, self.running_pulse]):
            self.update_measurement()
        if self.running_experiment:
           # self.ax_experiment_plot.set_xdata(self.get_time_vector_seconds())
           # self.ax_experiment_plot.set_ydata(self.get_power_experiment())
           # self.ax_experiment.set_ylim([min(self.get_power_experiment())*0.9, max(self.get_power_experiment())*1.1])
            
            self.current_pow_exp_var.config(text=str(round(self.current_power_experiment,3)))
            self.maximum_pow_exp_var.config(text=str(round(self.max_power_experiment,3)))
        
        if self.running_pulse:
           # self.ax_pulse_plot.set_xdata(self.time_vector_seconds)
           # self.ax_pulse_plot.set_ydata(self.power_pulse_vector)
           # self.ax_pulse.set_ylim([min(self.power_pulse_vector)*0.9, max(self.power_pulse_vector)*1.1])

            self.current_pow_pul_var.config(text=str(round(self.current_power_pulse,3)))
            self.maximum_pow_pul_var.config(text=str(round(self.max_power_pulse,3)))
        
        #self.canvas.draw()
        
        #self.ax_experiment.set_xlim([0, max(self.time_vector_seconds)])
        if  any([self.running_experiment, self.running_pulse]):
            self.gui_window.update_idletasks()
            if self.force_gui_update:
                self.gui_window.after(self.refresh_rate,self.update_gui)
            
    
    def start_gui(self):
        self.gui_window.mainloop()
        
if __name__ == '__main__':
    import numpy as np
    from matplotlib import pyplot as plt
    from scipy.io import savemat
    import pyaceqd.pulsegenerator as pg 
    import os as os
    from Pulse_v2 import fake_spectrometer, fake_motor, pulse_shaper_obj, simulator, attenuator, half_wave_plate, fake_attenuator, motor, spectrometer, load_pulse_device, create_experiment, save_pulse, save_temp, save_device, excecute_folder, power_meter, fake_power_meter 
    
    
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
    
    lab_power_meter = fake_power_meter()
    
    power_meter_A = power_meter(device=lab_power_meter,name='fake')
    
    
    
    t_0 = 0 
    t_end = 30
    chirp = 0 #45
    
    initial_pulse = pg.PulseGenerator(t_0,t_end,0.02,calibration_file=qd_calibration)
    initial_pulse.add_gaussian_time(unit = 'nm', central_f= 779.89, width_t=0.1,t0=t_end/2,area_time=10, polarisation=[1,1])
    
    lab_motor = motor(fake_motor(),name='motor1')
    lab_att = fake_attenuator()
    
    pulse_shaper = pulse_shaper_obj(device=lab_motor, calibration_file=ps_calibration, name='A: Pulse shaper')
    att_object = attenuator(lab_att, name = 'B: Attenuator')
    
    
    ps_controller = pulse_shaper.open_control(pulse_object=initial_pulse)
    att_controller = att_object.open_control(previous_control=ps_controller)
    
    
    
    pm_controller_A = power_meter_A.open_control(previous_control=att_controller, open_gui=False)
    # pm_controller_A.update_measurement()
    # att_controller.set_control_value(0.1)
    # sleep(1)
    # pm_controller_A.update_measurement()
    #pm_controller_A.toggle_running_experiment()
    #pm_controller_A.toggle_running_pulse()
    pm_controller_A.set_refresh_rate(300)
    pm_controller_A.set_record_time(60)
    
    pm_controller_A.gui()
    
    pm_controller_A.start_gui()
    