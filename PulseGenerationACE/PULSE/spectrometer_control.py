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

class spectrometer_control:
    def __init__(self,spectrometer_object,pulse_object = None, simulation_control = None, open_gui = True,  parent_window = None, display_experiment = True, display_pulse = False, display_simulation = False, autoscale = True, previous_control = None, simulation_counts = 100)->None:
        
        self.spectrometer_object = spectrometer_object
        self.name = spectrometer_object.name
        self.pulse_object = pulse_object
        self.simulation_control = simulation_control
        self.simulation_object = None
        self.open_gui = open_gui
        self.display_experiment = display_experiment
        self.display_pulse = display_pulse
        self.display_simulation = display_simulation
        self.autoscale = autoscale
        self.previous_control = previous_control
        self.simulation_counts = simulation_counts
        self.running = False
        self.view = False
        self.pulse_scale = 1
        
        self.simulation_background = 0 
        self.simulation_gaussian_noise = 0
        self.poissonian_noise = False
        
        self.simulation_mode_combined = False
        
        self.date_str_create = datetime.now().strftime("%Y_%m_%d_%H_%M")
        self.save_folder = 'spectra_'+self.date_str_create
        
        if pulse_object is not None:
            if type(pulse_object) is str:
                 pulse_object = pg.load_pulse(pulse_object)
            out_pulse_object = pulse_object.copy_pulse()
            out_pulse_object.clear_filter()
            self.pulse_object = out_pulse_object
        else:
            self.pulse_object = None
        pass
    
        if spectrometer_object.device is not None:
            self.excecute = True
        else:
            self.excecute = False
        
        self.force_gui_update = True
        
        self.print_info()
        self.update_previous_control()
        self.set_measurement_method(method='single line')
        self.set_measurement_arguments()
        if self.open_gui:
            print('hier könnte ihre GUI stehen!')
            self.gui(parentWindow=parent_window)
        
        
  
    def toggle_autoscale(self):
        self.autoscale = not self.autoscale
        pass
    
    def toggle_display_experiment(self):
        self.display_experiment = not self.display_experiment
        pass
    
    def toggle_display_pulse(self):
        self.display_pulse = not self.display_pulse
        pass
    
    def toggle_display_simulation(self):
        self.display_simulation = not self.display_simulation
        pass
    
    def toggle_force_gui_update(self):
        self.force_gui_update = not self.force_gui_update
        pass
    
    def change_view(self):
        # change through experiment, 
        self.toggle_display_experiment()
        self.toggle_display_pulse()
        self.toggle_display_simulation()
    
    def toggle_excecute(self):
        self.excecute = not self.excecute
        pass
    
    def toggle_running(self):
        self.running = not self.running
        pass
    
    def toggle_simulation_mode_combined(self):
        self.simulation_mode_combined = not self.simulation_mode_combined
    
    def set_save_folder(self,folder):
        self.save_folder = folder
    
    def set_simulation_counts(self, simulation_counts_max):
        self.simulation_counts = simulation_counts_max     
    
    def set_pulse_scale(self, pulse_scale=1):
        self.pulse_scale = pulse_scale
        
    def set_simulation_background(self, background=0):
        self.simulation_background = background
        
    def set_simulation_gaussian_noise(self, gaussian_noise=0):
        self.simulation_gaussian_noise = gaussian_noise
        
    def set_poissonian_noise(self, poissonian_noise=False):
        self.poissonian_noise = poissonian_noise
        
    def set_simulation_mode_combined(self,simulation_mode_combined):
        self.simulation_mode_combined = simulation_mode_combined
    
    def get_central_wavelength(self):
        return self.spectrometer_object.get_central_wavelength()
    
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
                self.previous_control[j].update_control()
                pulse_list.append(self.previous_control[j].get_pulse_object())
                #print('shaper: '+ self.previous_control[j].pulse_shaper_object.name)
                #print(self.previous_control[j].get_pulse_object().pulse_power)
            if pulse_list[0] is not None:
                self.pulse_object = pulse_list[0].copy_pulse()
            else:
                self.pulse_object = None
        
            for i in range(len(pulse_list)-1): 
                #print(i+1)
                self.pulse_object.merge_pulses(pulse_list[i+1])
                
        elif self.previous_control is not None:
            if self.open_gui:
                #self.previous_control.update_gui()
                pass
            self.previous_control.update_previous_control()
            self.previous_control.update_control()
            self.pulse_object = self.previous_control.get_pulse_object()
            
        if self.simulation_control is not None:
            self.simulation_object = []
            if type(self.simulation_control) is not list:
                self.simulation_control = [self.simulation_control]
            for j in range(len(self.simulation_control)):
                self.simulation_control[j].update_previous_control()
                if self.display_simulation:
                    if self.simulation_control[j].open_gui:
                        self.simulation_control[j].running = True
                        self.simulation_control[j].force_gui_update = False
                        self.simulation_control[j].display_simulation = True
                        self.simulation_control[j].update_gui()
                    else:
                        self.simulation_control[j].update_simulation()
                self.simulation_object.append(self.simulation_control[j].get_simulator_object())
    
    def get_wavelength_experiment(self):
        return self.wavelength_experiment
    
    def get_intensity_experiment(self):
        return self.intensity_experiment
    
    def get_wavelength_pulse(self):
        return self.wavelength_pulse
    
    def get_intensity_pulse(self):
        return self.intensity_pulse
    
    def get_wavelength_simulation(self):
        return self.wavelength_simulation
    
    def get_intensity_simulation(self):
        return self.intensity_simulation
    
    def get_combined_pulse_simulation(self):
        return self.combined_pulse_simulation
    
    
    def update_measurement(self):
        self.update_previous_control()
        self.spectrometer_object.get_spectrum(pulse_object=self.pulse_object, sim_object=self.simulation_object,excecute = self.excecute)
        
        self.wavelength_experiment = self.spectrometer_object.output[0][0]
        self.intensity_experiment = self.spectrometer_object.output[0][1]
        
        if self.spectrometer_object.output[1][1] is not None:
            self.wavelength_pulse = self.spectrometer_object.output[1][0]
            self.intensity_pulse = self.spectrometer_object.output[1][1]*self.pulse_scale
        else:
            self.wavelength_pulse = self.wavelength_experiment
            self.intensity_pulse = np.zeros_like(self.wavelength_experiment)
        
        if self.spectrometer_object.output[2][1] is not None:
            self.wavelength_simulation = self.spectrometer_object.output[2][0]
            self.intensity_simulation = self.spectrometer_object.output[2][1]*self.simulation_counts 
        else:
            self.wavelength_simulation = self.wavelength_experiment
            self.intensity_simulation = np.zeros_like(self.wavelength_experiment)
        
        
        
        if self.poissonian_noise:
            self.intensity_simulation = np.random.poisson(self.intensity_simulation).astype(float)
            self.intensity_pulse = np.random.poisson(self.intensity_pulse).astype(float)
        
        self.combined_pulse_simulation = self.intensity_pulse + self.intensity_simulation
        
        self.intensity_pulse += self.simulation_background*np.ones_like(self.intensity_simulation)
        self.intensity_simulation += self.simulation_background*np.ones_like(self.intensity_simulation)
        self.combined_pulse_simulation += self.simulation_background*np.ones_like(self.intensity_simulation)
        
        if self.simulation_gaussian_noise != 0:
            noise = np.random.normal(0,self.simulation_gaussian_noise,len(self.intensity_pulse))
            self.intensity_pulse += noise
            self.intensity_simulation += noise
            self.combined_pulse_simulation += noise
            
        
        
    def single_line_measurement(self,arguments =[]):
        # measurment_kind = 0: experiment, 1: pulse, 2: simulation, 3: combined 
        center_wl =  float(arguments[0])
        width_wl = float(arguments[1])
        method = str(arguments[2])
        
        if len(arguments) == 4 and arguments[3] != '':
            target = float(arguments[3])
        else:
            target = None
        
        self.update_measurement()
        if self.display_experiment:
            wavelength = self.wavelength_experiment
            intensity = self.intensity_experiment
        elif self.display_pulse and self.display_simulation:
            wavelength = self.wavelength_pulse
            intensity = self.combined_pulse_simulation
        elif self.display_pulse:
            wavelength = self.wavelength_pulse
            intensity = self.intensity_pulse
        elif self.display_simulation:
            wavelength = self.wavelength_simulation
            intensity = self.intensity_simulation
        
        if width_wl == 0:
            wavelength_window = np.argmin(np.abs(wavelength-center_wl))
        else:
            wavelength_window = np.where(np.abs(wavelength-center_wl)<width_wl/2)
        
        if method == 'max':
            measurement = np.max(intensity[wavelength_window])
        elif method == 'mean':
            measurement = np.mean(intensity[wavelength_window])
        elif method == 'sum':
            measurement = np.sum(intensity[wavelength_window])
        
        if target is not None:
            measurement = target - np.abs(measurement - target)
        
        self.measurment = measurement
        print('Measurement: '+str(measurement))
        return self.measurment
    
    def get_measurement_result(self):
        return self.measurment
        
    def single_line_measurement_gui(self, control = None, row_offset = 0, column_offset = 0):
    
        control.measurement_args_spec = [None,None, None]
        num_rows = 5
        
        tk.Label(control.gui_window, text= self.name + ': Single line measurement').grid(row =row_offset +1, column=0+column_offset)
        tk.Label(control.gui_window, text= 'Center wavelength: ').grid(row = row_offset+2,column=0+column_offset)
        
        control.measurement_args_spec[0] = tk.Entry(control.gui_window)
        control.measurement_args_spec[0].grid(row=row_offset+2,column=1+column_offset)
        control.measurement_args_spec[0].insert(0,str(self.measurement_arguments[0])) 
        
        tk.Label(control.gui_window, text='Width wavelength: ').grid(row=row_offset+3,column=0+column_offset)
        control.measurement_args_spec[1] = tk.Entry(control.gui_window)
        control.measurement_args_spec[1].grid(row=row_offset+3,column=1+column_offset)
        control.measurement_args_spec[1].insert(0,str(self.measurement_arguments[1]))
        
        tk.Label(control.gui_window, text='Method: ').grid(row=row_offset+4,column=0+column_offset)
        control.measurement_args_spec[2] = tk.Entry(control.gui_window)
        control.measurement_args_spec[2].grid(row=row_offset+4,column=1+column_offset)
        control.measurement_args_spec[2].insert(0,self.measurement_arguments[2])
        
        tk.Label(control.gui_window, text='Target: ').grid(row=row_offset+5,column=0+column_offset)
        control.measurement_args_spec.append(tk.Entry(control.gui_window))
        control.measurement_args_spec[-1].grid(row=row_offset+5,column=1+column_offset)
        if len(self.measurement_arguments) == 4:
            control.measurement_args_spec[-1].insert(0,self.measurement_arguments[3])
        
        return num_rows
    
    def multi_line_measurement(self,arguments = []):
        output_list = []
        for i in range(int(arguments[0])):
            center_wl =  float(arguments[1+i*4])
            width_wl = float(arguments[2+i*4])
            method = str(arguments[3+i*4])
            target = str(arguments[4+i*4])
            output_list.append(self.single_line_measurement([center_wl,width_wl,method,target]))
        
        self.measurment = output_list
        return self.measurment
    
    def multi_line_measurement_gui(self, control = None, row_offset = 0, column_offset = 0):
        control.measurement_args_spec = [None]
        num_rows = 5
        
        tk.Label(control.gui_window, text= self.name + ': Multi line measurement - ' + str(self.measurement_arguments[0])+' lines').grid(row =row_offset +1, column=0+column_offset)
        
        control.measurement_args_spec[0] = tk.Entry(control.gui_window)
        # control.measurement_args_spec[0].grid(row=row_offset+3,column=0+column_offset)
        control.measurement_args_spec[0].insert(0,str(self.measurement_arguments[0])) 
        
        tk.Label(control.gui_window, text= 'Center wavelength: ').grid(row = row_offset+2,column=column_offset)
        tk.Label(control.gui_window, text='Width wavelength: ').grid(row=row_offset+3,column=column_offset)
        tk.Label(control.gui_window, text='Method: ').grid(row=row_offset+4,column=column_offset)
        tk.Label(control.gui_window, text='Target: ').grid(row=row_offset+5,column=column_offset)
        
        for i in range(int(self.measurement_arguments[0])):
            tk.Label(control.gui_window, text= 'Line '+str(i+1)).grid(row =row_offset +1, column=1+column_offset+i)
            control.measurement_args_spec.append(tk.Entry(control.gui_window))
            control.measurement_args_spec[-1].grid(row=row_offset+2,column=1+column_offset+i)
            control.measurement_args_spec[-1].insert(0,str(self.measurement_arguments[1+i*4]))
            
            control.measurement_args_spec.append(tk.Entry(control.gui_window))
            control.measurement_args_spec[-1].grid(row=row_offset+3,column=1+column_offset+i)
            control.measurement_args_spec[-1].insert(0,str(self.measurement_arguments[2+i*4]))
            
            control.measurement_args_spec.append(tk.Entry(control.gui_window))
            control.measurement_args_spec[-1].grid(row=row_offset+4,column=1+column_offset+i)
            control.measurement_args_spec[-1].insert(0,str(self.measurement_arguments[3+i*4]))
            
            control.measurement_args_spec.append(tk.Entry(control.gui_window))
            control.measurement_args_spec[-1].grid(row=row_offset+5,column=1+column_offset+i)
            control.measurement_args_spec[-1].insert(0,str(self.measurement_arguments[4+i*4]))
            
        return num_rows 
        
    
    def set_measurement_arguments(self, arguments = None):
        if self.measurement_method == self.single_line_measurement and arguments is None:
            self.measurement_arguments = [str(self.get_central_wavelength()),'0.2','max']
        elif self.measurement_method == self.multi_line_measurement and arguments is None:
            for i in range(self.measurement_arguments[0]):
                self.measurement_arguments.append(str(self.get_central_wavelength()))
                self.measurement_arguments.append('0.2')
                self.measurement_arguments.append('max')
        else:
            for i in range(len(arguments)):
                self.measurement_arguments[i] = arguments[i]
    
    def set_measurement_method(self, method = 'single line', arg = 2):
        if method == 'single line':
            self.measurement_method = self.single_line_measurement
            self.measurement_method_gui = self.single_line_measurement_gui
            self.number_measurement_outputs = 1
        elif method == 'multi line':
            self.measurement_method = self.multi_line_measurement
            self.measurement_method_gui = self.multi_line_measurement_gui
            self.measurement_arguments = [arg]
            self.number_measurement_outputs = arg
            for i in range(self.measurement_arguments[0]):
                self.measurement_arguments.append(str(self.get_central_wavelength()))
                self.measurement_arguments.append('0.2')
                self.measurement_arguments.append('max')
                self.measurement_arguments.append('')

    def get_number_measurement_outputs(self):
        return self.number_measurement_outputs
        
    def get_measurement_args(self,control):
        args = []
        for i in range(len(control.measurement_args_spec)):
            args.append(control.measurement_args_spec[i].get())
        return args
    
    def print_info(self):
        print('Spectrometer Control')
        print('Spectrometer Object: '+self.spectrometer_object.name)
        #print('Pulse Object: '+self.pulse_object.name)
        #print('Simulation Object: '+self.simulation_control.name)
        print('Open GUI: '+str(self.open_gui))
        print('Display Experiment: '+str(self.display_experiment))
        print('Display Pulse: '+str(self.display_pulse))
        print('Display Simulation: '+str(self.display_simulation))
        print('Autoscale: '+str(self.autoscale))
        pass
    
    def close(self):
        self.spectrometer_object.close()
    
    def button_color(self, button_object, state):
            if state:
                button_object.config(bg='green')
            else:
                button_object.config(bg='white')
    
    def update_gui(self): 
            
            # self.update_previous_control()
            # self.spectrometer_object.get_spectrum(pulse_object=self.pulse_object, sim_object=self.simulation_control,excecute = self.excecute)
        self.update_measurement()
        self.folder_entry.delete(0,tk.END)
        self.folder_entry.insert(0,self.save_folder)
        # ------
        #self.ax_experiment.set_xlim(float(self.min_wavelength_entry.get()),float(self.max_wavelength_entry.get()))
        self.ax_experiment.set_xlim(self.x_lim_min,self.x_lim_max)
        # -----
        
        if self.display_experiment: 
            ## add colouring of buttons here 
            
            
            self.ax_experiment_plot.set_xdata(self.spectrometer_object.output[0][0])
            self.ax_experiment_plot.set_ydata(self.intensity_experiment)
        
            #self.ax_experiment.set_xlim(self.spectrometer_object.output[0][0][0],self.spectrometer_object.output[0][0][-1])
            
            if self.autoscale:
                self.y_lim_min_exp = int(np.floor(np.min(self.intensity_experiment)))
                self.y_lim_max_exp = int(np.ceil(np.max(self.intensity_experiment)))
                self.ax_experiment.set_ylim(self.y_lim_min_exp,self.y_lim_max_exp)
                
                #self.min_exp_entry.delete(0,tk.END)
                #self.min_exp_entry.insert(0,str(self.y_lim_min_exp))
                
                
                
                # self.max_exp_entry.delete(0,tk.END)
                # self.max_exp_entry.insert(0,str(self.y_lim_max_exp))
                
            else:
                self.ax_experiment.set_ylim(int(self.y_lim_min_exp),int(self.y_lim_max_exp))

            self.min_exp_text.config(text=str(int(self.y_lim_min_exp)))
            self.max_exp_text.config(text=str(int(self.y_lim_max_exp)))
            
        else:
            self.ax_experiment_plot.set_xdata([])
            self.ax_experiment_plot.set_ydata([])
        
        if self.display_pulse:
            self.ax_pulse_plot.set_xdata(self.wavelength_pulse)
            self.ax_pulse_plot.set_ydata(self.intensity_pulse)
            
            #self.ax_simulation.set_xlim(self.spectrometer_object.output[1][0][0],self.spectrometer_object.output[1][0][-1])
            
            if self.autoscale:
                if self.display_simulation:
                    self.y_lim_min = int(np.floor(np.min(np.concatenate((self.intensity_pulse,self.intensity_simulation)))))
                    self.y_lim_max = int(np.ceil(np.max(np.concatenate((self.intensity_pulse,self.intensity_simulation)))))

                else:
                    self.y_lim_min = int(np.floor(np.min(self.intensity_pulse)))
                    self.y_lim_max = int(np.ceil(np.max(self.intensity_pulse)))
                self.ax_simulation.set_ylim(np.floor(self.y_lim_min),np.ceil(self.y_lim_max))
                
                # self.min_sim_entry.delete(0,tk.END)
                # self.min_sim_entry.insert(0,str(int(np.floor(self.y_lim_min))))
                
                # self.max_sim_entry.delete(0,tk.END)
                # self.max_sim_entry.insert(0,str(int(np.ceil(self.y_lim_max))))
                
            else:
                self.ax_simulation.set_ylim(int(self.y_lim_min),int(self.y_lim_max))
            
            self.min_sim_text.config(text=str(int(np.floor(self.y_lim_min))))
            self.max_sim_text.config(text=str(int(np.ceil(self.y_lim_max))))
        else:
            self.ax_pulse_plot.set_xdata([])
            self.ax_pulse_plot.set_ydata([])
            
        if self.display_simulation:
            #self.set_simulation_counts(int(self.simulation_counts_entry.get()))
            # self.simulation_counts_entry.delete(0,tk.END)
            # self.simulation_counts_entry.insert(0,str(self.simulation_counts))
            self.ax_simulation_plot.set_xdata(self.wavelength_simulation)
            if not self.simulation_mode_combined:
                self.ax_simulation_plot.set_ydata(self.intensity_simulation)
            else:
                self.ax_simulation_plot.set_ydata(self.combined_pulse_simulation)
            
            if self.autoscale:
                if self.simulation_mode_combined:
                    self.y_lim_min = int(np.floor(np.min(self.combined_pulse_simulation)))
                    self.y_lim_max = int(np.ceil(np.max(self.combined_pulse_simulation)))
                elif self.display_pulse:
                    self.y_lim_min = int(np.floor(np.min(np.concatenate((self.intensity_pulse,self.intensity_simulation)))))
                    self.y_lim_max = int(np.ceil(np.max(np.concatenate((self.intensity_pulse,self.intensity_simulation)))))
                else:
                    self.y_lim_min = int(np.floor(np.min(self.intensity_simulation)))
                    self.y_lim_max = int(np.ceil(np.max(self.intensity_simulation)))
                self.ax_simulation.set_ylim(np.floor(self.y_lim_min),np.ceil(self.y_lim_max))
                
                # self.min_sim_entry.delete(0,tk.END)
                # self.min_sim_entry.insert(0,str(int(np.floor(self.y_lim_min))))
                
                # self.max_sim_entry.delete(0,tk.END)
                # self.max_sim_entry.insert(0,str(int(np.ceil(self.y_lim_max))))
                
            else:
                self.ax_simulation.set_ylim(self.y_lim_min,self.y_lim_max)
            # fix this to something smart!!!
            self.min_sim_text.config(text=str(int(np.floor(self.y_lim_min))))
            self.max_sim_text.config(text=str(int(np.ceil(self.y_lim_max))))
            
        else:
            self.ax_simulation_plot.set_xdata([])
            self.ax_simulation_plot.set_ydata([])
            
        if self.running:
            self.canvas_fig_spec.draw()
            
            if self.force_gui_update:
                #self.set_simulation_counts(int(self.simulation_counts_entry.get()))
                self.gui_window.after(200,self.update_gui)
            
            #self.simulation_counts_entry.delete(0,tk.END)
            #self.simulation_counts_entry.insert(0,str(self.simulation_counts))
    
    
    def gui(self,parentWindow=None):
        #self.xlim = [self.spectrometer_object.exp_wl_0[0],self.spectrometer_object.exp_wl_0[-1]]
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
        
        # if self.previous_control is None:
        #     self.gui_window = tk.Tk()
        # else:
        #     if type(self.previous_control) is list:
        #         self.gui_window = tk.Toplevel(self.previous_control[0].gui_window)
        #     else:
        #         self.gui_window = tk.Toplevel(self.previous_control.gui_window)
        def close_gui():
            self.running = False
            if self.spectrometer_object.device is not None:
                self.close()
            self.gui_window.destroy()
            print(self.name+' closed')
        
        self.gui_window.protocol("WM_DELETE_WINDOW", close_gui)
        self.gui_window.title('Spectrometer Control: '+self.spectrometer_object.name)
        #self.gui_window.geometry('800x600')
        
        fig_spectrometer = Figure(figsize=(5, 4), dpi=100)
        # twin x axis for spectrum
        self.ax_experiment = fig_spectrometer.add_subplot()
        self.ax_simulation = self.ax_experiment.twinx()
        
        self.ax_experiment_plot, = self.ax_experiment.plot([],[],'k')
        self.ax_pulse_plot, = self.ax_simulation.plot([],[],'r')
        self.ax_simulation_plot, = self.ax_simulation.plot([],[],'b')
        
        
        self.ax_experiment.set_xlabel('Wavelength [nm]')
        self.x_lim_min = self.spectrometer_object.exp_wl_0[0]
        self.x_lim_max = self.spectrometer_object.exp_wl_0[-1]
        self.ax_experiment.set_xlim(self.x_lim_min,self.x_lim_max)
        
        self.y_lim_min_exp = 0
        self.y_lim_max_exp = 1
        self.ax_experiment.set_ylim(self.y_lim_min_exp,self.y_lim_max_exp)
        
        self.y_lim_min = 0
        self.y_lim_max = 1
        self.ax_simulation.set_ylim(self.y_lim_min,self.y_lim_max)
        
        self.canvas_fig_spec = FigureCanvasTkAgg(fig_spectrometer, master=self.gui_window)
        self.canvas_fig_spec.draw()
        
        
                
                #print('debugging')
                #print('lol')
                #print(str(self.pulse_object.pulse_power))
            
        
        
        self.canvas_fig_spec.get_tk_widget().grid(row=0,column=0, columnspan=6)
        
        
        self.run_button = tk.Button(self.gui_window, text='Run', command= lambda: [self.toggle_running(),self.button_color(self.run_button,self.running),self.update_gui()],bg='white')
        self.run_button.grid(row=1,column=0)
        
        self.autoscale_button = tk.Button(self.gui_window, text='Autoscale',bg='white')
        self.autoscale_button.config(command= lambda: [self.toggle_autoscale(),self.button_color(self.autoscale_button,self.autoscale)])
        self.autoscale_button.grid(row=1,column=1)
        
        self.settings_button = tk.Button(self.gui_window, text='Settings')
        self.settings_button.config(command= self.settings_gui)   
        self.settings_button.grid(row=1,column=2)    
        
        #simulation_counts_label = tk.Label(self.gui_window, text='Simulation Counts').grid(row=1,column=2)
        
        #self.simulation_counts_entry = tk.Entry(self.gui_window, width=10)
        #self.simulation_counts_entry.insert(0,str(self.simulation_counts))
        #self.simulation_counts_entry.grid(row=1,column=3)
        
        
        
        
        self.min_exp_label = tk.Label(self.gui_window, text='Min Exp.').grid(row=2,column=0)
        
        # self.min_exp_entry = tk.Entry(self.gui_window, width=10)
        # self.min_exp_entry.insert(0,str(0))
        # self.min_exp_entry.grid(row=2,column=1)
        
        self.min_exp_text = tk.Label(self.gui_window, text=str(self.y_lim_min_exp))
        self.min_exp_text.grid(row=2,column=1)
        
        self.max_exp_label = tk.Label(self.gui_window, text='Max Exp.').grid(row=2,column=2)
        
        # self.max_exp_entry = tk.Entry(self.gui_window, width=10)
        # self.max_exp_entry.insert(0,str(1))
        # self.max_exp_entry.grid(row=2,column=3)
        
        self.max_exp_text = tk.Label(self.gui_window, text=str(self.y_lim_max_exp))
        self.max_exp_text.grid(row=2,column=3)
        
        self.min_sim_label = tk.Label(self.gui_window, text='Min Sim.').grid(row=3,column=0)
        
        # self.min_sim_entry = tk.Entry(self.gui_window, width=10)
        # self.min_sim_entry.insert(0,str(0))
        # self.min_sim_entry.grid(row=3,column=1)
        
        self.min_sim_text = tk.Label(self.gui_window, text=str(self.y_lim_min))
        self.min_sim_text.grid(row=3,column=1)
        
        self.max_sim_label = tk.Label(self.gui_window, text='Max Sim.').grid(row=3,column=2)
        
        # self.max_sim_entry = tk.Entry(self.gui_window, width=10)
        # self.max_sim_entry.insert(0,str(1))
        # self.max_sim_entry.grid(row=3,column=3)
        
        self.max_sim_text = tk.Label(self.gui_window, text=str(self.y_lim_max))
        self.max_sim_text.grid(row=3,column=3)
        
        
        self.toggle_display_experiment_button = tk.Button(self.gui_window, text='Display Experiment')
        self.toggle_display_experiment_button.config(command= lambda: [self.toggle_display_experiment(), self.button_color(self.toggle_display_experiment_button,self.display_experiment)])
        
        self.toggle_display_experiment_button.grid(row=4,column=0)
        
        if self.previous_control[0].get_pulse_object() is not None:
            button_state = tk.NORMAL
        else:
            button_state = tk.DISABLED
            self.display_pulse = False
        
        self.toggle_display_pulse_button = tk.Button(self.gui_window, text='Display Pulse',state=button_state)
        self.toggle_display_pulse_button.config(command= lambda: [self.toggle_display_pulse(),self.button_color(self.toggle_display_pulse_button,self.display_pulse)])
        self.toggle_display_pulse_button.grid(row=4,column=1)
        
        if self.simulation_control is not None:
            button_state = tk.NORMAL
        else:
            button_state = tk.DISABLED
            self.display_simulation = False
        
        self.toggle_display_simulation_button = tk.Button(self.gui_window, text='Display Simulation', state=button_state)
        self.toggle_display_simulation_button.config(command= lambda: [self.toggle_display_simulation(),self.button_color(self.toggle_display_simulation_button,self.display_simulation)])
        self.toggle_display_simulation_button.grid(row=4,column=2)
        
        
        # color buttons 
        self.button_color(self.toggle_display_experiment_button,self.display_experiment)
        self.button_color(self.toggle_display_pulse_button,self.display_pulse)
        self.button_color(self.toggle_display_simulation_button,self.display_simulation)
        self.button_color(self.autoscale_button,self.autoscale)
        self.button_color(self.run_button,self.running)
        
        # self.min_wavelength_label = tk.Label(self.gui_window, text='Min Wavelength').grid(row=5,column=0)
        
        # self.min_wavelength_entry = tk.Entry(self.gui_window, width=10)
        # self.min_wavelength_entry.insert(0,str(self.x_lim_min))
        # self.min_wavelength_entry.grid(row=5,column=1)
        
        # self.max_wavelength_label = tk.Label(self.gui_window, text='Max Wavelength').grid(row=5,column=2)
        
        # self.max_wavelength_entry = tk.Entry(self.gui_window, width=10)
        # self.max_wavelength_entry.insert(0,str(self.x_lim_max))
        # self.max_wavelength_entry.grid(row=5,column=3)
        
        self.save_spectrum_button = tk.Button(self.gui_window, text='Save Spectrum')
        self.save_spectrum_button.config(command= lambda: self.save_spectrum_gui())
        self.save_spectrum_button.grid(row=6,column=0) 
        
        tk.Label(self.gui_window, text='Folder: ').grid(row=6,column=1)
        
        self.folder_entry = tk.Entry(self.gui_window, width=10)
        self.folder_entry.insert(0,self.save_folder)
        self.folder_entry.grid(row=6,column=2)
        
        self.change_folder_button = tk.Button(self.gui_window, text='Change Folder')
        self.change_folder_button.config(command= lambda: [self.change_folder_gui(),self.update_gui()])
        self.change_folder_button.grid(row=6,column=3)
        
        if self.running:
            self.update_gui()
    
    def settings_gui(self):
        self.settings_window = tk.Toplevel(self.gui_window)
        self.settings_window.title('Spectrometer settings')
        self.set_settings_button = tk.Button(self.settings_window, text='Set Settings')
        self.set_settings_button.config(command= lambda: [self.set_settings(float(self.min_wavelength_entry.get()),float(self.max_wavelength_entry.get()),float(self.simulation_background_entry.get()),float(self.simulation_gaussian_noise_entry.get()),self.poissonian_noise, float(self.simulation_counts_entry.get()), float(self.pulse_scale_entry.get()) ,float(self.min_exp_entry.get()), float(self.max_exp_entry.get()), float(self.min_sim_entry.get()), float(self.max_sim_entry.get()), self.lp_angle_entry.get())])
        self.set_settings_button.grid(row=0,column=0)
        
        self.noise_off_button = tk.Button(self.settings_window, text='Noise Off')
        self.noise_off_button.config(command= lambda: [self.set_simulation_background(0),self.set_simulation_gaussian_noise(0),self.set_poissonian_noise(False),self.button_color(self.poissonian_noise_button,self.poissonian_noise)])
        self.noise_off_button.grid(row=0,column=1)
        
        self.min_wavelength_label = tk.Label(self.settings_window, text='Min Wavelength').grid(row=1,column=0)
        
        self.min_wavelength_entry = tk.Entry(self.settings_window, width=10)
        self.min_wavelength_entry.insert(0,str(self.x_lim_min))
        self.min_wavelength_entry.grid(row=1,column=1)
        
        self.max_wavelength_label = tk.Label(self.settings_window, text='Max Wavelength').grid(row=2,column=0)
        
        self.max_wavelength_entry = tk.Entry(self.settings_window, width=10)
        self.max_wavelength_entry.insert(0,str(self.x_lim_max))
        self.max_wavelength_entry.grid(row=2,column=1)
        
        self.simulation_background_label = tk.Label(self.settings_window, text='Simulation Background').grid(row=3,column=0)
        
        self.simulation_background_entry = tk.Entry(self.settings_window, width=10)
        self.simulation_background_entry.insert(0,str(self.simulation_background))
        self.simulation_background_entry.grid(row=3,column=1)
        
        self.simulation_gaussian_noise_label = tk.Label(self.settings_window, text='Simulation Gaussian Noise').grid(row=4,column=0)
        
        self.simulation_gaussian_noise_entry = tk.Entry(self.settings_window, width=10)
        self.simulation_gaussian_noise_entry.insert(0,str(self.simulation_gaussian_noise))
        self.simulation_gaussian_noise_entry.grid(row=4,column=1)
        
        self.poissonian_noise_label = tk.Label(self.settings_window, text='Poissonian Noise').grid(row=5,column=0)
        
        
        self.poissonian_noise_button = tk.Button(self.settings_window, text='Toggle')
        self.poissonian_noise_button.config(command= lambda: [self.set_poissonian_noise(not self.poissonian_noise),self.button_color(self.poissonian_noise_button,self.poissonian_noise)])
        self.poissonian_noise_button.grid(row=5,column=1)
        self.button_color(self.poissonian_noise_button,self.poissonian_noise)
        
        self.simulation_counts_label = tk.Label(self.settings_window, text='Simulation Counts').grid(row=6,column=0)
        
        self.simulation_counts_entry = tk.Entry(self.settings_window, width=10)
        self.simulation_counts_entry.insert(0,str(self.simulation_counts))
        self.simulation_counts_entry.grid(row=6,column=1)
        
        self.pulse_scale_label = tk.Label(self.settings_window, text='Pulse Scale').grid(row=7,column=0)
        
        self.pulse_scale_entry = tk.Entry(self.settings_window, width=10)
        self.pulse_scale_entry.insert(0,str(self.pulse_scale))
        self.pulse_scale_entry.grid(row=7,column=1)
        
        self.simulation_mode_combined_label = tk.Label(self.settings_window,text='Combined view').grid(row=8,column=0)
        
        self.simulation_mode_combined_button = tk.Button(self.settings_window, text='Toggle')
        self.simulation_mode_combined_button.config(command= lambda: [self.toggle_simulation_mode_combined(), self.button_color(self.simulation_mode_combined_button,self.simulation_mode_combined)])
        self.simulation_mode_combined_button.grid(row=8, column=1)
        self.button_color(self.simulation_mode_combined_button,self.simulation_mode_combined)
        
        tk.Label(self.settings_window, text='LP angle: ').grid(row=9,column=0)
        
        self.lp_angle_entry = tk.Entry(self.settings_window, width=10)
        if self.spectrometer_object.polarisation_angle is None:
            self.lp_angle_entry.insert(0,'')
        else:
            self.lp_angle_entry.insert(0,str(self.spectrometer_object.polarisation_angle))
        self.lp_angle_entry.grid(row=9,column=1)
        
        tk.Label(self.settings_window, text='-----').grid(row=10,column=0)
        
        self.min_exp_label = tk.Label(self.settings_window, text='Min Exp.').grid(row=11,column=0)
        
        self.min_exp_entry = tk.Entry(self.settings_window, width=10)
        self.min_exp_entry.insert(0,str(self.y_lim_min_exp))
        self.min_exp_entry.grid(row=11,column=1)
        
        self.max_exp_label = tk.Label(self.settings_window, text='Max Exp.').grid(row=11,column=2)
        
        self.max_exp_entry = tk.Entry(self.settings_window, width=10)
        self.max_exp_entry.insert(0,str(self.y_lim_max_exp))
        self.max_exp_entry.grid(row=11,column=3)
        
        self.min_sim_label = tk.Label(self.settings_window, text='Min Sim.').grid(row=12,column=0)
        
        self.min_sim_entry = tk.Entry(self.settings_window, width=10)
        self.min_sim_entry.insert(0,str(self.y_lim_min))
        self.min_sim_entry.grid(row=12,column=1)
        
        self.max_sim_label = tk.Label(self.settings_window, text='Max Sim.').grid(row=12,column=2)
        
        self.max_sim_entry = tk.Entry(self.settings_window, width=10)
        self.max_sim_entry.insert(0,str(self.y_lim_max))
        self.max_sim_entry.grid(row=12,column=3)
        
        
        
        
        pass
    
    def set_settings(self, min_wl, max_wl, background = 0, gaussian_noise = 0, poissonian_noise = False,simulation_counts = 1, pulse_scale = 1,y_lim_min_exp = 0, y_lim_max_exp = 1, y_lim_min = 0, y_lim_max = 1, polarisation_angle = None):
        self.x_lim_min = min_wl
        self.x_lim_max = max_wl
        self.simulation_background = background
        self.simulation_gaussian_noise = gaussian_noise
        self.poissonian_noise = poissonian_noise
        self.simulation_counts = simulation_counts
        self.pulse_scale = pulse_scale
        self.y_lim_min_exp = y_lim_min_exp
        self.y_lim_max_exp = y_lim_max_exp
        self.y_lim_min = y_lim_min
        self.y_lim_max = y_lim_max
        if polarisation_angle == '':
            polarisation_angle = None
        else:
            polarisation_angle = float(polarisation_angle)
        self.spectrometer_object.set_polarisation(polarisation_angle)
        
        pass
        
        
    
    def save_spectrum_gui(self):
        self.set_save_folder(self.folder_entry.get())
        self.save_spectrum()
    
    def save_spectrum(self):
        now = datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
        save_name = self.save_folder+'/'+'spectrum_'+now+'.txt'
        
        if not os.path.exists(self.save_folder):
            os.makedirs(self.save_folder)
            
        with open(save_name, 'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(['Wavelength (nm)','Experiment','Pulse','Simulation'])
            for i in range(len(self.wavelength_experiment)):
                writer.writerow([self.wavelength_experiment[i],self.intensity_experiment[i],self.intensity_pulse[i],self.intensity_simulation[i]])
        file.close()
    
    
    
    def change_folder_gui(self):
        self.save_folder = tk.filedialog.askdirectory()
        pass
    
    def start_gui(self):
        self.gui_window.mainloop()
        pass
    
if __name__ == '__main__':
    import numpy as np
    from matplotlib import pyplot as plt
    from scipy.io import savemat
    import pyaceqd.pulsegenerator as pg 
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
    lab_motor3 = motor(fake_motor(),name='motor3')
    
    lab_att= fake_attenuator()
    att_object = attenuator(lab_att, name = 'F Att')
    att_controller = att_object.open_control(pulse_object=initial_pulse,previous_control=None)
    
    pulse_shaper = pulse_shaper_obj(device=lab_motor, calibration_file=ps_calibration, name='heisl')
    ps_controller = pulse_shaper.open_control(pulse_object=None,previous_control = att_controller)
    
    pulse_shaper2 = pulse_shaper_obj(device=lab_motor2, calibration_file=ps_calibration, name='bert')
    ps_controller2 = pulse_shaper2.open_control(pulse_object = initial_pulse,previous_control = ps_controller)
    
    pulse_shaper3 = pulse_shaper_obj(device=lab_motor3, calibration_file=ps_calibration, name='gregor')
    ps_controller3 = pulse_shaper3.open_control(pulse_object = initial_pulse,previous_control = None) 
    
    lab_spectromter = fake_spectrometer(start_wl=774, end_wl = 790, n_wl = 1024)

    spec = spectrometer(device=lab_spectromter, name='MF_spec') 
    sm_controller = spec.open_control(pulse_object=None,previous_control =[ps_controller2, ps_controller3]) #, ps_controller3#ps_controller2.get_pulse_object()
    #sm_controller.toggle_display_pulse()
    sm_controller.print_info()
    ps_controller.start_gui()
        


    