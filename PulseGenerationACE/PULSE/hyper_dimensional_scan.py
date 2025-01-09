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

import pandas as pd


class hyper_scan():
    def __init__(self, device_control: list = [], measururement_control: list = [], open_gui = True, name = 'Hyper_scan') -> None: 
        
        if type(device_control) != list:
            device_control = [device_control]
            
        if type(measururement_control) != list:
            measururement_control = [measururement_control]
        
        self.name = name
        self.device_control_initial = device_control
        self.measururement_control_initial = measururement_control
        
        self.scan_limits = []
        self.scan_num = []
        for i in range(len(self.device_control_initial)):
            cur_control_value = self.device_control_initial[i].get_control_value()
            self.scan_limits.append([cur_control_value-0.5, cur_control_value+0.5])
            
            self.scan_num.append(10)
        
        self.action_multiplier = 1
        self.device_control = []
        self.scan_vectors = []
        self.measururement_control = []
        self.device_string = 'Devices: \n'
        self.measurement_string = 'Measurements: \n'
        self.open_gui = open_gui
        self.save_name = 'scan_name.txt'
        self.snapshot_name = 'snapshot_name.txt'
        self.date_str_create = datetime.now().strftime("%Y_%m_%d_%H_%M")
        self.snapshot_folder = 'snapshots_'+self.date_str_create
        self.hyper_scan_folder = 'hyper_scans_'+self.date_str_create
        if self.open_gui:
            self.gui()
        
    
    def set_save_name(self, save_name):
        self.save_name = save_name
    
    def set_snapshot_folder(self, snapshot_folder = ''):
        self.snapshot_folder = snapshot_folder
        
    def set_hyper_scan_folder(self, hyper_scan_folder):
        self.hyper_scan_folder = hyper_scan_folder
    
    def reset_scan(self):
        self.device_control = []
        self.measururement_control = []
        self.scan_vectors = []
        self.device_string = 'Devices: \n'
        self.measurement_string = 'Measurements: \n'
        self.action_multiplier = 1
        self.save_name = 'scan_name.txt'
        
        
    def set_scan_limits(self, device_index = 0, scan_limits = [0,1]):
        self.scan_limits[device_index] = scan_limits 
        
    def set_scan_num(self, device_index = 0, scan_num = 10):
        self.scan_num[device_index] = scan_num
    
    def get_scan_limits(self, device_index = 0):
        return self.scan_limits[device_index]    
    
    def create_scan_vector(self,device_index = 0):
        return np.linspace(self.scan_limits[device_index][0],self.scan_limits[device_index][1],self.scan_num[device_index])
    
    
    def add_device_to_scan(self, device_index = 0):
        self.scan_vectors.append(self.create_scan_vector(device_index))
        self.device_control.append(self.device_control_initial[device_index])
        self.action_multiplier *= self.scan_num[device_index]
        self.device_string += self.device_control_initial[device_index].name +': '+ str(self.action_multiplier) +' actions \n' 
        
    def add_measurement_to_scan(self, measurement_index = 0):
        self.measururement_control.append(self.measururement_control_initial[measurement_index])
        self.measurement_string += self.measururement_control_initial[measurement_index].name + '\n'
        
    def compile_scan(self):
        self.total_scan_num = 1
        self.total_scan = []
        for i in range(len(self.scan_vectors)):
            self.total_scan_num = self.total_scan_num*len(self.scan_vectors[i])
        
        self.device_string += 'Total moves: ' + str(self.total_scan_num)
        self.measurement_string += 'Total measurements: ' + str(len(self.measururement_control)*self.total_scan_num)
        
        now = datetime.now()
        now_str = now.strftime("%Y_%m_%d_%H_%M")
        self.save_name = 'Hyper_scan_'+str(len(self.device_control))+'D_'+now_str+'.txt'
        
        self.save_header = 'Hyper scan: '+now_str+'\nNumber of devices: '+str(len(self.device_control_initial))+ '\nNumber of active devices: '+str(len(self.device_control)) + '\nNumber of measurement devices: '+str(len(self.measururement_control)) + '\nTotal moves: '+str(self.total_scan_num)+ '\nActive devices: '+str(self.device_string) + '\nActive measurements: '+str(self.measurement_string) + '\n'
        
        device_name_str = 'Number\t'
        for i in range(len(self.device_control_initial)):
            device_name_str += self.device_control_initial[i].name + '\t'
            
        for i in range(len(self.measururement_control)):
            for j in range(self.measururement_control[i].get_number_measurement_outputs()): #### changed
                device_name_str += self.measururement_control[i].name +'_'+str(j)+'\t'
        self.save_header += device_name_str + '\n'
        
        def compile_helper(device_index = 0, scan_index = 0):
            if device_index == len(self.scan_vectors):
                return 
            for i in range(len(self.scan_vectors[device_index])):
                self.total_scan.append([0,self.device_control[device_index],self.scan_vectors[device_index][i]])
                if device_index == len(self.scan_vectors)-1:
                    for measurment_index in range(len(self.measururement_control)):
                        self.total_scan.append([1,self.measururement_control[measurment_index],None])
                        if measurment_index == len(self.measururement_control)-1:
                            self.total_scan.append([2,None,None])
                    pass
                
                compile_helper(device_index+1)
        
        compile_helper()
    
    
    def get_total_scan_num(self):
        return self.total_scan_num
    
    def get_total_scan(self):
        return self.total_scan
    
    def get_device_state(self):
        device_state = []
        for device in self.device_control_initial:
            device_state.append(device.get_control_value())
        for measurement in self.measururement_control:
            if type(measurement.get_measurement_result()) == list:
                for result in measurement.get_measurement_result():
                    device_state.append(result)
            else:
                device_state.append(measurement.get_measurement_result())
        return device_state
    
    
    def snapshot(self):
        for measurement in self.measururement_control:
            args = measurement.get_measurement_args(self)
            measurement.measurement_method(arguments = args)
            
        
            
        now = datetime.now()
        now_str = now.strftime("%Y_%m_%d_%H_%M_%S")
        self.snapshot_name = 'Snapshot_'+now_str+'.txt'
        
        self.save_header = 'Snapshot: '+now_str+'\nNumber of devices: '+str(len(self.device_control_initial))+ '\nNumber of measurement devices: '+str(len(self.measururement_control)) + '\n'
        
        device_name_str = 'Number \t'
        for i in range(len(self.device_control_initial)):
            device_name_str += self.device_control_initial[i].name + '\t'
            
        for i in range(len(self.measururement_control)):
            for j in range(self.measururement_control[i].get_number_measurement_outputs()):
                device_name_str += self.measururement_control[i].name +'_'+str(j)+'\t'
        self.save_header += device_name_str + '\n'
        
        device_state = self.get_device_state()
        
        # check if directory for snapshots exists
        if not os.path.exists(self.snapshot_folder):
            os.makedirs(self.snapshot_folder)
        
        with open(self.snapshot_folder+'/'+self.snapshot_name, 'w') as file:
            file.write(self.save_header)
            file.write(str(0)+'\t')
            for device in device_state:
                file.write(str(device)+'\t')
            file.write('\n')
        
            
        
    
    def run_scan(self):
        counter = 0
        save_counter = 0
        
        # check if directory for hyper scans exists
        if not os.path.exists(self.hyper_scan_folder):
            os.makedirs(self.hyper_scan_folder)
        
        # open.txt file and write header
        with open(self.hyper_scan_folder+'/'+self.save_name, 'w') as file:
            file.write(self.save_header)

            # record run time 
            start_time = time.time()
            for i in range(len(self.total_scan)):
                if self.total_scan[i][0] == 0:
                    self.total_scan[i][1].set_control_value(self.total_scan[i][2])
                    counter += 1
                elif self.total_scan[i][0] == 1:
                    measurement_args = self.total_scan[i][1].get_measurement_args(self)
                    self.total_scan[i][1].measurement_method(arguments = measurement_args)
                    
                elif self.total_scan[i][0] == 2:
                    device_state = self.get_device_state()
                    file.write(str(save_counter)+'\t')
                    for device in device_state:
                        file.write(str(device)+'\t')
                    file.write('\n')
                    save_counter += 1
                    
                    
                if self.live_plot.get():
                    time_since_lp = time.time()-start_time
                    if time_since_lp >= self.live_plot_update_rate or i == len(self.total_scan)-1:
                        # set to read mode 
                        file.close()
                        #file = open(self.hyper_scan_folder+'/'+self.save_name, 'r')
                        self.update_live_plotting()
                        #file.close()
                        file = open(self.hyper_scan_folder+'/'+self.save_name, 'a')
                        start_time = time.time()    
                
                
                self.update_device_control()
                self.update_measurement_control()
                if self.open_gui:
                    self.scan_progress_bar['value'] = counter/self.total_scan_num
                    self.gui_window.update_idletasks()
        
        # save and close .txt file
            file.close()
            
            
        pass
                
    
    def update_device_control(self):
        for device in self.device_control:
            device.update_control()
            if device.open_gui:
                device.update_gui()
                #device.gui_window.update_idletasks()
    
    def update_measurement_control(self):
        for device in self.measururement_control_initial:
            if device.open_gui:
                #device.update_gui()
                #device.gui_window.update_idletasks()
                device.gui_window.update()
    
    def update_gui(self):
        for i in range(len(self.device_control_initial)):
            self.scan_limits[i][0] = float(self.scan_limit_entries_min[i].get())
            self.scan_limits[i][1] = float(self.scan_limit_entries_max[i].get())
            self.scan_num[i] = int(self.scan_num_entry[i].get())

        self.compile_button.config(state=tk.NORMAL)
        
    def button_color(self, button_object, state):
            if state:
                button_object.config(bg='green')
            else:
                button_object.config(bg='white')
    
    
    def snapshot_gui(self):
        self.set_snapshot_folder(self.snapshot_folder_entry.get())
        self.snapshot()
        self.snapshot_name_label.config(text = self.snapshot_name)
    
    def run_scan_gui(self):
        self.set_hyper_scan_folder(self.scan_folder_entry.get())
        self.save_name = self.save_name_entry.get()
        self.run_scan()
    
    def add_device_gui(self,device_index = 0):
        self.update_gui()
        self.add_device_to_scan(device_index)
        self.device_string_label.config(text = self.device_string)
    
    def set_hyper_scan_folder_gui(self):
        self.hyper_scan_folder = tk.filedialog.askdirectory()
        self.scan_folder_entry.delete(0, tk.END)
        self.scan_folder_entry.insert(0, self.hyper_scan_folder)
    
    def set_snapshot_folder_gui(self):
        self.snapshot_folder = tk.filedialog.askdirectory()
        self.snapshot_folder_entry.delete(0, tk.END)
        self.snapshot_folder_entry.insert(0, self.snapshot_folder)
    
    def add_measurement_to_scan_gui(self,measurement_index = 0):
        self.add_measurement_to_scan(measurement_index)
        self.measurement_string_label.config(text = self.measurement_string)
    
    def compile_scan_gui(self):
        self.run_button.config(state=tk.NORMAL)
        self.compile_scan()
        self.device_string_label.config(text = self.device_string)
        self.measurement_string_label.config(text = self.measurement_string)
        self.save_name_entry.config(state=tk.NORMAL)
        self.save_name_entry.delete(0, tk.END)
        self.save_name_entry.insert(0, self.save_name)
        #up to 3 dimensional scans can be live plotted 
        if len(self.device_control) <= 3:
            self.live_plot_checkbox.config(state=tk.DISABLED)
            self.live_plot_settings.config(state=tk.NORMAL)
        else:
            self.live_plot_checkbox.config(state=tk.DISABLED)
            self.live_plot.set(False)
            self.live_plot_settings.config(state=tk.DISABLED)
        
        
        #self.live_plot_settings.config(state=tk.NORMAL)
        
    def reset_scan_gui(self):
        self.reset_scan()
        self.run_button.config(state=tk.DISABLED)
        self.compile_button.config(state=tk.DISABLED)
        self.device_string_label.config(text = self.device_string)
        self.measurement_string_label.config(text = self.measurement_string)
        self.save_name_entry.config(state=tk.DISABLED)
        self.save_name_entry.delete(0, tk.END)
        self.save_name_entry.insert(0, self.save_name)
        #
        self.live_plot_checkbox.config(state=tk.DISABLED)
        self.live_plot_settings.config(state=tk.DISABLED)
        self.lp_data = []
        self.lp_windows = []
        
    def gui(self):
        self.open_gui = True
        if self.measururement_control_initial[0].open_gui:
            self.gui_window = tk.Toplevel(self.measururement_control_initial[0].gui_window)
        #self.gui_window = tk.Toplevel(self.measururement_control[0].gui_window)
        else:
            self.gui_window = tk.Tk()
            
        def close_gui():
            self.gui_window.destroy()
            print('Hyper Scan closed')
        
        self.gui_window.protocol("WM_DELETE_WINDOW", close_gui)
        self.gui_window.title('Hyper Dimensional Scan')
        
        self.snapshot_button = tk.Button(self.gui_window, text='Snapshot', command = self.snapshot_gui)
        self.snapshot_button.grid(row=0, column=0)
        
        self.run_button = tk.Button(self.gui_window, text='Run Scan', command = self.run_scan_gui,state=tk.DISABLED)
        self.run_button.grid(row=0, column=1)
        
        self.compile_button = tk.Button(self.gui_window, text='Compile Scan', command = self.compile_scan_gui, state=tk.DISABLED)
        self.compile_button.grid(row=0, column=2)
        
        self.reset_button = tk.Button(self.gui_window, text='Reset Scan', command = self.reset_scan_gui)
        self.reset_button.grid(row=0, column=3)
        
        self.scan_progress_bar = ttk.Progressbar(self.gui_window, orient='horizontal', length=200, mode='determinate', value=0,maximum=1)
        self.scan_progress_bar.grid(row=0, column=4, columnspan=2)
        
        tk.Label(self.gui_window, text='Device Control').grid(row=1, column=0)
        tk.Label(self.gui_window, text='Min').grid(row=1, column=1)
        tk.Label(self.gui_window, text='Max').grid(row=1, column=2)
        tk.Label(self.gui_window, text='Num').grid(row=1, column=3)
        device_row_offset = 2
        
        self.scan_limit_entries_min = []
        self.scan_limit_entries_max = []
        self.scan_num_entry = []
        self.scan_add_buttons = []
        
        for i, device in enumerate(self.device_control_initial):
            tk.Label(self.gui_window, text=device.name).grid(row=device_row_offset+i, column=0)
            
            self.scan_limit_entries_min.append(tk.Entry(self.gui_window))
            self.scan_limit_entries_min[i].grid(row=device_row_offset+i, column=1)
            self.scan_limit_entries_min[i].insert(0, str(self.scan_limits[i][0]))
            
            self.scan_limit_entries_max.append(tk.Entry(self.gui_window))
            self.scan_limit_entries_max[i].grid(row=device_row_offset+i, column=2)
            self.scan_limit_entries_max[i].insert(0, str(self.scan_limits[i][1]))
            
            self.scan_num_entry.append(tk.Entry(self.gui_window))
            self.scan_num_entry[i].grid(row=device_row_offset+i, column=3)
            self.scan_num_entry[i].insert(0, str(self.scan_num[i]))

            self.scan_add_buttons.append(tk.Button(self.gui_window, text='Add to Scan', command = lambda i=i: self.add_device_gui(i)))
            self.scan_add_buttons[i].grid(row=device_row_offset+i, column=4)
        
        tk.Label(self.gui_window, text='Snapshot folder: ').grid(row=device_row_offset+len(self.device_control_initial), column=0)
        
        self.snapshot_folder_entry = tk.Entry(self.gui_window, width=40)
        self.snapshot_folder_entry.grid(row=device_row_offset+len(self.device_control_initial), column=1, columnspan=2)
        self.snapshot_folder_entry.insert(0, self.snapshot_folder)
        
        set_snapshot_folder_button = tk.Button(self.gui_window, text='Search folder', command = self.set_snapshot_folder_gui)
        set_snapshot_folder_button.grid(row=device_row_offset+len(self.device_control_initial), column=3)
        
        self.snapshot_name_label = tk.Label(self.gui_window, text = self.snapshot_name)
        self.snapshot_name_label.grid(row=device_row_offset+len(self.device_control_initial), column=4, columnspan=2)
        
        tk.Label(self.gui_window, text='Scan folder: ').grid(row=device_row_offset+len(self.device_control_initial)+1, column=0)
        
        set_scan_folder_button = tk.Button(self.gui_window, text='Search folder', command = self.set_hyper_scan_folder_gui)
        set_scan_folder_button.grid(row=device_row_offset+len(self.device_control_initial)+1, column=3)
        
        self.scan_folder_entry = tk.Entry(self.gui_window, width=40)
        self.scan_folder_entry.grid(row=device_row_offset+len(self.device_control_initial)+1, column=1, columnspan=2)
        self.scan_folder_entry.insert(0, self.hyper_scan_folder)
        
        
        self.save_name_entry = tk.Entry(self.gui_window)
        self.save_name_entry.grid(row=device_row_offset+len(self.device_control_initial)+1, column=4, columnspan=2)
        self.save_name_entry.insert(0, self.save_name)
        self.save_name_entry.config(state=tk.DISABLED)
        
        tk.Label(self.gui_window, text='~~~~~~~Measurement Control~~~~~~~~').grid(row=device_row_offset+len(self.device_control_initial)+2, column=0, columnspan=1)
        
        
        self.live_plot = tk.BooleanVar(value=False)
        self.live_plot_checkbox = tk.Checkbutton(self.gui_window, text='Live Plot', variable = self.live_plot, state=tk.DISABLED)
        self.live_plot_checkbox.grid(row=device_row_offset+len(self.device_control_initial)+2, column=3)
        
        self.live_plot_settings = tk.Button(self.gui_window, text='Live Plot Settings', command = self.live_plot_settings_gui, state=tk.DISABLED)
        self.live_plot_settings.grid(row=device_row_offset+len(self.device_control_initial)+2, column=4)
        
        measurment_row_offset = device_row_offset + len(self.device_control_initial)+2
        
        measurment_row = 0
        self.measurement_add_buttons = []
        for i in range(len(self.measururement_control_initial)):
            rows = self.measururement_control_initial[i].measurement_method_gui(control = self, row_offset = measurment_row_offset+ measurment_row)
            measurment_row += rows
            tk.Label(self.gui_window, text='------------------------').grid(row=measurment_row_offset+measurment_row+1, column=0)
            
            self.measurement_add_buttons.append(tk.Button(self.gui_window, text='Add to Scan', command = lambda i=i: self.add_measurement_to_scan_gui(i)))
            self.measurement_add_buttons[i].grid(row=measurment_row_offset+measurment_row+1, column=1)
            
            
            
            measurment_row += 1

        self.device_string_label = tk.Label(self.gui_window, text = self.device_string)
        self.device_string_label.grid(row=measurment_row_offset+measurment_row+2, column=0)
        
        self.measurement_string_label = tk.Label(self.gui_window, text = self.measurement_string)
        self.measurement_string_label.grid(row=measurment_row_offset+measurment_row+2, column=2)
        
    def live_plot_settings_gui(self):
        self.lp_dimensions = np.min([len(self.device_control),3])
        self.lp_windows = []
        self.lp_data = []
        self.live_plot_update_rate = 3
        
        #self.live_plot_settings.config(state=tk.DISABLED)
        self.live_plot_settings_window = tk.Toplevel(self.gui_window)
        self.live_plot_settings_window.title('Live Plot Settings')
        #self.live_plot_settings_window.protocol("WM_DELETE_WINDOW", lambda: [self.live_plot_settings_window.withdraw])
        
        self.live_plot_values = []
        
        for i in range(len(self.device_control)):
            self.live_plot_values.append(self.device_control[i].name)
        
        for i in range(len(self.measururement_control)):
            for j in range(self.measururement_control[i].get_number_measurement_outputs()):
                self.live_plot_values.append(self.measururement_control[i].name+'_'+str(j))
        self.live_plot_values.append('Number')
        #self.live_plot_values.append('--None--')
        
        self.live_plot_axis_menu = []
        for i in range(self.lp_dimensions):
            tk.Label(self.live_plot_settings_window, text='--Axis '+str(i+1)+'--').grid(row=0, column=i)
            self.live_plot_axis_menu.append(ttk.Combobox(self.live_plot_settings_window, values = self.live_plot_values, state='readonly'))
            self.live_plot_axis_menu[i].grid(row=1, column=i)
            
            self.live_plot_axis_menu[i].current(i)

        tk.Label(self.live_plot_settings_window, text ='--Data--').grid(row=0, column=i+1)
        self.live_plot_data_menu = ttk.Combobox(self.live_plot_settings_window, values = self.live_plot_values, state='readonly')
        self.live_plot_data_menu.grid(row=1, column=i+1)
        self.live_plot_data_menu.current(len(self.device_control))
        
        tk.Label(self.live_plot_settings_window, text = 'Update rate (s):').grid(row=2, column=0)
        self.live_plot_update_rate_entry = tk.Entry(self.live_plot_settings_window)
        self.live_plot_update_rate_entry.grid(row=2, column=1)
        self.live_plot_update_rate_entry.insert(0, str(self.live_plot_update_rate))
        
        self.live_plot_add_plot_button = tk.Button(self.live_plot_settings_window, text='Add Plot', command = self.add_live_plot)
        self.live_plot_add_plot_button.grid(row=2, column=2)
        
        pass
    
    def add_live_plot(self):
        self.lp_windows.append(tk.Toplevel(self.gui_window))
        lp_current = self.lp_windows[-1]
        lp_current.title('Live Plot')
        
        def close_windows():
            lp_current.destroy()
            self.lp_windows.remove(lp_current)
            self.lp_data.remove(self.lp_data[-1])
            if len(self.lp_windows) == 0:
                self.live_plot.set(False)
                self.live_plot_checkbox.config(state=tk.DISABLED)
        
        lp_current.protocol("WM_DELETE_WINDOW", close_windows)
        
        lp_current.fig = Figure(figsize=(5, 4), dpi=100)
        lp_current.ax = lp_current.fig.add_subplot()
        
        lp_current.axis = []
        for i in range(self.lp_dimensions):
            lp_current.axis.append(self.live_plot_axis_menu[i].get()) # .current()
        lp_current.axis.append(self.live_plot_data_menu.get()) # .current()
        
        for i in range(self.lp_dimensions+1):
            self.lp_data.append(list(-1*np.ones(self.lp_dimensions+1)))
            
        lp_current.ax.set_xlabel(self.live_plot_axis_menu[0].get())
        
        print(self.lp_dimensions)
        print(self.lp_data)
        
        #print(self.lp_data[0])
        #print(self.lp_data[1])
        if self.lp_dimensions == 1:
            lp_current.ax.set_ylabel(self.live_plot_data_menu.get())
            lp_current.lp_plot_handle, = lp_current.ax.plot(self.lp_data[-1][0],self.lp_data[-1][1],'ko')
            lp_current.ax.set_xlim([0,1])
            lp_current.ax.set_ylim([0,1])
        elif self.lp_dimensions == 2:
            lp_current.ax.set_ylabel(self.live_plot_axis_menu[1].get())
            lp_current.lp_plot_handle = lp_current.ax.scatter(self.lp_data[-1][0],self.lp_data[-1][1],c=self.lp_data[-1][2],cmap='viridis')
            # colorbar
            lp_current.cbar = lp_current.fig.colorbar(lp_current.lp_plot_handle)
            lp_current.ax.set_xlim([0,1])
            lp_current.ax.set_ylim([0,1])
        elif self.lp_dimensions == 3:
            
            lp_current.lp_plot_handle = lp_current.ax.scatter(self.lp_data[-1][0],self.lp_data[-1][1],self.lp_data[-1][2],c=self.lp_data[-1][3],cmap='viridis')
            lp_current.ax.set_ylabel(self.live_plot_axis_menu[1].get())
            lp_current.ax.set_ylabel(self.live_plot_axis_menu[1].get())
            lp_current.ax.set_zlabel(self.live_plot_axis_menu[2].get())
            lp_current.ax.set_xlim([0,1])
            lp_current.ax.set_ylim([0,1])
            lp_current.ax.set_zlim([0,1])
            # colorbar
            lp_current.cbar = lp_current.fig.colorbar(lp_current.lp_plot_handle)
            
        
        lp_current.canvas = FigureCanvasTkAgg(lp_current.fig, master=lp_current)  # A tk.DrawingArea.
        lp_current.canvas.draw()
        lp_current.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        
        self.live_plot.set(True)
        self.live_plot_checkbox.config(state=tk.NORMAL)
        
        pass
    def initialise_live_plotting(self):
        self.lp_dimensions = np.min([len(self.device_control),3])
        self.lp_data = []
        self.live_plot_update_rate = 3
    
    def update_live_plotting(self):
        self.live_plot_update_rate = float(self.live_plot_update_rate_entry.get())
        # access date from the saved .txt file
        #print('updating live plotting')
        #with open(self.hyper_scan_folder+'/'+self.save_name, 'r') as file:
        
        ######
        # lines = file.readlines()
        # header_found = False
        # #self.lp_data = []
        # for line in lines:
        #     if line[0] == '0':
        #         header_found = True
            
        #     if header_found:
        #         data = line.split('\t')
        #         print(data)
                
        #         for lp_current in self.lp_windows:
        #             for i in range(self.lp_dimensions+1):
        #                 print(lp_current.axis[i])
        #                 self.lp_data[i].append(float(data[lp_current.axis[i]]))
        # print(self.lp_data)
        # if self.lp_dimensions == 1:
        #     lp_current.lp_plot_handle.set_xdata(self.lp_data[0])
        #     lp_current.lp_plot_handle.set_ydata(self.lp_data[1])
        self.lp_data = []
        for j, lp_current in enumerate(self.lp_windows):
            #print(lp_current.axis)
            header_number = 9 + len(self.device_control) + len(self.measururement_control)
            #print(header_number)
            data_file = pd.read_csv(self.hyper_scan_folder+'/'+self.save_name, delimiter = '\t',header=header_number, usecols=lp_current.axis,dtype=float)
            
            dat = []
            for i in range(self.lp_dimensions+1):
                dat.append(data_file[lp_current.axis[i]])
            self.lp_data.append(dat)
            
            lp_current.ax.set_xlim([np.min(self.lp_data[j][0]),np.max(self.lp_data[j][0])])
            lp_current.ax.set_ylim([np.min(self.lp_data[j][1]),np.max(self.lp_data[j][1])])
            
            if self.lp_dimensions == 1:
                lp_current.lp_plot_handle.set_xdata(self.lp_data[j][0])
                lp_current.lp_plot_handle.set_ydata(self.lp_data[j][1])
            elif self.lp_dimensions == 2:
                #print(self.lp_data[j][2])
                lp_current.lp_plot_handle.set_offsets(np.c_[self.lp_data[j][0],self.lp_data[j][1]])
                lp_current.lp_plot_handle.set_array(self.lp_data[j][2])
                lp_current.lp_plot_handle.set_clim([np.min(self.lp_data[j][2]),np.max(self.lp_data[j][2])])
                lp_current.cbar.update_normal(lp_current.lp_plot_handle)
            elif self.lp_dimensions == 3:
                lp_current.lp_plot_handle._offsets3d = (self.lp_data[j][0],self.lp_data[j][1],self.lp_data[j][2])
                lp_current.lp_plot_handle.set_array(self.lp_data[j][3])
                lp_current.lp_plot_handle.set_clim([np.min(self.lp_data[j][3]),np.max(self.lp_data[j][3])])
                lp_current.cbar.update_normal(lp_current.lp_plot_handle)
                
                
            
            lp_current.canvas.draw()
            lp_current.update()
                # for i in range(len(self.lp_windows)):
                #     lp_current = self.lp_windows[i]
                #     for j in range(self.lp_dimensions+1):
                #         self.lp_data[j] = data[j]
                    
                #     if self.lp_dimensions == 1:
                #         lp_current.lp_plot_handle.set_xdata(self.lp_data[0])
                #         lp_current.lp_plot_handle.set_ydata(self.lp_data[1])
                #     elif self.lp_dimensions == 2:
                #         lp_current.lp_plot_handle.set_xdata(self.lp_data[0])
                #         lp_current.lp_plot_handle.set_ydata(self.lp_data[1])
                #     elif self.lp_dimensions == 3:
                #         lp_current.lp_plot_handle.set_xdata(self.lp_data[0])
                #         lp_current.lp_plot_handle.set_ydata(self.lp_data[1])
                #         lp_current.lp_plot_handle.set_3d_properties(self.lp_data[2])
                    
                #     lp_current.ax.set_xlim([np.min(self.lp_data[0]),np.max(self.lp_data[0])])
                #     lp_current.ax.set_ylim([np.min(self.lp_data[1]),np.max(self.lp_data[1])])
                #     lp_current.canvas.draw()
                #     lp_current.canvas.flush_events()
                #     lp_current.update()
                #     lp_current.deiconify()
                #     lp_current.focus_force()
                    
                
           
                
                
        
       
        
        pass
    
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
    qd_calibration = 'QD_Iker_April_high_fss.txt'

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

    pulse_shaper = pulse_shaper_obj(device=lab_motor, calibration_file=ps_calibration, name='heisl')
    att_object = attenuator(lab_att, name = 'F Att')    
            
    
    pulse_shaper2 = pulse_shaper_obj(device=lab_motor2, calibration_file=ps_calibration, name='heisl2')
    att_object2 = attenuator(lab_att, name = 'F Att2')
    
    delay_stage = time_delay(device=lab_motor3, name='delay_stage')
    
    
    
    ps_controller = pulse_shaper.open_control(pulse_object=initial_pulse)
    att_controller = att_object.open_control(previous_control=ps_controller)
    delay_controller = delay_stage.open_control(previous_control=att_controller)

    simulator_object = simulator(qd_calibration=qd_calibration, sim_kind='ace', temp_dir='sim_dump/')
    sim_controller = simulator_object.open_control(previous_control = [att_controller,delay_controller], open_gui=True) 
      
    spec = spectrometer(device=lab_spectromter, name='MF_spec') 
    sm_controller = spec.open_control(pulse_object=None,previous_control = [delay_controller], simulation_control= sim_controller,open_gui=True) 
    
    
    lab_power_meter = fake_power_meter()
    
    power_meter_A = power_meter(device=lab_power_meter,name='fakeO')
    pm_controller = power_meter_A.open_control(previous_control=[att_controller,delay_controller])
    
    scan_controller = hyper_scan(device_control=[ps_controller,att_controller,delay_controller], measururement_control=[sm_controller, pm_controller], open_gui=True)
    
    # scan_controller.set_scan_limits(1, [0.3,0.7])
    # print(scan_controller.get_scan_limits(0))
    # print(scan_controller.get_scan_limits(1))
    # print(scan_controller.create_scan_vector(0))
    # scan_controller.reset_scan()
    # scan_controller.set_scan_num(0, 20)
    # scan_controller.set_scan_num(1, 30)
    # scan_controller.set_scan_num(2, 20)
    # scan_controller.add_device_to_scan(2)
    # scan_controller.add_device_to_scan(0)
    # scan_controller.add_device_to_scan(1)
    # scan_controller.add_measurement_to_scan(0)
    # scan_controller.compile_scan()
    # scan_controller.run_scan()
    # print(scan_controller.get_total_scan_num())
    #print(scan_controller.get_total_scan())
    
    scan_controller.start_gui()