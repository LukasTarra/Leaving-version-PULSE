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

import time

import threading

class simulator_control():
    def __init__(self,simulator_object,pulse_object = None, open_gui = True, display_simulation = True, previous_control = None)->None:
    
        self.simulator_object = simulator_object
        self.pulse_object = pulse_object
        self.open_gui = open_gui
        self.display_simulation = display_simulation
        self.previous_control = previous_control
        self.running = False
        self.simulator_kind = self.simulator_object.sim_kind
        self.dipole_moment = 1
        self.simulation_time = 0
        self.decay = self.simulator_object.decay
        if pulse_object is not None:
            if type(pulse_object) is str:
                    pulse_object = pg.load_pulse(pulse_object)
            out_pulse_object = pulse_object.copy_pulse()
            out_pulse_object.clear_filter()
            self.pulse_object = out_pulse_object
        else:
            self.pulse_object = None
        pass
        
        self.date_str_create = datetime.now().strftime("%Y_%m_%d_%H_%M")
        self.save_folder = 'simulation_'+self.date_str_create
        self.name = self.simulator_object.name
        
        self.force_gui_update = True
        self.dressed_states_flag = False
        
        self.print_info()
        self.update_previous_control()
        if self.open_gui:
            self.gui()
    
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
                #print(pulse_list)
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
    
    def set_simulation_kind(self,kind):
        self.simulator_kind = kind
        pass
    
    def set_dipole_moment(self,dipole_moment=1):
        self.dipole_moment = dipole_moment
    
    def set_display_simulation(self,display_simulation= True):
        self.display_simulation = display_simulation
        pass
    
    def set_force_gui_update(self,force_gui_update = True):
        self.force_gui_update = force_gui_update
        pass
    
    def set_save_folder(self,folder):
        self.save_folder = folder
    
    def toggle_dressed_states(self):
        self.dressed_states_flag = not self.dressed_states_flag
        pass
    
    def toggle_force_gui_update(self):
        self.force_gui_update = not self.force_gui_update
        pass
    
    def toggle_display_simulation(self):
        self.display_simulation = not self.display_simulation
        pass
    
    def toggle_running(self):
        self.running = not self.running
        pass
    
    def toggle_decay(self):
        self.decay = not self.decay
        self.simulator_object.set_decay(self.decay)
        
        pass
    
    def get_photon_emission_simulation(self):
        self.photon_wavelength = self.simulator_object.photon_wavelength
        self.photon_emission = self.simulator_object.photon_emission
        
    def get_simulation_time(self):
        return self.simulation_time
    
    def get_simulator_object(self): 
        return self.simulator_object
    
    def get_simulation_results(self):
        return self.simulation_result
    
    def update_simulation(self):
        start_time = time.time()
        self.update_previous_control()
        self.simulator_object.sim_kind = self.simulator_kind
        self.simulator_object.simulate(pulse_object = self.pulse_object,dipole_moment = self.dipole_moment)
        
        self.simulation_result = self.simulator_object.get_simulation_results()
        self.get_photon_emission_simulation()
        self.simulation_time = round(time.time()-start_time,3)
       #self.update_previous_control()
       # self.spectrometer_object.get_spectrum(pulse_object=self.pulse_object, sim_object=self.simulation_object,excecute = self.excecute)
        pass
    
    def update_simulation_thread(self):
        threading.Thread(target=self.update_simulation).start()
    
    def print_info(self):
        
        print('Simulator Control')
        print('Simulator Object: '+self.simulator_object.name)
        print(self.simulator_object.print_info())

        pass
    
    def button_color(self, button_object, state):
            if state:
                button_object.config(bg='green')
            else:
                button_object.config(bg='white')
    
    def update_gui(self):
        #self.folder_entry.delete(0,tk.END)
        #self.folder_entry.insert(0,self.save_folder)
        
        #self.dipole_moment = float(self.dipole_moment_entry.get())
        self.dipole_moment_text.config(text=str(self.dipole_moment))
        #self.update_simulation()
        self.pulse_x_plot.set_xdata(self.pulse_object.time)
        self.pulse_x_plot.set_ydata(np.real(self.pulse_object.temporal_representation_x))
        self.pulse_y_plot.set_xdata(self.pulse_object.time)
        self.pulse_y_plot.set_ydata(np.real(self.pulse_object.temporal_representation_y))
        
        # update x and y limits for pulses
        self.ax_pulses.set_xlim([np.min(self.pulse_object.time),np.max(self.pulse_object.time)])
        
        self.ax_pulses.set_ylim([1.1*np.min([np.real(self.pulse_object.temporal_representation_x),np.real(self.pulse_object.temporal_representation_y)]),1.1*np.max([np.abs(self.pulse_object.temporal_representation_x),np.abs(self.pulse_object.temporal_representation_y)])])
        
        self.button_color(self.run_button,self.running)
        
        # self.simulation_g_plot.set_xdata(self.simulation_result[0])
        # self.simulation_g_plot.set_ydata(self.simulation_result[1])
        # self.simulation_x_plot.set_xdata(self.simulation_result[0])
        # self.simulation_x_plot.set_ydata(self.simulation_result[2])
        # self.simulation_y_plot.set_xdata(self.simulation_result[0])
        # self.simulation_y_plot.set_ydata(self.simulation_result[3])
        # self.simulation_b_plot.set_xdata(self.simulation_result[0])
        # self.simulation_b_plot.set_ydata(self.simulation_result[4])
        
        if self.dressed_states_flag:
            self.update_states_gui()
        
        
        if self.display_simulation:
            self.update_simulation()
            if self.simulator_object.get_num_states() != self.old_num_states:
                self.old_num_states = self.simulator_object.get_num_states()
                self.plot_states_vec = []
                self.ax_simulation.clear()
                self.ax_simulation = self.ax_pulses.twinx()
                for i in range(self.simulator_object.get_num_states()):
                    cur_plot, = self.ax_simulation.plot([],[],label='State '+str(i))
                    if i == 0:
                        cur_plot.set_color('k')
                        cur_plot.set_alpha(0.75)
                    self.plot_states_vec.append(cur_plot)
                self.ax_simulation.set_ylabel('Population')
                self.ax_simulation.legend()
                self.ax_simulation.set_ylim([0,1])
            if self.running:
                for i in range(len(self.plot_states_vec)):
                    self.plot_states_vec[i].set_xdata(self.simulation_result[0])
                    self.plot_states_vec[i].set_ydata(self.simulation_result[i+1])
        else:
            self.update_previous_control()  
            if self.running:
                for i in range(len(self.plot_states_vec)):
                    self.plot_states_vec[i].set_xdata([])
                    self.plot_states_vec[i].set_ydata([])
        
        
        
        
        if self.running:
            #wait_time = self.simulation_time_entry.get()
            self.set_simulation_kind(self.simulation_kind_entry.get())
            #if int(wait_time) < int(np.round(1e3*self.simulation_time*1.1)):
            #    wait_time = int(np.round(1e3*self.simulation_time*1.1))
            
            # self.simulation_time_entry.delete(0,'end')
            # self.simulation_time_entry.insert(0,str(int(self.simulation_time*1e3)))
            self.simulation_time_text.config(text=str(int(self.simulation_time*1e3))+'ms')
            self.canvas.draw()
            
            if self.force_gui_update:
                #self.gui_window.after(int(np.round(1e3*self.simulation_time*1.1)), self.update_gui)
                self.gui_window.after(200,self.update_gui)
        else:
            #wait_time = 0
            #self.simulation_time_entry.delete(0,'end')
            #self.simulation_time_entry.insert(0,str(int(0)))
            self.simulation_time_text.config(text='off')
    
        
        
    
    def gui(self):
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
            
        
        
        # else:
        #     if type(self.previous_control) is list:
        #         self.gui_window = tk.Toplevel(self.previous_control[0].gui_window)
        #     else:
        #         self.gui_window = tk.Toplevel(self.previous_control.gui_window)
        
        self.gui_window.title('Simulator Control: '+self.simulator_object.name)
        #self.gui_window.geometry("500x520")
        
        self.set_display_simulation()
        self.update_previous_control()
        
        fig_simulator = Figure(figsize=(5, 4), dpi=100)
        
        self.ax_pulses = fig_simulator.add_subplot()
        self.ax_simulation = self.ax_pulses.twinx()
        
        self.pulse_x_plot, = self.ax_pulses.plot(self.pulse_object.time,np.real(self.pulse_object.temporal_representation_x),'b',alpha=0.3)
        self.pulse_y_plot, = self.ax_pulses.plot(self.pulse_object.time,np.real(self.pulse_object.temporal_representation_y),'r',alpha=0.3)
        
        self.ax_pulses.set_xlabel('Time (ps)')
        
        # self.simulation_g_plot, = self.ax_simulation.plot([],[],'k',alpha=0.75,label='g')
        # self.simulation_x_plot, = self.ax_simulation.plot([],[],'b',alpha=1,label='x')
        # self.simulation_y_plot, = self.ax_simulation.plot([],[],'r',alpha=1,label='y')
        # self.simulation_b_plot, = self.ax_simulation.plot([],[],'m',alpha=1,label='b')
        
        self.old_num_states = self.simulator_object.get_num_states()
        self.plot_states_vec = []
        for i in range(self.simulator_object.get_num_states()):
            cur_plot, = self.ax_simulation.plot([],[],label='State '+str(i))
            if i == 0:
                cur_plot.set_color('k')
                cur_plot.set_alpha(0.75)
            self.plot_states_vec.append(cur_plot)
        
        self.ax_simulation.set_ylabel('Population')
        self.ax_simulation.legend()
        self.ax_simulation.set_ylim([0,1])

        self.canvas = FigureCanvasTkAgg(fig_simulator, master=self.gui_window)  # A tk.DrawingArea.
        self.canvas.draw()
        
        self.canvas.get_tk_widget().grid(row=0,column=0,columnspan=3)
        
        self.run_button = tk.Button(self.gui_window, text="Run") 
        self.run_button.config(command= lambda: [self.toggle_running(),self.button_color(self.run_button,self.running),self.set_force_gui_update(),self.update_gui()]) 
        self.run_button.grid(row=1,column=0) 
        
        display_simulation_button = tk.Button(self.gui_window, text="Display Simulation")
        display_simulation_button.config(command= lambda: [self.toggle_display_simulation(),self.button_color(display_simulation_button,self.display_simulation)])
        display_simulation_button.grid(row=1,column=1)
        
        decay_button = tk.Button(self.gui_window, text="Decay" )
        decay_button.config(command= lambda: [self.toggle_decay(),self.button_color(decay_button,self.decay)])
        decay_button.grid(row=1,column=2)
        
        # color buttons 
        self.button_color(self.run_button,self.running)
        self.button_color(decay_button,self.decay)
        self.button_color(display_simulation_button,self.display_simulation)
        
        # self.simulation_time_entry = tk.Entry(master = self.gui_window, width = 10, state='normal')
        # self.simulation_time_entry.insert(0,str(self.simulation_time))
        # self.simulation_time_entry.grid(row=2,column=0)
        
        self.simulation_time_text = tk.Label(master = self.gui_window, text = str(self.simulation_time)+'ms')
        self.simulation_time_text.grid(row=2,column=0)
        
        tk.Label(master=self.gui_window, text = 'Simulation kind:').grid(row=2,column=1)
        
        # self.simulation_kind_entry = tk.Entry(master = self.gui_window, width = 10)
        # self.simulation_kind_entry.insert(0,self.simulator_kind)
        # self.simulation_kind_entry.grid(row=2,column=2)
        
        self.simulation_kind_entry = ttk.Combobox(self.gui_window, values=['ace','qutip','ace_6ls'],state='readonly')
        self.simulation_kind_entry.set(self.simulator_kind)
        self.simulation_kind_entry.grid(row=2,column=2)
    
        dipole_moment_label = tk.Label(master = self.gui_window, text = 'Dipole Moment:').grid(row=3,column=0)
        
        # self.dipole_moment_entry = tk.Entry(master = self.gui_window, width = 10)
        # self.dipole_moment_entry.insert(0,str(self.dipole_moment))
        # self.dipole_moment_entry.grid(row=3,column=1)
        
        self.dipole_moment_text = tk.Label(master = self.gui_window, text = str(self.dipole_moment))
        self.dipole_moment_text.grid(row=3,column=1)
        
        tk.Label(master = self.gui_window, text = 'Quantum dot file:').grid(row=4,column=0)
        self.qd_calibration_label =  tk.Label(master = self.gui_window, text = self.simulator_object.qd_calibration)
        self.qd_calibration_label.grid(row=4,column=1)
        
        self.change_qd_calibration_button = tk.Button(self.gui_window, text="Change QD Calibration")
        self.change_qd_calibration_button.config(command= lambda: [self.change_qd_calibration()])
        self.change_qd_calibration_button.grid(row=4,column=2)
        
        self.save_simulation_button = tk.Button(self.gui_window, text='Save simulation')
        self.save_simulation_button.config(command= lambda: self.save_simulation_gui())
        self.save_simulation_button.grid(row=5,column=0) 
        
        tk.Label(self.gui_window, text='Folder: ').grid(row=5,column=1)
        
        self.folder_entry = tk.Entry(self.gui_window, width=30)
        self.folder_entry.insert(0,self.save_folder)
        self.folder_entry.grid(row=5,column=2)
        
        self.change_folder_button = tk.Button(self.gui_window, text='Change Folder')
        self.change_folder_button.config(command= lambda: [self.change_folder_gui(),self.update_gui()])
        self.change_folder_button.grid(row=6,column=2)
        
        self.advanced_options_button = tk.Button(self.gui_window, text='Advanced Options')
        self.advanced_options_button.config(command= lambda: [self.advanced_options_gui()])
        self.advanced_options_button.grid(row=6,column=0)
        
        self.dressed_states_button = tk.Button(self.gui_window, text='Energy view')
        self.dressed_states_button.config(command= lambda: [self.toggle_dressed_states(),self.dressed_states_gui(),self.button_color(self.dressed_states_button,self.dressed_states_flag)])
        self.dressed_states_button.grid(row=6,column=1)
        self.button_color(self.dressed_states_button,self.dressed_states_flag)
        
        if self.running:
            self.update_gui()

    def dressed_states_gui(self):
        states_window = tk.Toplevel(self.gui_window)
        states_window.title('Dressed States (Beta)')
        
        def close_window():
            self.dressed_states_flag = False
            self.button_color(self.dressed_states_button,self.dressed_states_flag)
            states_window.destroy()
        
        states_window.protocol("WM_DELETE_WINDOW", close_window)
        
        states_fig = Figure(figsize=(5, 4), dpi=100)
        self.ax_states = states_fig.add_subplot()
        
        self.states_plot = []
        for i in range(self.simulator_object.get_num_states()):
            cur_plot, = self.ax_states.plot([],[],label='State '+str(i))
            if i == 0:
                cur_plot.set_color('k')
                cur_plot.set_alpha(0.75)
            self.states_plot.append(cur_plot)
        
        self.ax_states.set_xlabel('Time (ps)')
        self.ax_states.set_ylabel('Energy (meV)')
        self.ax_states.legend()
        
        self.canvas_states = FigureCanvasTkAgg(states_fig, master=states_window)  # A tk.DrawingArea.
        self.canvas_states.draw()
        self.canvas_states.get_tk_widget().grid(row=0,column=0,rowspan=self.simulator_object.get_num_states())
        
        self.states_plot_vec = []
        self.states_plot_button_vec = []
        
        def toggle_state_plot(i):
            self.states_plot_vec[i] = not self.states_plot_vec[i]
        
        for i in range(self.simulator_object.get_num_states()):
            self.states_plot_vec.append(True)
            cur_button = tk.Button(states_window, text='State '+str(i))
            cur_button.config(command= lambda i=i: [toggle_state_plot(i),self.button_color(self.states_plot_button_vec[i],self.states_plot_vec[i])])
            cur_button.grid(row=i,column=1)
            self.button_color(cur_button,self.states_plot_vec[i])
            self.states_plot_button_vec.append(cur_button)
            
        self.update_states_gui()
    
    def update_states_gui(self):
        plot_time, plot_energy = self.simulator_object.light_dressed_states(self.pulse_object)
        min_energy = 0
        max_energy = 0
        for i in range(self.simulator_object.get_num_states()):
            if self.states_plot_vec[i]:
                self.states_plot[i].set_xdata(plot_time)
                self.states_plot[i].set_ydata(plot_energy[i])
                if np.min(plot_energy[i]) < min_energy:
                    min_energy = np.min(plot_energy[i])
                if np.max(plot_energy[i]) > max_energy:
                    max_energy = np.max(plot_energy[i])
            else:
                self.states_plot[i].set_xdata([])
                self.states_plot[i].set_ydata([])
        
        self.ax_states.set_xlim([np.min(plot_time),np.max(plot_time)])
        self.ax_states.set_ylim([min_energy,max_energy])
        self.canvas_states.draw()
    
    def change_qd_calibration(self):
        config = configparser.ConfigParser()
        config.read(self.simulator_object.qd_calibration) 
        exciton_wavelength = config['EMISSION']['exciton_wavelength']
        biexciton_wavelength = config['EMISSION']['biexciton_wavelength']
        dark_wavelength = config['EMISSION']['dark_wavelength']
        
        fss_bright = config['SPLITTING']['fss_bright']
        fss_dark = config['SPLITTING']['fss_dark']
        
        exciton_lt = config['LIFETIMES']['exciton']
        biexciton_lt = config['LIFETIMES']['biexciton']
        
        g_ex = config['G_FACTORS']['g_ex']
        g_ez = config['G_FACTORS']['g_ez']
        g_hx = config['G_FACTORS']['g_hx']
        g_hz = config['G_FACTORS']['g_hz']
        
        self.qd_window = tk.Toplevel(self.gui_window)
        self.qd_window.title('Quantum Dot Calibration')
        #self.qd_window.protocol("WM_DELETE_WINDOW", self.live_plot_settings_window.withdraw)
        
        tk.Label(master=self.qd_window, text= 'File name: ').grid(row=0,column=0)
        self.qd_calibration_entry = tk.Entry(master=self.qd_window, width = 20)
        self.qd_calibration_entry.insert(0,self.simulator_object.qd_calibration.split('/')[-1])
        self.qd_calibration_entry.grid(row=0,column=1)
        
        tk.Label(master=self.qd_window, text= '----').grid(row=1,column=0) 
    
        tk.Label(master=self.qd_window, text= 'Exciton Wavelength: ').grid(row=2,column=0)
        self.exciton_wavelength_entry = tk.Entry(master=self.qd_window, width = 20)
        self.exciton_wavelength_entry.insert(0,exciton_wavelength)
        self.exciton_wavelength_entry.grid(row=2,column=1)
        
        tk.Label(master=self.qd_window, text= 'Biexciton Wavelength: ').grid(row=3,column=0)
        self.biexciton_wavelength_entry = tk.Entry(master=self.qd_window, width = 20)
        self.biexciton_wavelength_entry.insert(0,biexciton_wavelength)
        self.biexciton_wavelength_entry.grid(row=3,column=1)
        
        tk.Label(master=self.qd_window, text= 'Dark Wavelength: ').grid(row=4,column=0)
        self.dark_wavelength_entry = tk.Entry(master=self.qd_window, width = 20)
        self.dark_wavelength_entry.insert(0,dark_wavelength)
        self.dark_wavelength_entry.grid(row=4,column=1)
        
        tk.Label(master=self.qd_window, text= '----').grid(row=5,column=0)
        
        tk.Label(master=self.qd_window, text= 'FSS Bright: ').grid(row=6,column=0)
        self.fss_bright_entry = tk.Entry(master=self.qd_window, width = 20)
        self.fss_bright_entry.insert(0,fss_bright)
        self.fss_bright_entry.grid(row=6,column=1)
        
        tk.Label(master=self.qd_window, text= 'FSS Dark: ').grid(row=7,column=0)
        self.fss_dark_entry = tk.Entry(master=self.qd_window, width = 20)
        self.fss_dark_entry.insert(0,fss_dark)
        self.fss_dark_entry.grid(row=7,column=1)
        
        tk.Label(master=self.qd_window, text= '----').grid(row=8,column=0)
        
        tk.Label(master=self.qd_window, text= 'Exciton Lifetime: ').grid(row=9,column=0)
        self.exciton_lt_entry = tk.Entry(master=self.qd_window, width = 20)
        self.exciton_lt_entry.insert(0,exciton_lt)
        self.exciton_lt_entry.grid(row=9,column=1)
        
        tk.Label(master=self.qd_window, text= 'Biexciton Lifetime: ').grid(row=10,column=0)
        self.biexciton_lt_entry = tk.Entry(master=self.qd_window, width = 20)
        self.biexciton_lt_entry.insert(0,biexciton_lt)
        self.biexciton_lt_entry.grid(row=10,column=1)
        
        tk.Label(master=self.qd_window, text= '----').grid(row=11,column=0)
        
        tk.Label(master=self.qd_window, text= 'g_ex: ').grid(row=12,column=0)
        self.g_ex_entry = tk.Entry(master=self.qd_window, width = 20)
        self.g_ex_entry.insert(0,g_ex)
        self.g_ex_entry.grid(row=12,column=1)
        
        tk.Label(master=self.qd_window, text= 'g_ez: ').grid(row=13,column=0)
        self.g_ez_entry = tk.Entry(master=self.qd_window, width = 20)
        self.g_ez_entry.insert(0,g_ez)
        self.g_ez_entry.grid(row=13,column=1)
        
        tk.Label(master=self.qd_window, text= 'g_hx: ').grid(row=14,column=0)
        self.g_hx_entry = tk.Entry(master=self.qd_window, width = 20)
        self.g_hx_entry.insert(0,g_hx)
        self.g_hx_entry.grid(row=14,column=1)
        
        tk.Label(master=self.qd_window, text= 'g_hz: ').grid(row=15,column=0)
        self.g_hz_entry = tk.Entry(master=self.qd_window, width = 20)
        self.g_hz_entry.insert(0,g_hz)
        self.g_hz_entry.grid(row=15,column=1)
        
        tk.Label(master=self.qd_window, text= '----').grid(row=16,column=0)
        
        tk.Label(master=self.qd_window, text='Dipole moment: ').grid(row=17,column=0)
        self.dipole_moment_entry = tk.Entry(master=self.qd_window, width = 20)
        self.dipole_moment_entry.insert(0,self.dipole_moment)
        self.dipole_moment_entry.grid(row=17,column=1)
        
        self.load_qd_calibration_button = tk.Button(self.qd_window, text="Load Calibration")
        self.load_qd_calibration_button.config(command= lambda: [self.load_qd_calibration()])
        self.load_qd_calibration_button.grid(row=18,column=0)
        
        self.save_qd_calibration_button = tk.Button(self.qd_window, text="Save Calibration")
        self.save_qd_calibration_button.config(command= lambda: [self.save_qd_calibration()])
        self.save_qd_calibration_button.grid(row=18,column=1)
        
    def load_qd_calibration(self):
        filepath = tk.filedialog.askopenfilename()
        # extract only the filename
        filepath = filepath.split('/')[-1]
        self.simulator_object.set_qd_calibration(filepath)
        self.qd_calibration_label.config(text=filepath)
        self.qd_window.destroy()
        self.change_qd_calibration()
        
    def save_qd_calibration(self):
        filepath = os.getcwd()+'/'+self.qd_calibration_entry.get()
        config = configparser.ConfigParser()
        config['EMISSION'] = {'exciton_wavelength': self.exciton_wavelength_entry.get(),
                             'biexciton_wavelength': self.biexciton_wavelength_entry.get(),
                             'dark_wavelength': self.dark_wavelength_entry.get()}
        
        config['SPLITTING'] = {'fss_bright': self.fss_bright_entry.get(),
                             'fss_dark': self.fss_dark_entry.get()}
        
        config['LIFETIMES'] = {'exciton': self.exciton_lt_entry.get(),
                             'biexciton': self.biexciton_lt_entry.get()}
        
        config['G_FACTORS'] = {'g_ex': self.g_ex_entry.get(),
                                'g_ez': self.g_ez_entry.get(),
                                'g_hx': self.g_hx_entry.get(),
                                'g_hz': self.g_hz_entry.get()}
        
        self.set_dipole_moment(float(self.dipole_moment_entry.get()))
        
        with open(filepath, 'w') as configfile:
            config.write(configfile)
        
        self.simulator_object.set_qd_calibration(filepath)
    
    def save_simulation_gui(self):
        self.set_save_folder(self.folder_entry.get())
        self.save_simulation()
    
    def save_simulation(self):
        now = datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
        save_name = self.save_folder+'/'+'simulation_'+now+'.txt'
        save_name_pulse = self.save_folder+'/'+'pulse_'+now+'.txt'
        
        if not os.path.exists(self.save_folder):
            os.makedirs(self.save_folder)
            
        with open(save_name, 'w', newline='') as file:
            writer = csv.writer(file)
            header_str = ['Time (ps)']
            for i in range(self.simulator_object.get_num_states()):
                header_str.append('state_'+str(i))
            writer.writerow(header_str)
            for i in range(len(self.simulation_result[0])):
                row = [self.simulation_result[0][i]]
                for j in range(self.simulator_object.get_num_states()):
                    row.append(self.simulation_result[j+1][i])
                writer.writerow(row)
        file.close()
                
        with open(save_name_pulse, 'w', newline='') as file:
            writer = csv.writer(file)
            header_str = ['Time (ps)','pulse_h','pulse_v']
            writer.writerow(header_str)
            for i in range(len(self.pulse_object.time)):
                row = [self.pulse_object.time[i],self.pulse_object.temporal_representation_x[i],self.pulse_object.temporal_representation_y[i]]
                writer.writerow(row)
        file.close()
                
                
    def change_folder_gui(self):
        self.save_folder = tk.filedialog.askdirectory()
        self.set_save_folder(self.save_folder)
        self.folder_entry.delete(0,tk.END)
        self.folder_entry.insert(0,self.save_folder)
        self.update_gui()
        pass
    
    def advanced_options_gui(self):
        
        def update_emission():
            #self.update_simulation()
            self.simulator_object.simulate(pulse_object = self.pulse_object,dipole_moment = self.dipole_moment)
            self.get_photon_emission_simulation()
            emission_label = ['X_H->G', 'X_V->G', 'B->X_H', 'B->X_V', 'D_H->G', 'D_V->G','B->D_H','B->D_V']
            tk.Label(master=self.advanced_window, text = 'Photon emission / nm').grid(row=0,column=3)
            tk.Label(master=self.advanced_window, text = 'Lifetime / ps').grid(row=0,column=4)
            for i in range(len(self.photon_wavelength)):
                tk.Label(master=self.advanced_window, text = emission_label[i]).grid(row=i+1,column=2)
                tk.Label(master=self.advanced_window, text = str(np.round(self.photon_wavelength[i],decimals=3))).grid(row=i+1,column=3)
                if i ==0:
                    tk.Label(master=self.advanced_window, text = str(np.round(self.pulse_object.lifetime_exciton*1/(1-self.simulator_object.decay_scale_x),decimals=3))).grid(row=i+1,column=4)
                elif i == 1:
                    tk.Label(master=self.advanced_window, text = str(np.round(self.pulse_object.lifetime_exciton*1/(1-self.simulator_object.decay_scale_y),decimals=3))).grid(row=i+1,column=4)
                elif i == 2:
                    tk.Label(master=self.advanced_window, text = str(np.round(self.pulse_object.lifetime_biexciton*1/(1-self.simulator_object.decay_scale_x),decimals=3))).grid(row=i+1,column=4)
                elif i == 3:
                    tk.Label(master=self.advanced_window, text = str(np.round(self.pulse_object.lifetime_biexciton*1/(1-self.simulator_object.decay_scale_y),decimals=3))).grid(row=i+1,column=4)
                elif i == 4:
                    if self.simulator_object.decay_scale_x == 0:
                        tk.Label(master=self.advanced_window, text = 'inf').grid(row=i+1,column=4)
                    else:
                        tk.Label(master=self.advanced_window, text = str(np.round(self.pulse_object.lifetime_exciton*1/(self.simulator_object.decay_scale_x),decimals=3))).grid(row=i+1,column=4)
                elif i == 5:
                    if self.simulator_object.decay_scale_y == 0:
                        tk.Label(master=self.advanced_window, text = 'inf').grid(row=i+1,column=4)
                    else:
                        tk.Label(master=self.advanced_window, text = str(np.round(self.pulse_object.lifetime_exciton*1/(self.simulator_object.decay_scale_y),decimals=3))).grid(row=i+1,column=4)
                elif i == 6:
                    if self.simulator_object.decay_scale_x == 0:
                        tk.Label(master=self.advanced_window, text = 'inf').grid(row=i+1,column=4)
                    else:
                        tk.Label(master=self.advanced_window, text = str(np.round(self.pulse_object.lifetime_biexciton*1/(self.simulator_object.decay_scale_x),decimals=3))).grid(row=i+1,column=4)
                elif i == 7:
                    if self.simulator_object.decay_scale_y == 0:
                        tk.Label(master=self.advanced_window, text = 'inf').grid(row=i+1,column=4)
                    else:
                        tk.Label(master=self.advanced_window, text = str(np.round(self.pulse_object.lifetime_biexciton*1/(self.simulator_object.decay_scale_y),decimals=3))).grid(row=i+1,column=4)
                        
            
        def update_advanced_options():
            self.simulator_object.set_dipole_orientation(float(self.dipole_orientation_entry.get()))
            self.simulator_object.set_temperature(float(self.temperature_entry.get()))
            
            if self.simulator_object.get_num_states() == 6:
                self.simulator_object.set_mag_field(b_x = float(self.mag_field_x_entry.get()),b_z = float(self.mag_field_z_entry.get()))
                
            for label in self.advanced_window.grid_slaves():
                if int(label.grid_info()['column']) >= 2:
                    label.grid_forget()
            update_emission()
            pass
        
        self.advanced_window = tk.Toplevel(self.gui_window)
        self.advanced_window.title('Advanced Options')
        
        self.update_advanced_options_button = tk.Button(master=self.advanced_window, text = 'Update')
        self.update_advanced_options_button.config(command= lambda: [update_advanced_options()])
        self.update_advanced_options_button.grid(row=0,column=0)
        
        tk.Label(master=self.advanced_window, text = 'Dipole orientation (deg): ').grid(row=1,column=0)
        
        self.dipole_orientation_entry = tk.Entry(master=self.advanced_window, width = 20)
        self.dipole_orientation_entry.insert(0,self.simulator_object.dipole_orientation)
        self.dipole_orientation_entry.grid(row=1,column=1)
        
        tk.Label(master=self.advanced_window, text = 'Toggle phonons: ').grid(row=2,column=0)
        
        self.toggle_phonons_button = tk.Button(master=self.advanced_window, text = 'Toggle')
        self.toggle_phonons_button.config(command= lambda: [self.simulator_object.toggle_phonons(),self.button_color(self.toggle_phonons_button,self.simulator_object.phonons)])
        self.toggle_phonons_button.grid(row=2,column=1)
        self.button_color(self.toggle_phonons_button,self.simulator_object.phonons)
        
        tk.Label(master=self.advanced_window, text = 'Temperature / K').grid(row=3,column=0)
        
        self.temperature_entry = tk.Entry(master=self.advanced_window, width = 20)
        self.temperature_entry.insert(0,self.simulator_object.temperature)
        self.temperature_entry.grid(row=3,column=1)
        
        
        update_emission()
        
        # enable when phonons are implemented
        if not self.simulator_object.phonons:
            self.toggle_phonons_button.config(state='disabled')
            self.temperature_entry.config(state='disabled')
        
        num_general_settings = 4
        if self.simulator_object.get_num_states() == 6:
            tk.Label(self.advanced_window, text = '~~~~ Six Level Settings ~~~~').grid(row=num_general_settings,column=0)
            tk.Label(master=self.advanced_window, text = 'Magnetic field Z (T): ').grid(row=num_general_settings+1,column=0)
            
            self.mag_field_z_entry = tk.Entry(master=self.advanced_window, width = 20)
            self.mag_field_z_entry.insert(0,self.simulator_object.b_z)
            self.mag_field_z_entry.grid(row=num_general_settings+1,column=1)
            
            tk.Label(master=self.advanced_window, text = 'Magnetic field X (T): ').grid(row=num_general_settings+2,column=0)
            
            self.mag_field_x_entry = tk.Entry(master=self.advanced_window, width = 20)
            self.mag_field_x_entry.insert(0,self.simulator_object.b_x)
            self.mag_field_x_entry.grid(row=num_general_settings+2,column=1)
            
            tk.Label(master=self.advanced_window, text = 'Transorm in mag. frame: ').grid(row=num_general_settings+3,column=0)
            
            self.toggle_mag_frame_button = tk.Button(master=self.advanced_window, text = 'Toggle')
            self.toggle_mag_frame_button.config(command= lambda: [self.simulator_object.toggle_b_field_frame(),self.button_color(self.toggle_mag_frame_button,self.simulator_object.b_field_frame)])
            self.toggle_mag_frame_button.grid(row=num_general_settings+3,column=1)
            self.button_color(self.toggle_mag_frame_button,self.simulator_object.b_field_frame)
        
        
        pass
    
    def start_gui(self):
        self.gui_window.mainloop()
        
        
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
    qd_calibration2 = 'QD_Iker_April_2.txt'
    ps_calibration2 = 'calibration_slit_38.txt' 
    ps_calibration = 'calibration_Cam_230424.txt' 
    
    t_0 = 0 
    t_end = 20
    
    initial_pulse = pg.PulseGenerator(t_0,t_end,0.2,calibration_file=qd_calibration)
    #initial_pulse.add_gaussian_time(unit = 'nm', central_f= 779.89, width_t=0.5,t0=t_end/2,area_time=10)
    #initial_pulse.plot_pulses(domain='nm')
    initial_pulse.add_gaussian_time(unit = 'nm', central_f= 779.89, width_t=0.1,t0=10,area_time=4)
    
    
    lab_motor = None #motor(fake_motor(),name='motor1')
    lab_motor2 = None # motor(fake_motor(),name='motor2')
    
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
    
    simulator_object = simulator(qd_calibration=qd_calibration, sim_kind='ACE', temp_dir='\sim_dump')
    sim_controller = simulator_object.open_control(previous_control = [att_controller]) 
    
    simulator_object_2 = simulator(qd_calibration=qd_calibration2, sim_kind='ACE', temp_dir='\sim_dump')
    sim_controller_2 = simulator_object_2.open_control(previous_control = [att_controller]) 
    
    lab_spectromter = fake_spectrometer(start_wl=774, end_wl = 790, n_wl = 1024)

    spec = spectrometer(device=lab_spectromter, name='MF_spec') 
    sm_controller = spec.open_control(pulse_object=None,previous_control = [att_controller], simulation_control= [sim_controller, sim_controller_2]) #, ps_controller3#ps_controller2.get_pulse_object()
    
    #sm_controller.print_info()
    ps_controller.start_gui()
    
    pass