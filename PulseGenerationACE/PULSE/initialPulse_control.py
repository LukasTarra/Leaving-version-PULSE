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


class Pulse_Generator():
    def __init__(self, t0 = 0, t_end = 100, dt = 0.05, central_wavelength = 800, open_gui = True, initial_pulse = None, name = 'Pulse_generator') -> None:
        self.pulse_list = []
        self.pulse_args = []
        self.pulse_pol = []
        self.pulse_powers = []
        self.pulse_phases = []
        self.pulse_kind = []
        self.name = name
        
        self.open_gui = open_gui
        
        if initial_pulse is None:
            self.central_wavelength = central_wavelength
            self.t0 = t0
            self.t_end = t_end
            self.dt = dt
            
            self.pulse_object = pg.PulseGenerator(t0=self.t0, tend=self.t_end, dt=self.dt, central_wavelength=self.central_wavelength)
        else:
            self.central_wavelength = initial_pulse.central_wavelength
            self.t0 = initial_pulse.t0
            self.t_end = initial_pulse.tend
            self.dt = initial_pulse.dt
            
            self.pulse_object = initial_pulse.copy_pulse()
        
        self.initial_pulse_object = initial_pulse
        
        
        self.pulse_kinds = ['gaussian_time',
                            'gaussian_frequency']
        
        if self.open_gui:
            self.gui()
        pass
    
    def save_pulse(self,name):
        self.save_name = self.pulse_object.save_pulse(save_name=name)
    
    def update_previous_control(self):
        pass
    
    def merge_pulse_list(self):
        self.reset_pulse_object()
        if self.initial_pulse_object is not None:
            self.pulse_object.merge_pulses(self.initial_pulse_object)
        for pulse in self.pulse_list:
            self.pulse_object.merge_pulses(pulse)
    
    def get_pulse_object(self):
        if self.pulse_object is None:
            return None
        else:
            return self.pulse_object.copy_pulse()
        
    def add_pulse(self, pulse, index = None):
        if index is None:
            self.pulse_list.append(pulse)
        else:
            self.pulse_list[index] = pulse
       
    
    def add_pulse_args(self, args, index = None):
        if index is None:
            self.pulse_args.append(args)
        else:
            self.pulse_args[index] = args
        
    
    def add_pulse_pol(self, pol, index = None):
        if index is None:
            self.pulse_pol.append(pol)
        else:
            self.pulse_pol[index] = pol
        
    
    def add_pulse_power(self, power, index = None):
        if index is None:
            self.pulse_powers.append(power)
        else:
            self.pulse_powers[index] = power
    
        
    def add_pulse_phase(self, phase, index = None):
        if index is None:
            self.pulse_phases.append(phase)
        else:
            self.pulse_phases[index] = phase
    
    def add_pulse_kind(self, kind, index = None):
        if index is None:
            self.pulse_kind.append(kind)
        else:
            self.pulse_kind[index] = kind
    
    def set_central_wavelength(self, central_wavelength):
        self.central_wavelength = central_wavelength
    
    def set_t0(self, t0):
        self.t0 = t0
        
    def set_t_end(self, t_end):
        self.t_end = t_end
        
    def set_dt(self, dt):
        self.dt = dt
    
    def remove_pulse(self, index):
        self.pulse_list.pop(index)
        self.pulse_args.pop(index)
        self.pulse_pol.pop(index)
        self.pulse_powers.pop(index)
        self.pulse_phases.pop(index)
        self.pulse_kind.pop(index)
    
    def reset_pulse_list(self):
        self.pulse_list = []
        self.pulse_args = []
        self.pulse_pol = []
        self.pulse_powers = []
        self.pulse_phases = []
        self.pulse_kind = []

    def reset_pulse_object(self):
        self.pulse_object = pg.PulseGenerator(t0=self.t0, tend=self.t_end, dt=self.dt, central_wavelength=self.central_wavelength)
        
    def clear_pulse_object(self):
        self.pulse_object.clear_all()
    
    
    def add_phase(self, pulse, central_wl, phase):
        pulse.add_filter_rectangle() 
        pulse.add_phase_filter(unit = 'nm', central_f = central_wl, phase_taylor = phase)
        pulse.apply_frequency_filter()
        pulse.clear_filter()
        
        return pulse
    
    def update_control(self):
        pass
    
    ### pulse kinds ###
    
    
    def add_gaussian_time(self, args, pol = [1,0], phase_args = [0,0], power = 1, index = None):
        pulse = pg.PulseGenerator(t0=self.t0, tend=self.t_end, dt=self.dt, central_wavelength=self.central_wavelength) 
        
        central_wl = args[0]
        t0 = args[1]
        fwhm_t = args[2]
        pulse.add_gaussian_time(unit = 'nm', central_f=central_wl, t0=t0, width_t=fwhm_t, field_or_intesity='int', sig_or_fwhm='fwhm', polarisation=pol)
        
        pulse = self.add_phase(pulse, phase_args[0], phase_args[1:])
        
        pulse.set_pulse_power(power)
        
        self.add_pulse(pulse,index)
        self.add_pulse_args(args,index)
        self.add_pulse_pol(pol,index)
        self.add_pulse_power(power,index)
        self.add_pulse_phase(phase_args,index)
        #self.pulse_kind.append('gaussian_time')
        self.add_pulse_kind('gaussian_time',index)
    
    def add_gaussian_frequency(self, args, pol = [1,0], phase_args = [0,0], power = 1, index = None):
        pulse = pg.PulseGenerator(t0=self.t0, tend=self.t_end, dt=self.dt, central_wavelength=self.central_wavelength)
        
        
        central_wl = args[0]
        t0 = args[1]
        fwhm_f = args[2]
        
        pulse.add_gaussian_freq(unit='nm', central_f=central_wl, shift_time= t0, width_f=fwhm_f, field_or_intesity='int', sig_or_fwhm='fwhm', polarisation=pol)
        
        pulse = self.add_phase(pulse, phase_args[0], phase_args[1:])
        
        pulse.set_pulse_power(power)
        
        self.add_pulse(pulse,index)
        self.add_pulse_args(args,index)
        self.add_pulse_pol(pol,index)
        self.add_pulse_power(power,index)
        self.add_pulse_phase(phase_args,index)
        self.pulse_kind.append('gaussian_frequency')
        pass
    
   ### Gui ###
    
    def update_gui(self):
        self.set_t0(float(self.t0_entry.get()))
        self.set_t_end(float(self.t_end_entry.get()))
        self.set_dt(float(self.dt_entry.get()))
        self.set_central_wavelength(float(self.central_wavelength_entry.get()))
        
        self.reset_pulse_object()
        self.merge_pulse_list()
        
        self.ax_temporal_x.set_xdata(self.pulse_object.time)
        self.ax_temporal_y.set_xdata(self.pulse_object.time)
        
        self.ax_temporal_x.set_ydata(np.real(self.pulse_object.temporal_representation_x))
        self.ax_temporal_y.set_ydata(np.real(self.pulse_object.temporal_representation_y))
        
        self.ax_temporal_x_abs_pos.set_xdata(self.pulse_object.time)
        self.ax_temporal_x_abs_neg.set_xdata(self.pulse_object.time)
        
        self.ax_temporal_x_abs_pos.set_ydata(np.abs(self.pulse_object.temporal_representation_x))
        self.ax_temporal_x_abs_neg.set_ydata(-np.abs(self.pulse_object.temporal_representation_x))
        
        self.ax_temporal_y_abs_pos.set_xdata(self.pulse_object.time)
        self.ax_temporal_y_abs_neg.set_xdata(self.pulse_object.time)
        
        self.ax_temporal_y_abs_pos.set_ydata(np.abs(self.pulse_object.temporal_representation_y))
        self.ax_temporal_y_abs_neg.set_ydata(-np.abs(self.pulse_object.temporal_representation_y))
        
        self.ax_spectral_x.set_xdata(self.pulse_object.wavelengths)
        self.ax_spectral_y.set_xdata(self.pulse_object.wavelengths)
        
        self.ax_spectral_x.set_ydata(np.abs(self.pulse_object.frequency_representation_x)**2)
        self.ax_spectral_y.set_ydata(np.abs(self.pulse_object.frequency_representation_y)**2)
        
        self.ax_temporal.set_xlim([self.t0, self.t_end])
        self.ax_spectral.set_xlim([min(self.pulse_object.wavelengths), max(self.pulse_object.wavelengths)])
        
        self.ax_temporal.set_ylim(np.min([np.min(-np.abs(self.pulse_object.temporal_representation_x)),np.min(-np.abs(self.pulse_object.temporal_representation_y))]), 
                                  np.max([np.max(np.abs(self.pulse_object.temporal_representation_x)),np.max(np.abs(self.pulse_object.temporal_representation_y))]))
        
        self.ax_spectral.set_ylim(np.min([np.min(np.abs(self.pulse_object.frequency_representation_x)**2),np.min(np.abs(self.pulse_object.frequency_representation_y)**2)]), 
                                  np.max([np.max(np.abs(self.pulse_object.frequency_representation_x)**2),np.max(np.abs(self.pulse_object.frequency_representation_y)**2)]))
        
        
        self.canvas_temporal.draw()
        self.canvas_spectral.draw()
        
        # add the pulse list to the gui
        
        for i in range(len(self.pulse_list_labels)):
            self.pulse_list_labels[i].destroy()
            self.pulse_list_destroy_buttons[i].destroy()
            self.pulse_list_modify_buttons[i].destroy()
            
        
        self.pulse_list_labels = []
        self.pulse_list_destroy_buttons = []
        self.pulse_list_modify_buttons = []
        for i in range(len(self.pulse_list)):
            self.pulse_list_labels.append(tk.Label(self.gui_window, text='Pulse '+str(i)+': '+self.pulse_kind[i]))
            self.pulse_list_labels[i].grid(row=9+i, column=0)
            #tk.Label(self.gui_window, text='power: '+str(self.pulse_powers[i])).grid(row=9+i, column=1)
            
            self.pulse_list_destroy_buttons.append(tk.Button(self.gui_window, text='Remove', command=lambda i=i: [self.remove_pulse(i), self.update_gui()]))
            self.pulse_list_destroy_buttons[i].grid(row=9+i, column=1)
            
            self.pulse_list_modify_buttons.append(tk.Button(self.gui_window, text='Modify', command=lambda i=i: [self.modify_pulse_gui(i), self.update_gui()]))
            self.pulse_list_modify_buttons[i].grid(row=9+i, column=2)
        pass
    
    def add_pulse_gui(self):
        self.pulse_gui(self.add_pulse_kind_str.get())
        
    def modify_pulse_gui(self, index):
        print(self.pulse_kind)
        self.pulse_gui(self.pulse_kind[index], args = self.pulse_args[index], pol = self.pulse_pol[index], phase_args = self.pulse_phases[index], power = self.pulse_powers[index],index=index)
    
    def pulse_gui(self,pulse_kind, args = None, pol = None, phase_args = None, power = None, index = None):
        self.pulse_window = tk.Toplevel(self.gui_window)
        
        if len(phase_args) < 5: 
            for i in range(5-len(phase_args)):
                phase_args.append(0)
            #general settings
        self.label4 = tk.Label(self.pulse_window, text='Polarisation [h,v]: ')
        self.label4.grid(row=0, column=2)
        self.pol_entry_x = tk.Entry(self.pulse_window)
        if pol is None:
            self.pol_entry_x.insert(0, '1')
        else:
            self.pol_entry_x.insert(0, str(pol[0]))
        self.pol_entry_x.grid(row=0, column=3)
        self.pol_entry_y = tk.Entry(self.pulse_window)
        if pol is None:
            self.pol_entry_y.insert(0, '0')
        else:
            self.pol_entry_y.insert(0, str(pol[1]))
        self.pol_entry_y.grid(row=0, column=4)
        
        self.label5 = tk.Label(self.pulse_window, text='Phase: ')
        self.label5.grid(row=1, column=2)
        
        self.label6 = tk.Label(self.pulse_window, text='Central Wavelength phase (nm): ')
        self.label6.grid(row=2, column=2)
        self.phase_entry_wl = tk.Entry(self.pulse_window)
        if phase_args is None:
            self.phase_entry_wl.insert(0, str(self.central_wavelength))
        else:
            self.phase_entry_wl.insert(0, str(phase_args[0]))
        self.phase_entry_wl.grid(row=2, column=3)
        
        number_taylor = 4
        self.phase_entry = []
        for i in range(number_taylor):
            tk.Label(self.pulse_window, text='taylor '+str(i)+': ').grid(row=3+i, column=2)
            self.phase_entry.append(tk.Entry(self.pulse_window))
            self.phase_entry[i].grid(row=3+i, column=3)
            if phase_args is None:
                self.phase_entry[i].insert(0, '0')
            else:
                self.phase_entry[i].insert(0, str(phase_args[i+1]))
        tk.Label(self.pulse_window, text='Power / a.u.: ').grid(row=3+number_taylor, column=2)
        self.power_entry = tk.Entry(self.pulse_window)
        if power is None:
            self.power_entry.insert(0, '1')
        else:
            self.power_entry.insert(0, str(power))
        self.power_entry.grid(row=3+number_taylor, column=3)
        
        self.add_pulse_kind_button = tk.Button(self.pulse_window, text = 'Update pulse')
        self.add_pulse_kind_button.grid(row=3+number_taylor+1, column=2)
        
        # pulse specific settings (args)
        args_list_entry = []
        if pulse_kind == 'gaussian_time':
            self.pulse_window.title('Add Gaussian Time')
            self.label1 = tk.Label(self.pulse_window, text='Central Wavelength (nm): ')
            self.label1.grid(row=0, column=0)
            args_list_entry.append(tk.Entry(self.pulse_window))
            args_list_entry[0].grid(row=0, column=1)
            # if args is not None:
            #     args_list_entry[0].insert(0, str(args[0]))
            
            
            self.label2 = tk.Label(self.pulse_window, text='t0 (ps): ')
            self.label2.grid(row=1, column=0)
            args_list_entry.append(tk.Entry(self.pulse_window))
            args_list_entry[1].grid(row=1, column=1)
            # if args is not None:
            #     args_list_entry[1].insert(0, str(args[1]))
            
            self.label3 = tk.Label(self.pulse_window, text='fwhm (ps): ')
            self.label3.grid(row=2, column=0)
            args_list_entry.append(tk.Entry(self.pulse_window))
            args_list_entry[2].grid(row=2, column=1)
            # if args is not None:
            #     args_list_entry[2].insert(0, str(args[2]))
            
        
            
            self.add_pulse_kind_button.config(command = lambda: [self.add_gaussian_time([float(args_list_entry[0].get()), float(args_list_entry[1].get()), float(args_list_entry[2].get())], 
                                                                  [float(self.pol_entry_x.get()), float(self.pol_entry_y.get())],
                                                                  [float(self.phase_entry_wl.get()), float(self.phase_entry[0].get()), float(self.phase_entry[1].get()), float(self.phase_entry[2].get()), float(self.phase_entry[3].get())],
                                                                  float(self.power_entry.get()),index = index), self.update_gui(), self.pulse_window.destroy()])
            
        elif pulse_kind == 'gaussian_frequency':
            # self.pulse_window.title('Add Gaussian Frequency')
            # self.label1 = tk.Label(self.pulse_window, text='Central Wavelength (nm): ')
            # self.label1.grid(row=0, column=0)
            # args_list_entry.append(tk.Entry(self.pulse_window))
            # args_list_entry[0].grid(row=0, column=1)
            
            # self.label2 = tk.Label(self.pulse_window, text='t0 (ps): ')
            # self.label2.grid(row=1, column=0)
            # args_list_entry.append(tk.Entry(self.pulse_window))
            # args_list_entry[1].grid(row=1, column=1)
            
            # self.label3 = tk.Label(self.pulse_window, text='fwhm (nm): ')
            # self.label3.grid(row=2, column=0)
            # args_list_entry.append(tk.Entry(self.pulse_window))
            # args_list_entry[2].grid(row=2, column=1)
            
            self.label1 = tk.Label(self.pulse_window, text='Central Wavelength (nm): ')
            self.label1.grid(row=0, column=0)
            args_list_entry.append(tk.Entry(self.pulse_window))
            args_list_entry[0].grid(row=0, column=1)
            # if args is not None:
            #     args_list_entry[0].insert(0, str(args[0]))
            
            
            self.label2 = tk.Label(self.pulse_window, text='t0 (ps): ')
            self.label2.grid(row=1, column=0)
            args_list_entry.append(tk.Entry(self.pulse_window))
            args_list_entry[1].grid(row=1, column=1)
            # if args is not None:
            #     args_list_entry[1].insert(0, str(args[1]))
            
            self.label3 = tk.Label(self.pulse_window, text='fwhm (nm): ')
            self.label3.grid(row=2, column=0)
            args_list_entry.append(tk.Entry(self.pulse_window))
            args_list_entry[2].grid(row=2, column=1)
            # if args is not None:
            #     args_list_entry[2].insert(0, str(args[2]))
            
        
            
            self.add_pulse_kind_button.config(command = lambda: [self.add_gaussian_frequency([float(args_list_entry[0].get()), float(args_list_entry[1].get()), float(args_list_entry[2].get())], 
                                                                  [float(self.pol_entry_x.get()), float(self.pol_entry_y.get())],
                                                                  [float(self.phase_entry_wl.get()), float(self.phase_entry[0].get()), float(self.phase_entry[1].get()), float(self.phase_entry[2].get()), float(self.phase_entry[3].get())],
                                                                  float(self.power_entry.get()),index = index), self.update_gui(), self.pulse_window.destroy()])
            
            
            
            
        
        if args is not None:
            for i in range(len(args)):
                args_list_entry[i].insert(0, str(args[i]))
    
        
        
            
            
            
    
    def gui(self):
        self.open_gui = True
        
        self.gui_window = tk.Tk()
        def close_gui():
            self.gui_window.destroy()
            print('Pulse Generator GUI closed')
        
        self.gui_window.protocol("WM_DELETE_WINDOW", close_gui)
        self.gui_window.title('Pulse Generator')
        
        fig_temporal = Figure(figsize=(5, 4), dpi=100)
        fig_spectral = Figure(figsize=(5, 4), dpi=100)
        
        self.ax_temporal = fig_temporal.add_subplot()
        self.ax_spectral = fig_spectral.add_subplot()
        
        self.ax_temporal.set_title('Temporal')
        self.ax_spectral.set_title('Spectral')
        
        self.ax_temporal.set_xlabel('Time (ps)')
        self.ax_temporal.set_ylabel('field (a.u.)')
        
        self.ax_temporal.set_xlim([self.t0, self.t_end])
        
        self.ax_spectral.set_xlabel('Wavelength (nm)')
        self.ax_spectral.set_ylabel('Intensity (a.u.)')
        
        self.ax_spectral.set_xlim([min(self.pulse_object.wavelengths), max(self.pulse_object.wavelengths)])
        
        self.ax_temporal_x, = self.ax_temporal.plot(self.pulse_object.time, np.real(self.pulse_object.temporal_representation_x), 'b-', alpha=0.5)
        self.ax_temporal_y, = self.ax_temporal.plot(self.pulse_object.time, np.real(self.pulse_object.temporal_representation_y), 'r-', alpha=0.5)
        
        self.ax_temporal_x_abs_pos, = self.ax_temporal.plot(self.pulse_object.time, np.abs(self.pulse_object.temporal_representation_x), 'b-')
        self.ax_temporal_x_abs_neg, = self.ax_temporal.plot(self.pulse_object.time, -np.abs(self.pulse_object.temporal_representation_x), 'b-')
        
        self.ax_temporal_y_abs_pos, = self.ax_temporal.plot(self.pulse_object.time, np.abs(self.pulse_object.temporal_representation_y), 'r-')
        self.ax_temporal_y_abs_neg, = self.ax_temporal.plot(self.pulse_object.time, -np.abs(self.pulse_object.temporal_representation_y), 'r-')
        
        self.ax_spectral_x, = self.ax_spectral.plot(self.pulse_object.wavelengths, np.abs(self.pulse_object.frequency_representation_x)**2, 'b-')
        self.ax_spectral_y, = self.ax_spectral.plot(self.pulse_object.wavelengths, np.abs(self.pulse_object.frequency_representation_y)**2, 'r-')
        
        self.canvas_temporal = FigureCanvasTkAgg(fig_temporal, master=self.gui_window)
        self.canvas_temporal.draw()
        
        self.canvas_spectral = FigureCanvasTkAgg(fig_spectral, master=self.gui_window)
        self.canvas_spectral.draw()
        
        self.canvas_temporal.get_tk_widget().grid(row=0, column=0,columnspan=1,rowspan=7)
        self.canvas_spectral.get_tk_widget().grid(row=0, column=1,columnspan=1,rowspan=7)
        
        tk.Label(self.gui_window, text='Pulse settings').grid(row=0, column=2)
        
        self.update_button = tk.Button(self.gui_window, text='Update', command=self.update_gui)
        self.update_button.grid(row=0, column=3)
        
        self.clear_button = tk.Button(self.gui_window, text='Clear', command= lambda: [self.clear_pulse_object(),self.reset_pulse_list(),self.update_gui()])
        self.clear_button.grid(row=0, column=4)
        
        tk.Label(self.gui_window, text='start: ').grid(row=1, column=2)
        self.t0_entry = tk.Entry(self.gui_window)
        self.t0_entry.insert(0, str(self.t0))
        self.t0_entry.grid(row=1, column=3)
        tk.Label(self.gui_window, text='ps').grid(row=1, column=4)
        
        tk.Label(self.gui_window, text='end: ').grid(row=2, column=2)
        self.t_end_entry = tk.Entry(self.gui_window)
        self.t_end_entry.insert(0, str(self.t_end))
        self.t_end_entry.grid(row=2, column=3)
        tk.Label(self.gui_window, text='ps').grid(row=2, column=4)
        
        tk.Label(self.gui_window, text='dt: ').grid(row=3, column=2)
        self.dt_entry = tk.Entry(self.gui_window)
        self.dt_entry.insert(0, str(self.dt))
        self.dt_entry.grid(row=3, column=3)
        tk.Label(self.gui_window, text='ps').grid(row=3, column=4)
        
        tk.Label(self.gui_window, text='central wavelength: ').grid(row=4, column=2)
        self.central_wavelength_entry = tk.Entry(self.gui_window)
        self.central_wavelength_entry.insert(0, str(self.central_wavelength))
        self.central_wavelength_entry.grid(row=4, column=3)
        tk.Label(self.gui_window, text='nm').grid(row=4, column=4)
        
        self.add_pulse_button = tk.Button(self.gui_window, text='Add pulse', command=self.add_pulse_gui)
        self.add_pulse_button.grid(row=8, column=0)
        
        self.add_pulse_kind_str = tk.StringVar()
        self.drop_down = ttk.Combobox(self.gui_window, values=self.pulse_kinds, textvariable=self.add_pulse_kind_str, width=50, state='readonly')
        self.drop_down.current(0)
        self.drop_down.grid(row=8, column=1)
        
        self.pulse_list_labels = []
        self.pulse_list_destroy_buttons = []
        self.pulse_list_modify_buttons = []
        self.update_gui()
        
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
    
    qd_calibration = 'QD_Iker_April_high_fss.txt' 
    initial_pulse = pg.PulseGenerator(0,50,0.2,calibration_file=qd_calibration)
    initial_pulse.add_gaussian_time(unit = 'nm', central_f= 779.89, width_t=0.1,t0=50/2,area_time=20, polarisation=[1,0.6]) # <-- more power for tpe 
    
    ic = Pulse_Generator(open_gui=False, t_end=250,initial_pulse=None)
    
    ic.add_gaussian_time([802, 30, 10],pol=[0,1],phase_args=[800.5,0,0,0],power=50)
    
    #ic.add_gaussian_time([800, 20, 10],pol=[1,0],phase_args=[778,0,0,0],power=100)
    
    #ic.merge_pulse_list()
    ic.gui()
    ic.start_gui()