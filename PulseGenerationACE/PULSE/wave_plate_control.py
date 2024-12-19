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

class waveplate_control:
    def __init__(self,rotor_object, pulse_object = None, step_size = 1, large_step = 10, initial_angle = None, open_gui = True, parent_window = None, previous_control = None):
        
        self.rotor_object = rotor_object
        self.name = rotor_object.name
        if rotor_object.get_unit() == 'deg':
            self.step_size = step_size
        else:
            self.step_size = step_size*np.pi/180
        self.large_step = large_step
        self.initial_angle = initial_angle
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
    
        if initial_angle is None:
            self.current_angle = 0
            
        if self.rotor_object.device is None:
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
        self.update_angle()
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
    
    def update_angle(self):
        #self.update_previous_control()
        self.rescale_step_size()
        self.current_angle = self.current_angle + self.step_size*self.button_step_small_up - self.step_size*self.button_step_small_down + self.step_size*self.button_step_large_up*self.large_step - self.step_size*self.button_step_large_down*self.large_step
        
        # modulo 360 if deg and 2*pi if rad
        if self.rotor_object.get_unit() == 'deg':
            self.current_angle = self.current_angle % 360
        else:
            self.current_angle = self.current_angle % (2*np.pi)
        
        self.pulse_object_out = self.rotor_object.rotate(angle = self.current_angle, pulse_object = self.pulse_object,excecute = self.excecute)
    
    def update_control(self):
        self.update_angle()
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
        
    def set_current_angle(self,angle):
        self.current_angle = angle
        
    def set_control_value(self,control_value):
        self.current_angle = control_value
    
    def get_control_value(self):
        return self.current_angle
    
    def get_current_angle(self):
        return self.current_angle
    
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
        print('Current angle: ', self.current_angle)
        print('Step size: ', self.step_size)
        print('Large step size: ', self.large_step)
        print('Initial angle: ', self.initial_angle)
        print('Open GUI: ', self.open_gui)
        print('Parent window: ', self.parent_window)
        print('Previous control: ', self.previous_control)
        
    def get_button_state(self):
        print('Wave plate buttons')
        print('Step increase: ', self.button_step_increase)
        print('Step decrease: ', self.button_step_decrease)
        print('Small step up: ', self.button_step_small_up)
        print('Small step down: ', self.button_step_small_down)
        print('Large step up: ', self.button_step_large_up)
    
    def close(self):
        self.rotor_object.close()
    
    def update_gui(self):
        #self.button_press(button)
        self.set_step_size(float(self.step_size_entry.get()))
        self.update_angle()
        #self.print_info()
        #self.get_button_state()
        self.reset_buttons()
        self.current_angle_entry.delete(0, 'end')
        self.current_angle_entry.insert(0,str(np.round(self.current_angle,decimals=3)))
        
        self.step_size_entry.delete(0, 'end')
        self.step_size_entry.insert(0,str(np.round(self.step_size,decimals=3)))
        
        self.angle_bar['value'] = self.current_angle
        
        
    def gui(self):
        self.open_gui = True
        if self.previous_control is None:
            self.gui_window = tk.Tk()
        else:
            self.gui_window = tk.Toplevel(self.previous_control.gui_window)
            
        def close_gui():
            self.gui_window.destroy()
            if self.rotor_object.device is not None:
                self.rotor_object.close()
            
            print('Window closed')
        
        def apply_angle():
            self.set_current_angle(float(self.current_angle_entry.get()))
            self.update_gui()

        # # update the gui after 100 ms to catch the
        # button pres
            
        self.gui_window.protocol("WM_DELETE_WINDOW", close_gui)
        self.gui_window.geometry('500x100')
        self.gui_window.title('Waveplate control: '+ self.rotor_object.name)
        
        tk.Button(self.gui_window, text = '<<', width= 10, command=lambda: [self.set_button_step_large_down(), self.update_gui()]).grid(row=1,column=0)
        
        tk.Button(self.gui_window, text = '<', width= 10, command=lambda: [self.set_button_step_small_down(), self.update_gui()]).grid(row=1,column=1)
        
        tk.Button(self.gui_window, text = '>', width= 10, command=lambda: [self.set_button_step_small_up(), self.update_gui()]).grid(row=1,column=2)
        
        tk.Button(self.gui_window, text = '>>', width= 10, command=lambda: [self.set_button_step_large_up(), self.update_gui()]).grid(row=1,column=3)
        
        tk.Button(self.gui_window, text = 'Set angle:', command= apply_angle).grid(row=2,column=0)
        
        self.current_angle_entry = tk.Entry(master= self.gui_window, width=10)
        self.current_angle_entry.insert(0,str(np.round(self.current_angle,decimals=3)))
        self.current_angle_entry.grid(row=2,column=1)
        
        tk.Button(self.gui_window, text='step +', width= 10, command=lambda: [self.set_button_step_increase(), self.update_gui()]).grid(row=2,column=2)
        
        tk.Button(self.gui_window, text='step -', width= 10, command=lambda: [self.set_button_step_decrease(), self.update_gui()]).grid(row=3,column=2)
        
        self.step_size_entry = tk.Entry(master= self.gui_window, width=10)
        self.step_size_entry.insert(0,str(np.round(self.step_size,decimals=3)))
        self.step_size_entry.grid(row=2,column=3)
        
        if self.rotor_object.get_unit() == 'deg':
            max_angle = 360
        else:
            max_angle = 2*np.pi
        
        self.angle_bar = ttk.Progressbar(master=self.gui_window, orient='horizontal', length=200, mode='determinate', maximum=max_angle, value=self.current_angle)
        self.angle_bar.grid(row=3,column=0)
        
    def start_gui(self):
        self.gui_window.mainloop()
        
if __name__ == '__main__':

    
    
    
    

    pass
        