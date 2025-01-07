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
import sys

# ----- Collaps of the GUI -----

class gui_manager():
    def __init__(self, device_control = [], inital_state_all = False, open_gui = True, num_coloums = 1) -> None:
        self.device_control = device_control
        self.inital_state_all = inital_state_all
        
        self.gui_window_name = []
        self.gui_window_list = []
        self.gui_window_state = []
        
        self.num_coloums = num_coloums
        
        self.open_gui = open_gui
        
        
        if self.open_gui:
            self.gui()
        
    def gui(self):
        self.initialise()
        self.open_gui = True
        self.gui_window = tk.Tk()
        
        self.gui_window.title('GUI manager')
        
        # when closing the window it closes all connected control windows 
        def close_windows():
            for i in self.device_control:
                # only destroy if the window is not already destroyed
                if hasattr(i,'close'):
                    i.close()
                i.gui_window.destroy
            self.gui_window.destroy()
            sys.exit()
            
        self.gui_window.protocol("WM_DELETE_WINDOW", close_windows)
        
        self.num_rows = int(np.ceil(self.num_devices/self.num_coloums))
        
        self.name_label_list = [] 
        self.toggle_button_list = []
        for gui_index in range(self.num_devices):
            self.name_label_list.append(tk.Label(self.gui_window, text = self.gui_window_name[gui_index]))
            self.name_label_list[gui_index].grid(row = int(gui_index/self.num_coloums), column = gui_index%self.num_coloums*2)
            
            self.toggle_button_list.append(tk.Button(self.gui_window, text = 'Toggle', command = lambda gui_index = gui_index: self.toggle_gui(gui_index)))
            self.toggle_button_list[gui_index].grid(row = int(gui_index/self.num_coloums), column = gui_index%self.num_coloums*2+1)
            
            button_color(self.toggle_button_list[gui_index], self.gui_window_state[gui_index])
            
        
        pass
    
    def add_controller(self, controller = []):
        for i in controller:
            self.device_control.append(i)
            
    def toggle_gui(self, controller_index):
        if self.gui_window_state[controller_index]:
            self.collapse_gui(controller_index)
        else:
            self.expand_gui(controller_index)
        
        button_color(self.toggle_button_list[controller_index], self.gui_window_state[controller_index])
        
    def set_initial_state_all(self, state = False):
        self.inital_state_all = state
        
    def collapse_gui(self, controller_index = []):
        self.gui_window_list[controller_index].withdraw()
        self.gui_window_state[controller_index] = False
    
    def expand_gui(self, controller_index = []):
        
        self.gui_window_list[controller_index].deiconify()
        self.gui_window_list[controller_index].focus_force()
        self.gui_window_state[controller_index] = True

    def initialise(self):
        for index, controller in enumerate(self.device_control):
            self.gui_window_state.append(self.inital_state_all)
            self.gui_window_list.append(controller.gui_window)
            self.gui_window_name.append(controller.name)
            if self.inital_state_all:
               self.expand_gui(index)
            else:
                self.collapse_gui(index)
            
            #device.gui_window.protocol("WM_DELETE_WINDOW",lambda: [device.gui_window.withdraw(), button_color(self.toggle_button_list[i], self.gui_window_state[i])])
            controller.gui_window.protocol("WM_DELETE_WINDOW", lambda controller_index = index: self.toggle_gui(controller_index))
        self.num_devices = len(self.device_control)
            
    
    def start_gui(self):
        self.gui_window.mainloop()
        
        
def button_color(button_object, state):
            if state:
                button_object.config(bg='green')
            else:
                button_object.config(bg='white')