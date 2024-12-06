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



class pulse_shaper_control(): 
    def __init__(self,pulse_shaper_object,pulse_object = None, step_size = 0.1, large_step = 10, initial_position = None, open_gui = True, parent_window = None, previous_control = None)->None:
        
        self.pulse_shaper_object = pulse_shaper_object
        self.name = pulse_shaper_object.name
        self.step_size = step_size
        self.large_step = large_step
        self.open_gui = open_gui
        self.pulse_object_out = None
        self.previous_control = previous_control
        
        if pulse_object is not None:
            if type(pulse_object) is str:
                 pulse_object = pg.load_pulse(pulse_object)
            self.pulse_object = pulse_object.copy_pulse()
            self.pulse_object.clear_filter()
            
        else:
            self.pulse_object = None
        
        if initial_position is None:
            self.current_position = np.mean([self.pulse_shaper_object.wavelength_min,self.pulse_shaper_object.wavelength_max])
            
        if self.pulse_shaper_object.device is None:
            self.excecute = False
        else:
            self.excecute = True
        
        
        # initiallize the buttons
        self.button_step_increase = 0
        self.button_step_decrease = 0
        
        self.button_step_small_up = 0
        self.button_step_small_down = 0
        
        self.button_step_large_up = 0
        self.button_step_large_down = 0 
        
        
        self.update_previous_control()
        self.update_shaper()
        self.print_info()
        
        if self.open_gui:
            #print('hier könnte ihre GUI stehen!')
            return self.gui(parentWindow=parent_window)
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
        
    def update_shaper(self):
        #self.update_previous_control()
        self.rescale_step_size()
        self.current_position = self.current_position + self.step_size*self.button_step_small_up - self.step_size*self.button_step_small_down + self.step_size*self.button_step_large_up*self.large_step - self.step_size*self.button_step_large_down*self.large_step

        self.pulse_object_out = self.pulse_shaper_object.move_slit(self.current_position,pulse_object = self.pulse_object,excecute = self.excecute)
        
        self.mask_wavelengths = self.pulse_shaper_object.mask_wavelengths
        self.mask = self.pulse_shaper_object.mask
        # self.reset_buttons()
        
    # def button_press(self,button):
    #     if button == 0:
    #         button = 1 
    #         print('Button pressed: '+ str(button))
    #     else:
    #         print('Button not pressed: '+ str(button))
    #     self.get_button_state()
    
    # button presses
    
    def update_control(self):
        self.update_shaper()
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
        
    def set_current_position(self,position):
        self.current_position = position
        
    def set_control_value(self,control_value):
        self.current_position = control_value
    
    def get_control_value(self):
        return self.current_position
    
    def get_pulse_object(self):
        if self.pulse_object_out is None:
            return None
        return self.pulse_object_out.copy_pulse()
    
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
        print('Pulseshaper: ' + self.pulse_shaper_object.name)
        print('Current position: '+ str(self.pulse_shaper_object.slit_center_wavelength))
        print('Current step size: '+ str(self.step_size))
        
    def get_button_state(self):
        print('PS buttons: ')
        print('Step increase: '+ str(self.button_step_increase))
        print('Step decrease: '+ str(self.button_step_decrease))
        print('Step small up: '+ str(self.button_step_small_up))
        print('Step small down: '+ str(self.button_step_small_down))
        print('Step large up: '+ str(self.button_step_large_up))
        print('Step large down: '+ str(self.button_step_large_down))
        
    def close(self):
        self.pulse_shaper_object.close()
        
    def update_gui(self):
        #self.button_press(button)
        self.set_step_size(float(self.step_size_entry.get()))
        self.update_shaper()
        #self.print_info()
        #self.get_button_state()
        self.reset_buttons()
        self.mask_plot.set_xdata(self.mask_wavelengths)
        self.mask_plot.set_ydata(self.mask)
        self.ax_mask.set_xlim([self.mask_wavelengths[0],self.mask_wavelengths[-1]])
        self.ax_mask.set_ylim([0, np.max(np.ceil(self.mask*10)/10)])
        self.current_position_entry.delete(0, 'end')
        self.current_position_entry.insert(0,str(np.round(self.current_position,decimals=3)))
        
        self.step_size_entry.delete(0, 'end')
        self.step_size_entry.insert(0,str(np.round(self.step_size,decimals=3)))
        
        self.canvas_fig_mask.draw()
        
    def gui(self, parentWindow = None): 
        self.open_gui = True
        if self.previous_control is None:
            self.gui_window = tk.Tk()
        else:
            self.gui_window = tk.Toplevel(self.previous_control.gui_window)
            
    
        
        
        def close_gui():
            self.gui_window.destroy()
            if self.pulse_shaper_object.device is not None:
                self.pulse_shaper_object.close()
            
            print(self.name + ' closed')
        
        def move_slit():
            self.set_current_position(float(self.current_position_entry.get()))
            self.update_gui()

        # # update the gui after 100 ms to catch the
        # button presses 
            
        self.gui_window.protocol("WM_DELETE_WINDOW", close_gui)
        #self.gui_window.geometry('500x400')
        self.gui_window.title('Pulse shaper control: '+ self.pulse_shaper_object.name)
       
        
        # mask plot
        fig_mask = Figure(figsize=(3, 2), dpi=100)
        self.ax_mask = fig_mask.add_subplot()
        self.mask_plot, = self.ax_mask.plot(self.mask_wavelengths,self.mask/np.max(self.mask),'r-')
        self.ax_mask.set_xlabel('Wavelength (nm)')        
        
        self.canvas_fig_mask = FigureCanvasTkAgg(fig_mask, master=self.gui_window)
        self.canvas_fig_mask.draw()
        
        self.canvas_fig_mask.get_tk_widget().grid(row=0,column=0, columnspan=4)
        # bottons 
        tk.Button(self.gui_window, text = '<<', width= 10, command=lambda: [self.set_button_step_large_down(), self.update_gui()]).grid(row=1,column=0)
        
        tk.Button(self.gui_window, text = '<', width= 10, command=lambda: [self.set_button_step_small_down(), self.update_gui()]).grid(row=1,column=1)
        
        tk.Button(self.gui_window, text = '>', width= 10, command=lambda: [self.set_button_step_small_up(), self.update_gui()]).grid(row=1,column=2)
        
        tk.Button(self.gui_window, text = '>>', width= 10, command=lambda: [self.set_button_step_large_up(), self.update_gui()]).grid(row=1,column=3)
        
        tk.Button(self.gui_window, text = 'Move slit:', command= move_slit).grid(row=2,column=0)
        
        self.current_position_entry = tk.Entry(master= self.gui_window, width=10)
        self.current_position_entry.insert(0,str(np.round(self.current_position,decimals=3)))
        self.current_position_entry.grid(row=2,column=1)
        
        tk.Button(self.gui_window, text='step +', width= 10, command=lambda: [self.set_button_step_increase(), self.update_gui()]).grid(row=2,column=2)
        
        tk.Button(self.gui_window, text='step -', width= 10, command=lambda: [self.set_button_step_decrease(), self.update_gui()]).grid(row=3,column=2)
        
        self.step_size_entry = tk.Entry(master= self.gui_window, width=10)
        self.step_size_entry.insert(0,str(np.round(self.step_size,decimals=3)))
        self.step_size_entry.grid(row=2,column=3)
        
        tk.Button(self.gui_window, text = 'Calibration file', command=self.calibration_file_gui).grid(row=3,column=0)
        
        self.update_gui()
        #self.gui_window.after(500, self.update_gui)
        
    def calibration_file_gui(self):
        config = configparser.ConfigParser()
        config.read(self.pulse_shaper_object.calibration_file)
        
        poly_0 = config['pulse_shaper']['poly_0']   
        poly_1 = config['pulse_shaper']['poly_1']
        slit_width = config['pulse_shaper']['slit_width']
        psf_width = config['pulse_shaper']['psf_width']
        slit_pos_min = config['pulse_shaper']['slit_pos_min']
        slit_pos_max = config['pulse_shaper']['slit_pos_max']
        
        self.calibration_file_window = tk.Toplevel(self.gui_window)
        self.calibration_file_window.title('Edit calibration file: '+ self.name)
        
        tk.Label(self.calibration_file_window, text='File name: ').grid(row=0,column=0)
        
        self.calibration_file_name_entry = tk.Entry(master= self.calibration_file_window, width=20)
        self.calibration_file_name_entry.insert(0,self.pulse_shaper_object.calibration_file)
        self.calibration_file_name_entry.grid(row=0,column=1)
        
        tk.Label(self.calibration_file_window, text = '------').grid(row=1,column=0)
        
        tk.Label(self.calibration_file_window, text = 'Slit width (nm): ').grid(row=2,column=0)
        
        self.slit_width_entry = tk.Entry(master= self.calibration_file_window, width=10)
        self.slit_width_entry.insert(0,slit_width)
        self.slit_width_entry.grid(row=2,column=1)
        
        tk.Label(self.calibration_file_window, text = 'PSF width (nm): ').grid(row=3,column=0)
        
        self.psf_width_entry = tk.Entry(master= self.calibration_file_window, width=10)
        self.psf_width_entry.insert(0,psf_width)
        self.psf_width_entry.grid(row=3,column=1)
        
        tk.Label(self.calibration_file_window, text = '------').grid(row=4,column=0)
        
        tk.Label(self.calibration_file_window, text = 'Polynom 0 (nm): ').grid(row=5,column=0)
        
        self.poly_0_entry = tk.Entry(master= self.calibration_file_window, width=10)
        self.poly_0_entry.insert(0,poly_0)
        self.poly_0_entry.grid(row=5,column=1)
        
        tk.Label(self.calibration_file_window, text = 'Polynom 1 (nm/motor unit): ').grid(row=6,column=0)
        
        self.poly_1_entry = tk.Entry(master= self.calibration_file_window, width=10)
        self.poly_1_entry.insert(0,poly_1)
        self.poly_1_entry.grid(row=6,column=1)
        
        tk.Label(self.calibration_file_window, text = '------').grid(row=7,column=0)
        
        tk.Label(self.calibration_file_window, text = 'calibration limit 1: ').grid(row=8,column=0)
        
        self.slit_pos_min_entry = tk.Entry(master= self.calibration_file_window, width=10)
        self.slit_pos_min_entry.insert(0,slit_pos_min)
        self.slit_pos_min_entry.grid(row=8,column=1)
        
        tk.Label(self.calibration_file_window, text = 'calibration limit 2: ').grid(row=9,column=0)
        
        self.slit_pos_max_entry = tk.Entry(master= self.calibration_file_window, width=10)
        self.slit_pos_max_entry.insert(0,slit_pos_max)
        self.slit_pos_max_entry.grid(row=9,column=1)
        
        
        self.load_calibration_button = tk.Button(self.calibration_file_window, text = 'Load')
        self.load_calibration_button.config(command=lambda: [self.load_calibration_file(), self.update_gui()])
        self.load_calibration_button.grid(row=10,column=0)
        
        self.save_calibration_file_button = tk.Button(self.calibration_file_window, text = 'Save')
        self.save_calibration_file_button.config(command=lambda: [self.save_calibration_file(), self.update_gui()])
        self.save_calibration_file_button.grid(row=10,column=1)
        
        
        
    
    def load_calibration_file(self):
        filepath = tk.filedialog.askopenfilename()
        # extract only the filename
        filepath = filepath.split('/')[-1]
        self.pulse_shaper_object.set_calibration_file(filepath)
        self.calibration_file_name_entry.delete(0, 'end')
        self.calibration_file_name_entry.insert(0,filepath)
        
        self.calibration_file_window.destroy()
        self.calibration_file_gui()
        
        
        
    def save_calibration_file(self):
        filepath = os.getcwd() + '/' + self.calibration_file_name_entry.get()
        config = configparser.ConfigParser()
        config['pulse_shaper'] = {'poly_0': self.poly_0_entry.get(),
                                   'poly_1': self.poly_1_entry.get(),
                                   'slit_width': self.slit_width_entry.get(),
                                   'psf_width': self.psf_width_entry.get(),
                                   'slit_pos_min': self.slit_pos_min_entry.get(),
                                   'slit_pos_max': self.slit_pos_max_entry.get()}
        
        with open(filepath, 'w') as configfile:
            config.write(configfile)
            # close the file
            configfile.close()
        
        self.pulse_shaper_object.set_calibration_file(self.calibration_file_name_entry.get())
        
    def start_gui(self):
        self.gui_window.mainloop()
    
     
     
     #########################################   
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
    ps_calibration = 'calibration_Cam_230424.txt' 
    ps_calibration2 = 'calibration_slit_38.txt' 
    qd_calibration = 'QD_Iker_April.txt'

    lab_motor = motor(fake_motor(),name='motor1')
    lab_motor2 = motor(fake_motor(),name='motor2')

    t_0 = 0 
    t_end = 100

    initial_pulse = pg.PulseGenerator(t_0,t_end,0.2,calibration_file=qd_calibration)
    initial_pulse.add_gaussian_time(width_t=0.5,t0=t_end/2)
    initial_pulse.plot_pulses(domain='nm')

    pulse_shaper = pulse_shaper_obj(device=lab_motor, calibration_file=ps_calibration, name='heisl')
    
    
    ps_controller = pulse_shaper.open_control(pulse_object=initial_pulse)
    
    pulse_shaper2 = pulse_shaper_obj(device=lab_motor2, calibration_file=ps_calibration2, name='bert')
    ps_controller2 = pulse_shaper2.open_control(pulse_object = ps_controller.get_pulse_object())
    
    ps_controller.start_gui()

    pass
    
        
        

    
   