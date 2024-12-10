import numpy as np
from matplotlib import pyplot as plt
from scipy.io import savemat
import pyaceqd.pulsegenerator as pg 

import tkinter as tk
from tkinter import ttk
from matplotlib.backend_bases import key_press_handler
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,
                                               NavigationToolbar2Tk)
from matplotlib.figure import Figure

from scipy.signal import savgol_filter
  
from Pulse_v2 import load_pulse_device, pulse_shaper_obj, motor, fake_motor, spectrometer, fake_spectrometer
import datetime
import os
### 

#### connect your devices (as you would usually do) -> Calibration needs one motor and one spectrometer 
lab_motor = fake_motor() 
lab_spectrometer = fake_spectrometer() 

#### transform your devices into PULSE devices
global MOTOR 
global SPEC 

MOTOR = motor(lab_motor, name= 'Fancy motor',motor_sleep = 0.01)
SPEC = spectrometer(lab_spectrometer, name= 'Fancy spectrometer')
#### 


#set your working directory
cwd = os.getcwd()
print(cwd)
os.chdir(cwd)

def close_function():
        MOTOR.close()
        SPEC.close()
        root.destroy()

#initialize window 
root = tk.Tk()
root.protocol("WM_DELETE_WINDOW", close_function)
root.resizable(True, True)
root.wm_title("Pulse shaper calibration")

#initialize figure 
wavelength_lab = np.linspace(790, 800, 10)
counts_envelope = np.zeros_like(wavelength_lab)
fig = Figure(figsize=(5, 4), dpi=100)
ax = fig.add_subplot() 
env, = ax.plot(wavelength_lab, counts_envelope,'k-')
slice, = ax.plot(wavelength_lab, counts_envelope,'g-')
center_slice, = ax.plot(wavelength_lab, counts_envelope,'r-')
ax.set_xlabel('Wavelength (nm)')
ax.set_ylabel('Counts')

canvas = FigureCanvasTkAgg(fig, master=root)  # A tk.DrawingArea.
canvas.draw()

# def connect_spectrometer():
#     global SPEC
#     filename = tk.filedialog.askopenfilename()
#     SPEC = load_pulse_device(filename)
#     if SPEC.suffix != 'sp':
#         print('Device is not a spectrometer!')
#         return
#     connect_spectrometer.configure(text="Spectrometer connected: "+SPEC.name)

# def connect_motor():
#     global MOTOR
#     filename = tk.filedialog.askopenfilename()
#     MOTOR = load_pulse_device(filename)
#     if MOTOR.suffix != 'mo':
#         print('Device is not a motor!')
#         return
#     connect_motor.configure(text="Motor connected: "+MOTOR.name)

def aquire_data():
    # if SPEC not in globals():
    #     print('Spectrometer not connected!')
    #     return
    # if not motor_connected:
    #     print('Motor not connected!')
    #     return
    #SPEC.get_spectrum()
    #MOTOR.set_position(0.5,excecute = True)
    aquire_window = tk.Toplevel(root)
    aquire_window.title('Aquire settings')
    tk.Label(aquire_window, text="Motor min").grid(row=0)
    tk.Label(aquire_window, text="Motor max").grid(row=1)
    tk.Label(aquire_window, text="Motor step").grid(row=2)
    motor_min = tk.Entry(aquire_window)
    motor_max = tk.Entry(aquire_window)
    motor_step = tk.Entry(aquire_window)

    motor_min.insert(0, '0')
    motor_max.insert(0, '5')
    motor_step.insert(0, '0.1')

    motor_min.grid(row=0, column=1)
    motor_max.grid(row=1, column=1)
    motor_step.grid(row=2, column=1)

    def load_folder():
        
        old_dir = os.getcwd()
        enable_gui()
        folder = tk.filedialog.askdirectory()
        os.chdir(folder)
        files = os.listdir(folder)
        files.sort(key=lambda x: os.path.getmtime(x))
        global MO_POS 
        global WL_MAT
        global COUNTS_MAT
        global COUNTS_ENV
        MO_POS = np.linspace(float(motor_min.get()), float(motor_max.get()), len(files))
        for i, file in enumerate(files):
            data = np.loadtxt(file)
            if i == 0:
                WL_MAT = np.zeros((len(files), len(data)))
                COUNTS_MAT = np.zeros((len(files), len(data)))
            WL_MAT[i,:] = data[:,0]
            COUNTS_MAT[i,:] = data[:,1]
        COUNTS_ENV = np.amax(COUNTS_MAT, axis=0)
        env.set_data(WL_MAT[0,:], COUNTS_ENV)
        slice.set_data(WL_MAT[0,:], COUNTS_MAT[0,:])
        center_slice.set_data(WL_MAT[0,:], COUNTS_MAT[0,:])
        ax.set_xlim([np.min(WL_MAT[0,:]), np.max(WL_MAT[0,:])])
        ax.set_ylim([0, np.max(COUNTS_ENV)*1.1])
        canvas.draw()
        slider_slice.config(from_=float(motor_min.get()), to=float(motor_max.get()), resolution=np.abs(MO_POS[1]-MO_POS[0]))
        slider_slice.set(np.mean([float(motor_min.get()),float(motor_max.get())]))
        
        os.chdir(old_dir)
        xlim_min.delete(0, tk.END)
        xlim_max.delete(0, tk.END)
        xlim_min.insert(0, str(np.round(np.min(WL_MAT[0,:]),decimals=2)))
        xlim_max.insert(0, str(np.round(np.max(WL_MAT[0,:]),decimals=2)))
        num_slices.config(from_=2, to = len(MO_POS), resolution=1)
        separation_slices.delete(0, tk.END)
        separation_slices.insert(0,string=str(np.round(6*np.abs(MO_POS[1]-MO_POS[0]),decimals=2)))
        
        aquire_window.destroy()

    def _aquire():
        enable_gui()
        min = float(motor_min.get())
        max = float(motor_max.get())
        step = float(motor_step.get())

        global MO_POS 
        global WL_MAT
        global COUNTS_MAT
        global COUNTS_ENV
        MO_POS = np.arange(min, max+step, step, dtype=float,)
        
        for i, pos in enumerate(MO_POS):
            MOTOR.set_position(pos, excecute = True)
            spec_out = SPEC.get_spectrum(excecute = True)
            wavelength, counts = spec_out[0]
            print(str(i+1) + '/' + str(len(MO_POS)) + ' done')
            if i == 0:
                WL_MAT = np.zeros((len(MO_POS), len(wavelength)))
                COUNTS_MAT = np.zeros((len(MO_POS), len(counts)))
            WL_MAT[i,:] = wavelength
            COUNTS_MAT[i,:] = counts

        COUNTS_ENV = np.amax(COUNTS_MAT, axis=0)
        env.set_data(wavelength, COUNTS_ENV)
        slice.set_data(WL_MAT[0,:], COUNTS_MAT[0,:])
        center_slice.set_data(WL_MAT[0,:], COUNTS_MAT[0,:])
        ax.set_xlim([np.min(wavelength), np.max(wavelength)])
        ax.set_ylim([0, np.max(COUNTS_ENV)*1.1])
        canvas.draw()

        slider_slice.config(from_=min, to=max, resolution=step)
        slider_slice.set(np.mean([min,max]))
        xlim_min.delete(0, tk.END)
        xlim_max.delete(0, tk.END)
        xlim_min.insert(0, str(np.round(np.min(WL_MAT[0,:]),decimals=2)))
        xlim_max.insert(0, str(np.round(np.max(WL_MAT[0,:]),decimals=2)))
        num_slices.config(from_=2, to = len(MO_POS), resolution=1)
        separation_slices.delete(0, tk.END)
        separation_slices.insert(0,string=str(np.round(6*step,decimals=2)))
        
        aquire_window.destroy()

    
    tk.Button(aquire_window, text='Aquire', command=_aquire).grid(row=3, column=0, sticky=tk.W, pady=4)
    tk.Button(aquire_window, text='Load folder', command=load_folder).grid(row=3, column=1, sticky=tk.W, pady=4)
def draw_slice(val):
    pos = slider_slice.get()
    global POS_VEC
    global IDX
    POS_VEC = [pos]
    IDX = [np.argmin(np.abs(MO_POS - pos))]
    cur_pos = pos
    for i in range(int(num_slices.get())-1):
        i += 1
        cur_pos = float(cur_pos) + i*float(separation_slices.get())*(-1)**i
        if cur_pos not in POS_VEC and cur_pos <= float(np.max(MO_POS)) and cur_pos >= float(np.min(MO_POS)):
            POS_VEC.append(cur_pos)
        IDX.append(np.argmin(np.abs(MO_POS - POS_VEC[-1])))
    slice.set_data(WL_MAT[IDX[1:],:], COUNTS_MAT[IDX[1:],:])
    center_slice.set_data(WL_MAT[IDX[0],:], COUNTS_MAT[IDX[0],:])
    ax.set_xlim([float(xlim_min.get()), float(xlim_max.get())])
    # ax.set_xlim([np.min(WL_MAT[idx,:]), np.max(WL_MAT[idx,:])])
    # ax.set_ylim([0, np.max(COUNTS_MAT[idx,:])*1.1])
    canvas.draw()

def start_calibration():
    # L_fit = WL_MAT[0] <= float(xlim_max.get())
    # L_fit = L_fit & (WL_MAT[0] >= float(xlim_min.get()))

    # wl_fit = WL_MAT[L_fit][:]
    # counts_fit = COUNTS_MAT[:][L_fit]
    # counts_ini = COUNTS_ENV[:][L_fit]
    if num_slices.get() < 2:
        print('Need at least 2 slices for calibration')
        return
        

    pos_fit = MO_POS[IDX]
    wl_fit = []
    counts_fit = []
    L_fit = np.array(WL_MAT[0]) <= float(xlim_max.get())
    L_fit = L_fit & (np.array(WL_MAT[0]) >= float(xlim_min.get()))
    counts_ini = np.array(COUNTS_ENV[L_fit])
    for i in range(len(IDX)):
        counts_fit.append(COUNTS_MAT[IDX[i]][L_fit])
        wl_fit.append(WL_MAT[IDX[i]][L_fit])
    
    ps = pulse_shaper_obj(name='Dummy_PS')
    ps.calibrate(pos_fit, wl_fit, counts_fit, counts_ini, reference_pulse = None)
    # values read -> form pulse object 
    #apply the calibrated filter

    #initial_pulse_obj = pg.PulseGenerator(0,200,0.01,central_wavelength=np.mean(wl_fit[0]))
    initial_pulse_obj = pg.PulseGenerator(0,central_wavelength=np.mean(wl_fit[0]),
                                          f0=np.min(wl_fit[0]),fend=np.max(wl_fit[0]),fN = len(wl_fit[0]))
                                    
    initial_pulse_obj.add_spectrum_frequ(wl_fit[0], counts_ini, power = None,plot=False)
    
    def _smoothing_fkt(var): 
        global counts_ini_sm
        global counts_fit_sm

        counts_ini_sm = np.copy(counts_ini)
        counts_fit_sm = np.copy(counts_fit)

        smooth_window = int(len(counts_ini)*smooth_slider.get()*0.2)
        counts_ini_sm[counts_ini < background_slider.get()] = 0
        if smooth_window > 3:
            counts_ini_sm = savgol_filter(counts_ini_sm,smooth_window,3)
            counts_ini_sm[counts_ini_sm < background_slider.get()] = 0
            #counts_ini_sm = savgol_filter(counts_ini_sm,smooth_window,3)
            #counts_ini_sm[counts_ini_sm < background_slider.get()] = 0
            
        # for i in range(len(counts_fit)):
        #     counts_fit_sm[i][counts_fit[i] < background_slider.get()] = 0
        #     if smooth_window > 3:
        #         counts_fit_sm[i] = savgol_filter(counts_fit_sm[i],smooth_window,3)
        #     counts_fit_sm[1][counts_fit_sm[i] < background_slider.get()] = 0
        
        initial_pulse_obj.clear_all()
        initial_pulse_obj.add_spectrum_frequ(wl_fit[0], counts_ini_sm, power = None,plot=False)

        shaped_pulse = ps.move_slit(float(move_slider.get()),pulse_object=initial_pulse_obj)
        slice_index = int(np.argmin(np.abs(MO_POS - ps.motor_position)))
        env_ca.set_data(wl_fit[0], counts_ini - background_slider.get())  
        env_sm.set_data(initial_pulse_obj.wavelengths, np.abs(initial_pulse_obj.frequency_representation_x)**2)
        ax_time_env.set_data(np.real(initial_pulse_obj.time), np.abs(initial_pulse_obj.temporal_representation_x))
        ax_time_shaped.set_data(np.real(shaped_pulse.time),np.abs(shaped_pulse.temporal_representation_x))
        center_slice_ca.set_data(WL_MAT[0], COUNTS_MAT[slice_index] - background_slider.get())
        center_slice_sm.set_data(shaped_pulse.wavelengths,np.abs(shaped_pulse.frequency_representation_x)**2)


        canvas_cal.draw()
        canvas_time.draw()

    def _save_calibration():

        def _save_helper():
            ps.save_calibration(file_name = save_name.get())
            save_cal_window.destroy()
        save_cal_window = tk.Toplevel(calibration_window)
        save_cal_window.title('Save calibration')
        save_name = tk.Entry(save_cal_window,width=50)
        save_name.insert(0, 'calibration_'+datetime.datetime.now().strftime("%Y%m%d%H%M"))
        save_button = tk.Button(save_cal_window, text="Save", command=_save_helper)
        save_name.pack(side=tk.TOP)
        save_button.pack(side=tk.TOP)

        

    def _save_initial_pulse():
        initial_pulse_save = initial_pulse_obj.copy_pulse()
        save_window = tk.Toplevel(root)
        save_window.title('Save initial pulse')

        fig_save = Figure(figsize=(5, 4), dpi=100)
        ax_save = fig_save.add_subplot()
        ax_save_env, = ax_save.plot(np.real(initial_pulse_save.time), np.abs(initial_pulse_save.temporal_representation_x),'k-',alpha=1)
        ax_save_field, = ax_save.plot(np.real(initial_pulse_save.time), np.real(initial_pulse_save.temporal_representation_x),'k-',alpha=0.5)
        ax_save.set_xlabel('Time (ps)')
        ax_save.set_xlim([np.min(np.real(initial_pulse_save.time)), np.max(np.real(initial_pulse_save.time))])

        canvas_save = FigureCanvasTkAgg(fig_save, master=save_window)  # A tk.DrawingArea.
        canvas_save.draw()
        canvas_save.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True) 
        # update function 
        def _update_save(var):
            initial_clone = initial_pulse_obj.copy_pulse()
            initial_pulse_save = pg.PulseGenerator(float(t0_entry.get()),float(tend_entry.get()),float(dt_entry.get()),central_wavelength=initial_pulse_obj.central_wavelength)
            initial_pulse_save.clear_all()
            initial_clone.clear_filter()
            initial_clone.add_filter_rectangle()
            time_shift = float(shift_slider.get())-(initial_pulse_obj.t0+initial_pulse_obj.tend)/2
            initial_clone.add_phase_filter(phase_taylor = [0,time_shift])
            initial_clone.apply_frequency_filter()

            initial_pulse_save.merge_pulses(initial_clone)
            if qd_calibration_text.cget("text") != 'None':
                initial_pulse_save.set_rotating_frame(qd_calibration_text.cget("text"))
            else:
                initial_pulse_save.set_rotating_frame(float(central_wavelength_entry.get()))

            #initial_pulse_save.add_filter_rectangle()
            #time_shift = float(shift_slider.get())-(initial_pulse_obj.t0+initial_pulse_obj.tend)/2
            #initial_pulse_save.add_phase_wedge(time_shift, central_f = 0, shift_time = True, polarisation = 'b',kind = 'r')
        
            

            ax_save_env.set_data(np.real(initial_pulse_save.time), np.abs(initial_pulse_save.temporal_representation_x))
            ax_save_field.set_data(np.real(initial_pulse_save.time), np.real(initial_pulse_save.temporal_representation_x))
            ax_save.set_xlim([np.min(np.real(initial_pulse_save.time)), np.max(np.real(initial_pulse_save.time))])
            shift_slider.config(from_=initial_pulse_save.t0, to=initial_pulse_save.tend, resolution=initial_pulse_save.dt)
            canvas_save.draw()

            if var == 1:

                def _save_helper():
                    initial_pulse_save.set_pulse_power(1)
                    initial_pulse_save.clear_filter()
                    initial_pulse_save.save_pulse(save_name = save_name.get())
                    print('Pulse saved as:' + save_name.get()+'.pulse')
                    save_cal_window.destroy()
                save_cal_window = tk.Toplevel(calibration_window)
                save_cal_window.title('Save Pulse')
                save_name = tk.Entry(save_cal_window,width=50)
                save_name.insert(0, 'pulse_'+datetime.datetime.now().strftime("%Y%m%d%H%M"))
                save_button = tk.Button(save_cal_window, text="Save", command=_save_helper)
                save_name.pack(side=tk.TOP)
                save_button.pack(side=tk.TOP)

                
                

        def load_calibration():
            old_dir = os.getcwd()
            calib_file = tk.filedialog.askopenfilename()
            initial_pulse_save.set_rotating_frame(calib_file)
            #dummy_pu = pg.PulseGenerator(0,2,1,calibration_file=calib_file)
            #central_wavelength_entry.insert(0, str(dummy_pu.central_wavelength))
            central_wavelength_entry.delete(0, tk.END)
            central_wavelength_entry.insert(0, str(initial_pulse_save.central_wavelength))
            os.chdir(old_dir)
            qd_calibration_text.config(text=calib_file)
            _update_save(0)

        # fields and sliders 
        t0_entry = tk.Entry(save_window,text="t0")
        t0_entry.insert(0, str(np.round(initial_pulse_save.t0,decimals=2)))
        tend_entry = tk.Entry(save_window,text="tend")
        tend_entry.insert(0, str(np.round(initial_pulse_save.tend,decimals=2)))
        dt_entry = tk.Entry(save_window,text="dt")
        dt_entry.insert(0, str(np.round(initial_pulse_save.dt,decimals = 2)))
        central_wavelength_entry = tk.Entry(save_window,text="central_wavelength")
        central_wavelength_entry.insert(0, str(np.round(initial_pulse_save.central_wavelength,decimals = 2)))


        shift_slider = tk.Scale(save_window, from_=initial_pulse_obj.t0, to=initial_pulse_obj.tend, orient=tk.HORIZONTAL,
                                command=_update_save,label="time shift",resolution=initial_pulse_obj.dt,length=300)
        shift_slider.set((initial_pulse_obj.t0+initial_pulse_obj.tend)/2)

        update_button = tk.Button(save_window, text="Update", command=lambda: _update_save(0))
        load_calib_button = tk.Button(save_window, text="Load QD calibration", command=load_calibration)
        save_button = tk.Button(save_window, text="Save", command=lambda: _update_save(1))
        qd_calibration_text = tk.Label(save_window, text="None")


        t0_entry.pack(side=tk.BOTTOM)
        tend_entry.pack(side=tk.BOTTOM)
        dt_entry.pack(side=tk.BOTTOM)
        central_wavelength_entry.pack(side=tk.BOTTOM)
        shift_slider.pack(side=tk.BOTTOM)
        qd_calibration_text.pack(side=tk.BOTTOM)
        load_calib_button.pack(side=tk.TOP)
        update_button.pack(side=tk.TOP)
        save_button.pack(side=tk.TOP)


        #initial_pulse_save.save_pulse()
        
        

    calibration_window = tk.Toplevel(root)

    shaped_pulse = ps.move_slit(794.2,pulse_object=initial_pulse_obj)

    # window design
    calibration_window.title('Calibration settings')
    fig_cal = Figure(figsize=(5, 4), dpi=100)
    ax_env = fig_cal.add_subplot() 
    ax_env.title.set_text('Spectral representation')
    env_ca, = ax_env.plot(wl_fit[0], counts_ini,'k-',alpha=0.5)
    env_sm, = ax_env.plot(initial_pulse_obj.wavelengths, np.abs(initial_pulse_obj.frequency_representation_x)**2,'k-',alpha=1)
    # for i in range(len(IDX)-1):
    #     ax_env.plot(wl_fit[i+1], counts_fit[i+1],'r-',alpha=0.5)
    center_slice_ca, = ax_env.plot(wl_fit[0], counts_fit[0],'r-',alpha=0.5)
    center_slice_sm, = ax_env.plot(shaped_pulse.wavelengths,np.abs(shaped_pulse.frequency_representation_x)**2,'r-',alpha=1)
    ax_env.set_xlabel('Wavelength (nm)')
    ax_env.set_ylabel('Counts')
    ax_env.set_xlim([np.min(wl_fit[0]), np.max(wl_fit[0])])

    canvas_cal = FigureCanvasTkAgg(fig_cal, master=calibration_window)  # A tk.DrawingArea.
    canvas_cal.draw()
    canvas_cal.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

    fig_time = Figure(figsize=(5, 4), dpi=100)
    ax_time = fig_time.add_subplot()
    ax_time.title.set_text('Temporal representation')
    ax_time_env, = ax_time.plot(np.real(initial_pulse_obj.time), np.abs(initial_pulse_obj.temporal_representation_x),'k-',alpha=1)
    ax_time_shaped, = ax_time.plot(np.real(shaped_pulse.time),np.abs(shaped_pulse.temporal_representation_x),'r-',alpha=1)
    ax_time.set_xlabel('Time (ps)')

    canvas_time = FigureCanvasTkAgg(fig_time, master=calibration_window)  # A tk.DrawingArea.
    canvas_time.draw()
    canvas_time.get_tk_widget().pack(side=tk.RIGHT, fill=tk.BOTH, expand=True) 

    background_slider = tk.Scale(calibration_window, from_=0, to=np.max(counts_ini), orient=tk.HORIZONTAL,
                                 command=_smoothing_fkt,label="Background",resolution=1,length=300)
    smooth_slider = tk.Scale(calibration_window, from_=0, to=1, orient=tk.HORIZONTAL,
                             command=_smoothing_fkt,label="Smoothing",resolution=0.01,length=300)
    move_slider = tk.Scale(calibration_window, from_=float(xlim_min.get()), to=float(xlim_max.get()), orient=tk.HORIZONTAL,
                                command=_smoothing_fkt,label="Move slit",resolution=np.abs(float(xlim_max.get())-float(xlim_min.get()))*1e-3,length=300)

    save_cal_button = tk.Button(calibration_window, text="Save calibration", command=_save_calibration)
    save_initioal_pulse_button = tk.Button(calibration_window, text="Save initial pulse", command=_save_initial_pulse)

    background_slider.set(0)
    move_slider.set(np.mean(wl_fit[0]))
    move_slider.pack(side=tk.BOTTOM)
    background_slider.pack(side=tk.BOTTOM)
    smooth_slider.pack(side=tk.BOTTOM)
    save_cal_button.pack(side=tk.RIGHT)
    save_initioal_pulse_button.pack(side=tk.RIGHT)

  
    pass

#connect_spectrometer = tk.Button(root, text="Connect to Spectrometer", command=connect_spectrometer)
#connect_motor = tk.Button(root, text="Connect to Motor", command=connect_motor)

connected_spectrometer_text = tk.Label(root, text="Spectrometer connected: "+SPEC.name)
connected_motor_text = tk.Label(root, text="Motor connected: "+MOTOR.name)

slider_slice = tk.Scale(root, from_=0, to=1, orient=tk.HORIZONTAL,
                               command=draw_slice,label="motor position",resolution=0.1,length=300)
num_slices = tk.Scale(root, from_=2, to = 10, orient=tk.HORIZONTAL,
                               command = draw_slice,resolution=1,length=300)
separation_slices = tk.Entry(root, text="Separation between slices")
#separation_slices.insert()
num_slices.set(2)

xlim_min = tk.Entry(root)
xlim_max = tk.Entry(root)


aquire_button = tk.Button(root, text="Aquire data", command=aquire_data)
start_cal_button = tk.Button(root, text="Start calibration", command=start_calibration)

# start_cal_button.pack(side=tk.TOP)
# #connect_spectrometer.pack(side=tk.TOP)
# #connect_motor.pack(side=tk.TOP)
# connected_spectrometer_text.pack(side=tk.TOP)
# connected_motor_text.pack(side=tk.TOP)
# aquire_button.pack(side=tk.TOP)
# canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
# slider_slice.pack(side=tk.BOTTOM)
# xlim_max.pack(side=tk.BOTTOM)
# xlim_min.pack(side=tk.BOTTOM)
# num_slices.pack(side=tk.BOTTOM)
# separation_slices.pack(side=tk.BOTTOM)

aquire_button.grid(row=0, column=0)
start_cal_button.grid(row=0, column=1)
start_cal_button.config(state=tk.DISABLED)
aquire_button.grid(row=0, column=0)

connected_spectrometer_text.grid(row=1, column=0,columnspan=2)
connected_motor_text.grid(row=2, column=0, columnspan=2)
canvas.get_tk_widget().grid(row=3, column=0, columnspan=2)

tk.Label(root, text = 'Calibration min (nm): ').grid(row=4, column=0)
xlim_min.grid(row=4, column=1)
xlim_min.config(state=tk.DISABLED)

tk.Label(root, text = 'Calibration max (nm): ').grid(row=5, column=0)
xlim_max.grid(row=5, column=1)
xlim_max.config(state=tk.DISABLED)

tk.Label(root, text = 'Slice separation: ').grid(row=6, column=0)
separation_slices.grid(row=6, column=1)
separation_slices.config(state=tk.DISABLED)

tk.Label(root, text = 'Number of slices: ').grid(row=7, column=0)
num_slices.grid(row=7, column=1)
num_slices.config(state=tk.DISABLED)

tk.Label(root, text = 'Center slice: ').grid(row=8, column=0)
slider_slice.grid(row=8, column=1)
slider_slice.config(state=tk.DISABLED)


def enable_gui():
    start_cal_button.config(state=tk.NORMAL)
    aquire_button.config(state=tk.NORMAL)
    xlim_min.config(state=tk.NORMAL)
    xlim_max.config(state=tk.NORMAL)
    separation_slices.config(state=tk.NORMAL)
    num_slices.config(state=tk.NORMAL)
    slider_slice.config(state=tk.NORMAL)
    
    


tk.mainloop()


