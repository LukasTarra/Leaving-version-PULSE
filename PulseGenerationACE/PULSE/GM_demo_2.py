import numpy as np
from matplotlib import pyplot as plt
from scipy.io import savemat
import pyaceqd.pulsegenerator as pg 
import os as os
from Pulse_v2 import fake_spectrometer, fake_motor, pulse_shaper_obj, simulator, attenuator, half_wave_plate, fake_attenuator, motor, spectrometer, load_pulse_device, create_experiment, save_pulse, save_temp, save_device, excecute_folder, power_meter, time_delay, fake_power_meter, polarising_beam_splitter, chirp_filter, fake_motor_alt, generic_wave_plate 

import control_optimizer as con_opt
import hyper_dimensional_scan as hds
import initialPulse_control as ip
import sys
from PULSE_gui_tools import gui_manager 

def blockPrint():
    sys.stdout = open(os.devnull, 'w')

#blockPrint()

cur_dir = os.getcwd()
if 'PulseGenerationACE' in cur_dir:
    pass
else:
# change the directory to the folder where the calibration files are stored
    os.chdir(cur_dir+'/PulseGenerationACE/PULSE')
print(cur_dir) 

# ----- Devices -----
ps_calibration = 'calibration_slit_13_900.txt' 
ps_calibration2 = 'calibration_slit_25_900.txt'
qd_calibration = 'NW_QD1.txt' 

motor_A = motor(fake_motor(),name='motor_A')
attenuator_A = attenuator(fake_attenuator(),name='attenuator_A_+45')
pulse_shaper_A = pulse_shaper_obj(device=motor_A, name='pulse_shaper_A_+45',calibration_file=ps_calibration)

motor_B = motor(fake_motor(),name='motor_B')
hwp = generic_wave_plate(motor_B,name='HWP', phase=np.pi)
qwp = generic_wave_plate(motor_B,name='QWP', phase=np.pi/2)


spectrometer_M1 = spectrometer(fake_spectrometer(start_wl=894, end_wl = 902, n_wl = 1340),name='spectrometer_M1', polarisation_angle= None)

power_meter_M2 = power_meter(fake_power_meter(),name='power_meter_A')


# ----- initial pulse -----
t_0 = 0
t_end = 300
dt = 0.2#0.02

chirp_A = +0
chirp_B = -0

initial_wavelength_A = 896.5 # tpe pulse
initial_wavelength_B = 898 # stim pulse

initial_pulse = pg.PulseGenerator(t_0,t_end,dt,calibration_file=qd_calibration)
initial_pulse.add_gaussian_time(unit = 'Hz', central_f= initial_pulse.tpe_resonance, width_t=0.2,t0=t_end/2,area_time=5, polarisation=[1,0]) # <-- more power for tpe 

initial_pulse.add_filter_double_erf(unit='nm', central_f= 896.5, width_f= 0.3, rise_f=0.05, invert=True)
initial_pulse.apply_frequency_filter()
#initial_pulse.clear_filter()
#initial_pulse.add_filter_double_erf(unit='nm', central_f= 898.1, width_f= 0.3, rise_f=0.05, invert=True)
#initial_pulse.apply_frequency_filter()
initial_pulse.clear_filter()
initial_pulse_A = initial_pulse


# ---- pulse generator controls -----

ipA_control = ip.Pulse_Generator(initial_pulse = initial_pulse_A, name='initial_pulse_A')

# ----- pulse shaper A -----
ps_control_A = pulse_shaper_A.open_control(previous_control=ipA_control,open_gui=False)
ps_control_A.set_control_value(initial_pulse._Units_inverse(initial_pulse.tpe_resonance,unit='nm'))
ps_control_A.update_control()
ps_control_A.set_step_size(0.025)
ps_control_A.set_control_value(896.5)
ps_control_A.gui()
att_control_A = attenuator_A.open_control(previous_control=ps_control_A)


# ----- wave plate -----
#hwp_control = hwp.open_control(previous_control=att_control_A)
#qwp_control = qwp.open_control(previous_control=hwp_control)


# ----- simulator -----
sim = simulator(qd_calibration=qd_calibration, sim_kind='ace', decay=True,name='QD1')
sim_control = sim.open_control(previous_control=[att_control_A], open_gui=False)
sim_control.set_dipole_moment(1.5)
sim_control.gui()


# ----- spectrometer -----
spectrometer_control_M1 = spectrometer_M1.open_control(previous_control=[att_control_A], open_gui=False, simulation_control= [sim_control])
spectrometer_control_M1.toggle_display_experiment()
spectrometer_control_M1.toggle_display_pulse()
spectrometer_control_M1.toggle_display_simulation()

spectrometer_control_M1.set_measurement_method('single line')
spectrometer_control_M1.set_measurement_arguments(arguments=[896.5, 0.05, 'max']) 

#spectrometer_control_M1.set_measurement_method('multi line', arg=2)
#spectrometer_control_M1.set_measurement_arguments(arguments=[2,896.5, 0.05, 'max','',896.8, 0.05, 'max','']) 

#spectrometer_control_M1.toggle_autoscale()
spectrometer_control_M1.set_simulation_counts(5000)



spectrometer_control_M1.set_pulse_scale(300)
spectrometer_control_M1.set_xlim(895, 898)
spectrometer_control_M1.toggle_autoscale()
spectrometer_control_M1.y_lim_min = 0
spectrometer_control_M1.y_lim_max = 5500

spectrometer_control_M1.gui()
#spectrometer_control_M1.ax_simulation.set_xlim([896, 900])

# ----- add spectrometer to ipA_control  and starting the gui-----
ipA_control.set_spectrometer_control(spectrometer_control_M1)
ipA_control.gui()

# ----- power meter -----
power_meter_control_M2 = power_meter_M2.open_control(previous_control=[att_control_A])


# ----- PULSE tools -----
# ----- control optimizer -----
optimizer = con_opt.control_optimizer(device_control=[ps_control_A,att_control_A], measururement_control=[spectrometer_control_M1])

# ----- hyper dimensional scan -----

scanner = hds.hyper_scan(device_control=[ps_control_A,att_control_A], measururement_control=[spectrometer_control_M1, power_meter_control_M2])

# ------ gui manager --------
manager = gui_manager(device_control=[ps_control_A, att_control_A, sim_control,  spectrometer_control_M1, power_meter_control_M2,  optimizer,scanner, ipA_control], num_coloums=3)





spectrometer_control_M1.start_gui()