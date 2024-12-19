import numpy as np
from matplotlib import pyplot as plt
from scipy.io import savemat
import pyaceqd.pulsegenerator as pg 
import os as os
from Pulse_v2 import fake_spectrometer, fake_motor, pulse_shaper_obj, simulator, attenuator, half_wave_plate, fake_attenuator, motor, spectrometer, load_pulse_device, create_experiment, save_pulse, save_temp, save_device, excecute_folder, power_meter, time_delay, fake_power_meter, polarising_beam_splitter, chirp_filter, generic_wave_plate 

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
qd_calibration_2  = 'NW_QD2.txt' 

motor_A = motor(fake_motor(),name='motor_A')
attenuator_A = attenuator(fake_attenuator(),name='attenuator_A_+45')
pulse_shaper_A = pulse_shaper_obj(device=motor_A, name='pulse_shaper_A_+45',calibration_file=ps_calibration)

motor_B = motor(fake_motor(),name='motor_B')
attenuator_B = attenuator(fake_attenuator(),name='attenuator_B_-45')
pulse_shaper_B = pulse_shaper_obj(device=motor_B, name='pulse_shaper_B_-45',calibration_file=ps_calibration2)

motor_C = motor(fake_motor(),name='motor_C')
delay_stage = time_delay(motor_C,name='delay_stage_C', offset=5)

motor_D = motor(fake_motor(),name='motor_D')
wave_plate_D = generic_wave_plate(motor_D,name='-45_HWP')

motor_E = motor(fake_motor(),name='motor_E')
wave_plate_E = generic_wave_plate(motor_E,name='+45_HWP')

motor_F = motor(fake_motor(),name='motor_F')
wave_plate_F = generic_wave_plate(motor_F,name='-45_QWP', phase=np.pi/2)

motor_G = motor(fake_motor(),name='motor_G')
wave_plate_G = generic_wave_plate(motor_G,name='+45_QWP', phase=np.pi/2)

motor_H = motor(fake_motor(),name='motor_H')
wave_plate_H = generic_wave_plate(motor_H,name='0_HWP')


spectrometer_M1 = spectrometer(fake_spectrometer(start_wl=894, end_wl = 902, n_wl = 1340),name='spectrometer_M1', polarisation_angle= None)

power_meter_M2 = power_meter(fake_power_meter(),name='power_meter_A')
power_meter_M3 = power_meter(fake_power_meter(),name='power_meter_B')

# ----- initial pulse -----
t_0 = 0
t_end = 150
dt = 0.2#0.02

chirp_A = +0
chirp_B = -0

initial_wavelength_A = 897.45 # tpe pulse
initial_wavelength_B = 898 # stim pulse

initial_pulse = pg.PulseGenerator(t_0,t_end,dt,calibration_file=qd_calibration)
initial_pulse.add_gaussian_time(unit = 'nm', central_f= 897.5, width_t=1,t0=t_end/2-30,area_time=20, polarisation=[2,1]) # <-- more power for tpe 

#wp_0 = wave_plate_H.open_control(pulse_object=initial_pulse)

initial_pulse_A, initial_pulse_B = polarising_beam_splitter(initial_pulse)

initial_pulse_A = chirp_filter(initial_pulse_A, chirp_A, central_wavelength=initial_wavelength_A)
initial_pulse_B = chirp_filter(initial_pulse_B, chirp_B, central_wavelength=initial_wavelength_B)

# ---- pulse generator controls -----

ipA_control = ip.Pulse_Generator(initial_pulse = initial_pulse_A, name='initial_pulse_A')
ipB_control = ip.Pulse_Generator(initial_pulse = initial_pulse_B, name='initial_pulse_B')
# ----- pulse shaper A -----
ps_control_A = pulse_shaper_A.open_control(previous_control=ipA_control,open_gui=False)
ps_control_A.set_control_value(initial_wavelength_A)
ps_control_A.update_control()
ps_control_A.set_step_size(0.025)
ps_control_A.gui()

wave_plate_control_A = wave_plate_E.open_control(previous_control=ps_control_A)
wave_plate_control_G = wave_plate_G.open_control(previous_control=wave_plate_control_A)

att_control_A = attenuator_A.open_control(previous_control=wave_plate_control_G)

# ----- pulse shaper B -----
ps_control_B = pulse_shaper_B.open_control(previous_control=ipB_control,open_gui=False)
ps_control_B.set_control_value(initial_wavelength_B)
ps_control_B.update_control()
ps_control_B.set_step_size(0.025)
ps_control_B.gui()
att_control_B = attenuator_B.open_control(previous_control=ps_control_B, open_gui=False)
att_control_B.set_control_value(0.0)
att_control_B.gui()

# ----- wave plate -----
wave_plate_control = wave_plate_D.open_control(previous_control=att_control_B)
wave_plate_control_qwp = wave_plate_F.open_control(previous_control=wave_plate_control)
# ----- delay stage -----
delay_control_C = delay_stage.open_control(previous_control=wave_plate_control_qwp)

# ----- simulator -----
sim = simulator(qd_calibration=qd_calibration, sim_kind='ace', decay=True,name='QD1',temp_dir='sim_dump/', phonons=False)
sim_control = sim.open_control(previous_control=[att_control_A,delay_control_C], open_gui=False)
sim_control.set_dipole_moment(1.3)
sim_control.gui()

sim2 = simulator(qd_calibration=qd_calibration_2, sim_kind=None, decay=True, name='QD2',temp_dir='sim_dump/',phonons=False)
sim_control2 = sim2.open_control(previous_control=[att_control_A,delay_control_C], open_gui=False)
sim_control2.set_dipole_moment(1.5)
sim_control2.gui()

# ----- spectrometer -----
spectrometer_control_M1 = spectrometer_M1.open_control(previous_control=[att_control_A,delay_control_C], open_gui=False, simulation_control= [sim_control,sim_control2])
spectrometer_control_M1.toggle_display_experiment()
spectrometer_control_M1.toggle_display_pulse()
spectrometer_control_M1.toggle_display_simulation()

#spectrometer_control_M1.set_measurement_method('single line')
#spectrometer_control_M1.set_measurement_arguments(arguments=[896.5, 0.05, 'max']) 

spectrometer_control_M1.set_measurement_method('multi line', arg=2)
spectrometer_control_M1.set_measurement_arguments(arguments=[2,896.5, 0.05, 'max','',896.8, 0.05, 'max','']) 

#spectrometer_control_M1.toggle_autoscale()
spectrometer_control_M1.set_simulation_counts(5000)
spectrometer_control_M1.set_simulation_background(60)
spectrometer_control_M1.set_simulation_gaussian_noise(10)
spectrometer_control_M1.set_poissonian_noise(True)
spectrometer_control_M1.set_pulse_scale(10)

spectrometer_control_M1.gui()
spectrometer_control_M1.ax_simulation.set_xlim([896, 900])
# ----- power meter -----
power_meter_control_M2 = power_meter_M2.open_control(previous_control=[att_control_A])
power_meter_control_M3 = power_meter_M3.open_control(previous_control=[delay_control_C])

# ----- PULSE tools -----
# ----- control optimizer -----
optimizer = con_opt.control_optimizer(device_control=[ps_control_A,att_control_A, ps_control_B, att_control_B, delay_control_C, wave_plate_control_A, wave_plate_control_G, wave_plate_control, wave_plate_control_qwp], measururement_control=[spectrometer_control_M1])

# ----- hyper dimensional scan -----

scanner = hds.hyper_scan(device_control=[ps_control_A,att_control_A, ps_control_B, att_control_B, delay_control_C], measururement_control=[spectrometer_control_M1, power_meter_control_M2, power_meter_control_M3])

# ------ gui manager --------
manager = gui_manager(device_control=[ps_control_A, att_control_A, sim_control,ps_control_B, att_control_B, sim_control2, delay_control_C,  spectrometer_control_M1, power_meter_control_M2, power_meter_control_M3, optimizer,scanner,ipA_control, ipB_control, wave_plate_control_A, wave_plate_control_G,wave_plate_control,wave_plate_control_qwp], num_coloums=3)





spectrometer_control_M1.start_gui()