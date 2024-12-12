import numpy as np
from matplotlib import pyplot as plt
from scipy.io import savemat
import pyaceqd.pulsegenerator as pg 
import os as os
from Pulse_v2 import fake_spectrometer, fake_motor, pulse_shaper_obj, simulator, attenuator, half_wave_plate, fake_attenuator, motor, spectrometer, load_pulse_device, create_experiment, save_pulse, save_temp, save_device, excecute_folder, power_meter, time_delay, fake_power_meter, polarising_beam_splitter, chirp_filter 

import control_optimizer as con_opt
import hyper_dimensional_scan as hds
import sys
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
#ps_calibration2 = 'calibration_slit_38.txt' 
qd_calibration = 'NW_QD1.txt' 
qd_calibration_2  = 'NW_QD2.txt' 

motor_A = motor(fake_motor(),name='motor_A')
attenuator_A = attenuator(fake_attenuator(),name='attenuator_A')
pulse_shaper_A = pulse_shaper_obj(device=motor_A, name='pulse_shaper_A',calibration_file=ps_calibration)

motor_B = motor(fake_motor(),name='motor_B')
attenuator_B = attenuator(fake_attenuator(),name='attenuator_B')
pulse_shaper_B = pulse_shaper_obj(device=motor_B, name='pulse_shaper_B',calibration_file=ps_calibration)

motor_C = motor(fake_motor(),name='motor_C')
delay_stage = time_delay(motor_C,name='delay_stage_C', offset=20)

spectrometer_M1 = spectrometer(fake_spectrometer(start_wl=894, end_wl = 902, n_wl = 1340),name='spectrometer_M1')

power_meter_M2 = power_meter(fake_power_meter(),name='power_meter_A')
power_meter_M3 = power_meter(fake_power_meter(),name='power_meter_B')

# ----- initial pulse -----
t_0 = 0
t_end = 200
dt = 0.2#0.02

chirp_A = +45
chirp_B = -45

initial_wavelength_A = 897.45 # tpe pulse
initial_wavelength_B = 898 # stim pulse

initial_pulse = pg.PulseGenerator(t_0,t_end,dt,calibration_file=qd_calibration)
initial_pulse.add_gaussian_time(unit = 'nm', central_f= 897, width_t=0.1,t0=t_end/2,area_time=20, polarisation=[2,1]) # <-- more power for tpe 

initial_pulse_A, initial_pulse_B = polarising_beam_splitter(initial_pulse)

initial_pulse_A = chirp_filter(initial_pulse_A, chirp_A, central_wavelength=initial_wavelength_A)
initial_pulse_B = chirp_filter(initial_pulse_B, chirp_B, central_wavelength=initial_wavelength_B)

# ----- pulse shaper A -----
ps_control_A = pulse_shaper_A.open_control(pulse_object=initial_pulse_A,open_gui=False)
ps_control_A.set_control_value(initial_wavelength_A)
ps_control_A.update_control()
ps_control_A.gui()
att_control_A = attenuator_A.open_control(previous_control=ps_control_A)

# ----- pulse shaper B -----
ps_control_B = pulse_shaper_B.open_control(pulse_object=initial_pulse_B,open_gui=False)
ps_control_B.set_control_value(initial_wavelength_B)
ps_control_B.update_control()
ps_control_B.gui()
att_control_B = attenuator_B.open_control(previous_control=ps_control_B, open_gui=False)
att_control_B.set_control_value(0.0)
att_control_B.gui()

# ----- delay stage -----
delay_control_C = delay_stage.open_control(previous_control=att_control_B)

# ----- simulator -----
sim = simulator(qd_calibration=qd_calibration, sim_kind=None, decay=True,name='QD1')
sim_control = sim.open_control(previous_control=[att_control_A,delay_control_C])

sim2 = simulator(qd_calibration=qd_calibration_2, sim_kind=None, decay=True, name='QD2')
sim_control2 = sim2.open_control(previous_control=[att_control_A,delay_control_C])

# ----- spectrometer -----
spectrometer_control_M1 = spectrometer_M1.open_control(previous_control=[att_control_A,delay_control_C], open_gui=False, simulation_control= [sim_control,sim_control2])
spectrometer_control_M1.toggle_display_experiment()
spectrometer_control_M1.toggle_display_pulse()
spectrometer_control_M1.toggle_display_simulation()

spectrometer_control_M1.set_measurement_method('single line')
spectrometer_control_M1.set_measurement_arguments(arguments=[896.5, 0.05, 'max']) 

spectrometer_control_M1.set_measurement_method('multi line', arg=2)
spectrometer_control_M1.set_measurement_arguments(arguments=[2,896.5, 0.05, 'max','',896.8, 0.05, 'max','']) 

spectrometer_control_M1.gui()
# ----- power meter -----
power_meter_control_M2 = power_meter_M2.open_control(previous_control=[att_control_A])
power_meter_control_M3 = power_meter_M3.open_control(previous_control=[delay_control_C])

# ----- PULSE tools -----
# ----- control optimizer -----
optimizer = con_opt.control_optimizer(device_control=[ps_control_A,att_control_A, ps_control_B, att_control_B, delay_control_C], measururement_control=[spectrometer_control_M1])

# ----- hyper dimensional scan -----

scanner = hds.hyper_scan(device_control=[ps_control_A,att_control_A, ps_control_B, att_control_B, delay_control_C], measururement_control=[spectrometer_control_M1, power_meter_control_M2, power_meter_control_M3])



spectrometer_control_M1.start_gui()