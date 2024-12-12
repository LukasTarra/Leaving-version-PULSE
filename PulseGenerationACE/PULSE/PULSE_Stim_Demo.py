import numpy as np
from matplotlib import pyplot as plt
from scipy.io import savemat
import pyaceqd.pulsegenerator as pg 
import os as os
from Pulse_v2 import fake_spectrometer, fake_motor, pulse_shaper_obj, simulator, attenuator, half_wave_plate, fake_attenuator, motor, spectrometer, load_pulse_device, create_experiment, save_pulse, save_temp, save_device, excecute_folder, power_meter, time_delay 

import control_optimizer as con_opt

import sys

# Disable
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
t_end = 50
dt = 0.2


###### set 1 +45 ps^2
chirp_1 = 0
chirp_2 = 0



initial_pulse = pg.PulseGenerator(t_0,t_end,dt,calibration_file=qd_calibration)
initial_pulse.add_gaussian_time(unit = 'nm', central_f= 779, width_t=0.1,t0=15,area_time=20, polarisation=[1,0])

initial_pulse2 = initial_pulse.copy_pulse()

initial_pulse.add_filter_rectangle()
initial_pulse.add_phase_filter(unit='nm', phase_taylor=[0,0,chirp_1],central_f=779)
initial_pulse.apply_frequency_filter()

initial_pulse2.add_filter_rectangle()
initial_pulse2.add_phase_filter(unit='nm', phase_taylor=[0,0,chirp_2],central_f=778)
initial_pulse2.apply_frequency_filter()




pulse_shaper = pulse_shaper_obj(device=lab_motor, calibration_file=ps_calibration, name='A: Pulse shaper')
att_object = attenuator(lab_att, name = 'B: Attenuator')

pulse_shaper2 = pulse_shaper_obj(device=lab_motor2, calibration_file=ps_calibration, name='C: Pulse shaper2')
att_object2 = attenuator(lab_att, name = 'D: Attenuator2')

delay_stage = time_delay(device=lab_motor3, name='F: Delay stage')



ps_controller = pulse_shaper.open_control(pulse_object=initial_pulse)
att_controller = att_object.open_control(previous_control=ps_controller)

ps_controller2 = pulse_shaper2.open_control(pulse_object=initial_pulse2)
att_controller2 = att_object2.open_control(previous_control=ps_controller2)


delay_controller = delay_stage.open_control(previous_control=att_controller2)


simulator_object = simulator(qd_calibration=qd_calibration, sim_kind='ace', temp_dir='sim_dump/')
sim_controller = simulator_object.open_control(previous_control = [att_controller,delay_controller],open_gui=False) #,delay_controller
sim_controller.gui()

spec = spectrometer(device=lab_spectromter, name='MF_spec') 
sm_controller = spec.open_control(pulse_object=None,previous_control = [att_controller,delay_controller], simulation_control= sim_controller,open_gui=False) #, ps_controller3#ps_controller2.get_pulse_object()
sm_controller.set_simulation_background(45)
sm_controller.set_simulation_gaussian_noise(5)#4.35
sm_controller.set_pulse_scale(0.08)
sm_controller.set_simulation_counts(600)
sm_controller.toggle_running()
sm_controller.gui()

sim_controller.toggle_running()
sim_controller.update_gui()

sm_controller.update_gui()
sm_controller.change_view()



co = con_opt.control_optimizer(device_control=[ ps_controller,att_controller, ps_controller2,att_controller2,delay_controller],measururement_control=sm_controller,measurement_kind='spectrometer',open_gui=False) # ,delay_controller , ps_controller,att_controller
#co.set_scan_limits()
co.set_spectrometer_measurement(779.89,0.1,weight='max')
co.gui()

sm_controller.start_gui()