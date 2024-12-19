# fake lab devices, can be replaced by real devices
from Pulse_v2 import fake_motor, fake_attenuator, fake_spectrometer, fake_power_meter

# PULSE devices and utilities
from Pulse_v2 import motor, attenuator, spectrometer, power_meter, pulse_shaper_obj
import control_optimizer as con_opt
import hyper_dimensional_scan as hds
from PULSE_gui_tools import gui_manager

# change the working directory --- maybe you don't need this
import os
cur_dir = os.getcwd()
if 'PulseGenerationACE' in cur_dir:
    pass
else:
# change the directory to the folder where the calibration files are stored
    os.chdir(cur_dir+'/PulseGenerationACE/PULSE')
print(cur_dir) 


# connect your lab devices here -> relpace the fake devices
lab_motor = fake_motor()
lab_attenuator = fake_attenuator()
lab_spectrometer = fake_spectrometer(start_wl=894, end_wl = 902, n_wl = 1340)
lab_power_meter = fake_power_meter()

# create PULSE devices
pulse_motor = motor(lab_motor, name='pulse_motor')
pulse_attenuator = attenuator(lab_attenuator, name='pulse_attenuator')
pulse_spectrometer = spectrometer(lab_spectrometer, name='pulse_spectrometer')
pulse_power_meter = power_meter(lab_power_meter, name='pulse_power_meter')

# create pulse shaper
calibration_file = 'calibration_slit_13_900.txt' # <- replace with your calibration file
pulse_shaper = pulse_shaper_obj(device=pulse_motor, name='pulse_shaper', calibration_file=calibration_file)

# open controlers (GUIs) 
ps_controler = pulse_shaper.open_control()
att_controler = pulse_attenuator.open_control(previous_control=ps_controler)
spectrometer_controler = pulse_spectrometer.open_control(previous_control=att_controler)
power_meter_controler = pulse_power_meter.open_control(previous_control=att_controler)

# open utilities
optimizer = con_opt.control_optimizer(device_control=[ps_controler, att_controler],measururement_control=[spectrometer_controler])
scanner = hds.hyper_scan(device_control=[ps_controler, att_controler],measururement_control=[spectrometer_controler, power_meter_controler])

# open gui manager
gm = gui_manager([ps_controler, att_controler, spectrometer_controler, power_meter_controler, optimizer, scanner])

# start mainloop
ps_controler.start_gui()