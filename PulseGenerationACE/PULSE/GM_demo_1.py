from Pulse_v2 import fake_spectrometer, fake_motor, fake_attenuator, motor, spectrometer, attenuator, generic_wave_plate, pulse_shaper_obj, simulator, generic_wave_plate, time_delay, polarising_beam_splitter, chirp_filter, power_meter, fake_power_meter 
from PULSE_gui_tools import gui_manager
import control_optimizer as con_opt
import os as os
import pyaceqd.pulsegenerator as pg
import hyper_dimensional_scan as hds

cur_dir = os.getcwd()
if 'PulseGenerationACE' in cur_dir:
    pass
else:
# change the directory to the folder where the calibration files are stored
    os.chdir(cur_dir+'/PulseGenerationACE/PULSE')
print(cur_dir) 

ps_calibration = 'calibration_slit_13.txt' 
ps_calibration_2 = 'calibration_slit_13.txt'
qd_calibration = 'QD_Iker_April.txt'

lab_motor = motor(fake_motor(),name='motor1')
lab_att = fake_attenuator()

lab_motor_2 = motor(fake_motor(),name='motor2')
lab_motor_3 = motor(fake_motor(),name='motor3')
lab_att_2 = fake_attenuator()

lab_power_meter = fake_power_meter()

lab_spectromter = fake_spectrometer(start_wl=774.15, end_wl = 793.422, n_wl = 1340)
lab_hwp = motor(fake_motor(),name='hwp')

pulse_shaper = pulse_shaper_obj(device=lab_motor,calibration_file=ps_calibration, name='A: Pulse Shaper')
att = attenuator(fake_attenuator(), name='A: Attenuator')

pulse_shaper_2 = pulse_shaper_obj(device=lab_motor_2,calibration_file=ps_calibration_2, name='B: Pulse Shaper')
att_2 = attenuator(fake_attenuator(), name='B: Attenuator')
time_delay_2 = time_delay(device=lab_motor_3, name='B: Time Delay', offset=30)

spec = spectrometer(device=lab_spectromter, name='Spectrometer')
hwp = generic_wave_plate(device=lab_hwp, name='B: HWP')
sim = simulator(qd_calibration=qd_calibration,name='Simulator',temp_dir="sim_dump/")
pow = power_meter(lab_power_meter, name='Power Meter')


initial_pulse = pg.PulseGenerator(t0=0, tend=200, dt=0.3, calibration_file=qd_calibration)

initial_pulse.add_gaussian_time(unit='meV', central_f= 0, sig_or_fwhm='fwhm', width_t=0.3, t0= 40, polarisation=[3,1])

initial_pulse.set_pulse_power(600)

initial_pulse_A, initial_pulse_B = polarising_beam_splitter(initial_pulse)

initial_pulse_A = chirp_filter(initial_pulse_A, 20,central_wavelength=780.8)
initial_pulse_B = chirp_filter(initial_pulse_B, -20, central_wavelength=781.8)
ps_control = pulse_shaper.open_control(pulse_object=initial_pulse_A, open_gui=False)
ps_control.set_step_size(0.05)
ps_control.set_control_value(780.965)
ps_control.gui()
att_control = att.open_control(previous_control=ps_control)

ps_control_2 = pulse_shaper_2.open_control(pulse_object=initial_pulse_B, open_gui=False)
ps_control_2.set_step_size(0.05)
ps_control_2.gui()
att_control_2 = att_2.open_control(previous_control=ps_control_2, open_gui=False)
att_control_2.set_current_attenuation(0)
att_control_2.gui()
time_delay_control = time_delay_2.open_control(previous_control=att_control_2)
hwp_control = hwp.open_control(previous_control=time_delay_control)



sim_control = sim.open_control(previous_control=[att_control,hwp_control],open_gui=False)
sim_control.toggle_decay()
sim_control.simulator_object.set_dipole_orientation(23)
sim_control.toggle_running()
sim_control.gui()
spec_control = spec.open_control(previous_control=[att_control,hwp_control],simulation_control=sim_control, open_gui=False)
spec_control.toggle_display_experiment()
spec_control.toggle_display_simulation()
spec_control.set_poissonian_noise(True)
spec_control.set_simulation_background(80)
spec_control.set_simulation_gaussian_noise(10)
spec_control.set_simulation_mode_combined(True)
spec_control.set_simulation_counts(4000)
spec_control.set_pulse_scale(3000)
spec_control.set_xlim(776,789)
spec_control.set_polarisation(85)
spec_control.set_measurement_arguments([779.89])
spec_control.toggle_autoscale()
spec_control.y_lim_min = 0
spec_control.y_lim_max = 4000
spec_control.gui()

pow_control = pow.open_control(previous_control=[att_control,hwp_control],open_gui=True)

opti_control = con_opt.control_optimizer([ps_control, att_control,ps_control_2,att_control_2,time_delay_control,hwp_control],measururement_control=spec_control)

hyp_scan = hds.hyper_scan([ps_control, att_control,ps_control_2,att_control_2,time_delay_control,hwp_control],measururement_control=[spec_control,pow_control])

gui = gui_manager([ps_control, att_control,ps_control_2,att_control_2,time_delay_control,hwp_control, spec_control,pow_control,sim_control,opti_control,hyp_scan],num_coloums=2)
gui.start_gui()
