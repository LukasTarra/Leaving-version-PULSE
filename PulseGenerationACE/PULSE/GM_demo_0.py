from Pulse_v2 import fake_spectrometer, fake_motor, fake_attenuator, motor, spectrometer, attenuator, generic_wave_plate, pulse_shaper_obj, simulator, generic_wave_plate 
from PULSE_gui_tools import gui_manager
import control_optimizer as con_opt
import os as os
import pyaceqd.pulsegenerator as pg

cur_dir = os.getcwd()
if 'PulseGenerationACE' in cur_dir:
    pass
else:
# change the directory to the folder where the calibration files are stored
    os.chdir(cur_dir+'/PulseGenerationACE/PULSE')
print(cur_dir) 

ps_calibration = 'calibration_slit_13.txt' 
qd_calibration = 'QD_Iker_April.txt'

lab_motor = motor(fake_motor(),name='motor1')
lab_att = fake_attenuator()
lab_spectromter = fake_spectrometer(start_wl=774.15, end_wl = 793.422, n_wl = 1340)
lab_hwp = motor(fake_motor(),name='hwp')

pulse_shaper = pulse_shaper_obj(device=lab_motor,calibration_file=ps_calibration, name='Pulse Shaper')
att = attenuator(fake_attenuator(), name='Attenuator')
spec = spectrometer(device=lab_spectromter, name='Spectrometer')
hwp = generic_wave_plate(device=lab_hwp, name='HWP')
sim = simulator(qd_calibration=qd_calibration,name='Simulator',temp_dir="sim_dump/")

initial_pulse = pg.PulseGenerator(t0=0, tend=50, dt=0.2, calibration_file=qd_calibration)

initial_pulse.add_gaussian_time(unit='meV', central_f= 0, sig_or_fwhm='fwhm', width_t=0.3, t0= 25, polarisation=[1,0])

initial_pulse.set_pulse_power(400)

ps_control = pulse_shaper.open_control(pulse_object=initial_pulse, open_gui=False)
ps_control.set_step_size(0.05)
ps_control.gui()
hwp_control = hwp.open_control(previous_control=ps_control)
att_control = att.open_control(previous_control=hwp_control)
sim_control = sim.open_control(previous_control=att_control)
spec_control = spec.open_control(previous_control=att_control,simulation_control=sim_control, open_gui=False)
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
spec_control.gui()

opti_control = con_opt.control_optimizer([ps_control, att_control,hwp_control],measururement_control=spec_control)

gui = gui_manager([ps_control, att_control,hwp_control, spec_control,sim_control,opti_control])
gui.start_gui()
