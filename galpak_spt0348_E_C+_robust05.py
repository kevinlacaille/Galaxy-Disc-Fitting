import numpy
import scipy
import astropy
import matplotlib
import bottleneck
import galpak
import asciitable
from astropy.io import fits
from galpak import run
import time


#Set the beam and check parameters
#restoring beam = 0.749", 0.665", 5.826deg
ALMA_b7 = galpak.Instrument(psf=galpak.GaussianPointSpreadFunction(fwhm=0.749,pa=5.826,ba=float(0.665/0.749)),lsf=galpak.GaussianLineSpreadFunction(fwhm=1.0))

'''
#The most simple run you can do
SPT0348_E = galpak.GalPaK3D('spt0348_C+_dirty_contsub_briggs_robust05_E.fits', instrument=ALMA_b7, redshift=5.656)
'''

#A more complicated run with setting min, max, and initial parameters
#using a guess from a previous run
min_bounds = galpak.GalaxyParameters(radius=0.5,inclination=15.0, velocity_dispersion=50.0, maximum_velocity=-500.0,turnover_radius=1e-5)
max_bounds = galpak.GalaxyParameters(radius=10.0,inclination=90.0, velocity_dispersion=500.0, maximum_velocity=500.0,turnover_radius=1.0)
initial_params = galpak.GalaxyParameters(x=22.57,y=27.07,z=32.45,flux=73.5,radius=7.80,inclination=63.3,pa=51.2,turnover_radius=0.005,maximum_velocity=177.0, velocity_dispersion=189.0)

#time the run
t_start = time.time()

SPT0348_E = run('../spt0348_C+_dirty_contsub_briggs_robust05_E.fits', instrument=ALMA_b7, flux_profile='exponential', redshift=5.652, min_boundaries=min_bounds, max_boundaries=max_bounds, initial_parameters = initial_params, random_scale = 2.0, max_iterations=int(1*15000))

#SPT0348_W = run('../spt0348_C+_dirty_contsub_briggs_robust05_W.fits', instrument=ALMA_b7, flux_profile='exponential', redshift=5.653, min_boundaries=min_bounds, max_boundaries=max_bounds, initial_parameters = initial_params,random_scale = 2.0,max_iterations=int(2*15000))


#measure total time
t_end = time.time()
t_tot = t_end-t_start

#tell me how long the run took
print 'run took: ' + str(int(t_tot/60.0)) + ' minutes'

'''
#A simple run setting only initial parameters
#using a guess from a previous run
SPT0348_E_simple.run_mcmc(random_scale=2.5,flux_profile='gaussian', redshift=5.656, initial_parameters=galpak.GalaxyParameters(x=22.55,y=27.05,z=32.55,flux=75.25,radius=7.85,inclination=67.5,pa=56.0,turnover_radius=0.015,maximum_velocity=173.0,velocity_dispersion=179.75))
'''


#Record data
print 'acceptance rate = ' + str(SPT0348_E.acceptance_rate) + ' %'  #should be ~30-50%
print 'dynamical mass = ' + str(float(SPT0348_E.dynamical_mass)*1e-10) + ' x10^10 Msun' 
SPT0348_E.save('galpak_SPT0348_E_run1_exp')
with open('galpak_SPT0348_E_run1_exp_chain.dat','r') as chain_file:
	data = asciitable.read(chain_file.read(),Reader=asciitable.FixedWidth)
	print 'min chi^2 = ' +str(min(data.reduced_chi))


#plot stuff
SPT0348_E.plot_images()
SPT0348_E.plot_correlations()
SPT0348_E.plot_obs_vfield()
SPT0348_E.plot_true_vfield()



