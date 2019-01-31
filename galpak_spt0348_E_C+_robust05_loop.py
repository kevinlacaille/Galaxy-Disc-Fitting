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

#A more complicated run with setting min, max, and initial parameters
#using a guess from a previous run
min_bounds = galpak.GalaxyParameters(radius=0.5,inclination=30.0, velocity_dispersion=100.0, maximum_velocity=-350.0,turnover_radius=1e-5)
max_bounds = galpak.GalaxyParameters(radius=10.0,inclination=90.0, velocity_dispersion=300.0, maximum_velocity=350.0,turnover_radius=0.03)
initial_params = galpak.GalaxyParameters(x=22.57,y=27.07,z=32.45,flux=73.5,radius=7.80,inclination=63.3,pa=51.2,turnover_radius=0.005,maximum_velocity=177.0, velocity_dispersion=189.0)

#measure total time for all loops
t_loop_start = time.time()

#loop N times
N = 100
loop = range(N)[1:]
for i in loop:
	print 'loop #' + str(i)

	#time the run
	t_start = time.time()

	SPT0348_E = run('../spt0348_C+_dirty_contsub_briggs_robust05_E.fits', instrument=ALMA_b7, flux_profile='gaussian', redshift=5.652, min_boundaries=min_bounds, max_boundaries=max_bounds, initial_parameters = initial_params, random_scale = 6.0, verbose=False)

	#measure total time
	t_end = time.time()
	t_tot = t_end-t_start

	#tell me how long the run took
	print 'run took: ' + str(int(t_tot)) + ' seconds'


	#Record data
	print 'acceptance rate = ' + str(SPT0348_E.acceptance_rate) + ' %'  #should be ~30-50%
	print 'dynamical mass = ' + str(float(SPT0348_E.dynamical_mass)*1e-10) + ' x10^10 Msun' 
	SPT0348_E.save('galpak_SPT0348_E_loop'+str(i)+'_gauss')
	with open('galpak_SPT0348_E_loop'+str(i)+'_gauss_chain.dat','r') as chain_file:
		data = asciitable.read(chain_file.read(),Reader=asciitable.FixedWidth)
		print 'min chi^2 = ' +str(min(data.reduced_chi))

	'''
	#plot stuff
	SPT0348_E.plot_images()
	SPT0348_E.plot_correlations()
	SPT0348_E.plot_obs_vfield()
	SPT0348_E.plot_true_vfield()
	'''

#measure total time for loops
t_loop_end = time.time()
t_loop_tot = t_loop_end - t_loop_start
#tell me how long loops took
print 'loops took: ' + str(int(t_loop_tot/60.0)) + ' minutes'


'''
for i in loop:
	open('galpak_SPT0348_E_loop'+str(i)+'_gauss_chain.dat','r')
'''



















































































































 
