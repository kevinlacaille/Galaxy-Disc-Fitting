from astropy.wcs import WCS
from astropy.io import fits
import matplotlib.pyplot as pl
import numpy as np
import math as m
from astropy.io.ascii.fastbasic import cparser
from astropy.table import np_utils



def arcsec_to_px(x):
    # wcs.wcs.crpix[0] = centre of image in pixel (x-coord)
    # wcs.wcs.cdelt[0] = size of pixel in arcsec
    return wcs.wcs.crpix[0]+x/wcs.wcs.cdelt[0]

def vel_to_px(v):
    # project the velocity to line of sight by *sin(inclination)
    # vel_proj=v*m.sin(m.radians(i))
    # shift image to the centre wcs.wcs.crpix[1]
    # divide data by velocity per pixel
    vel_proj=wcs.wcs.crpix[1]-v/wcs.wcs.cdelt[1]
    return vel_proj


dir = '/Users/kevinlacaille/Documents/University/PhD/Research/SPT0348/Data/ALMA/Cycle_4/GalPaK/galpak_august20_2018/'

E_gauss_arctan = np.loadtxt(dir + 'galpak_SPT0348_E_run1_gauss_arctan_true_Vrot.dat',unpack=True) #x,y,z,flux,radius,incl,pa,rV,Vmax,vel disp,reduced chi
dx_E_gauss_arctan = E_gauss_arctan[0] #arcsec
v_E_gauss_arctan = E_gauss_arctan[1] #km/s

E_gauss_tanh = np.loadtxt(dir + 'galpak_SPT0348_E_run1_gauss_tanh_true_Vrot.dat',unpack=True)
dx_E_gauss_tanh = E_gauss_tanh[0]
v_E_gauss_tanh = E_gauss_tanh[1]
#
E_gauss_exp = np.loadtxt(dir + 'galpak_SPT0348_E_run1_gauss_exp_true_Vrot.dat',unpack=True)
dx_E_gauss_exp = E_gauss_exp[0]
v_E_gauss_exp = E_gauss_exp[1]
# Exponential flux profile
E_exp_arctan = np.loadtxt(dir + 'galpak_SPT0348_E_run1_exp_arctan_true_Vrot.dat',unpack=True)
dx_E_exp_arctan = E_exp_arctan[0]
v_E_exp_arctan = E_exp_arctan[1]
#
E_exp_tanh = np.loadtxt(dir + 'galpak_SPT0348_E_run1_exp_tanh_true_Vrot.dat',unpack=True)
dx_E_exp_tanh = E_exp_tanh[0]
v_E_exp_tanh = E_exp_tanh[1]
#
E_exp_exp = np.loadtxt(dir + 'galpak_SPT0348_E_run1_exp_exp_true_Vrot.dat',unpack=True)
dx_E_exp_exp = E_exp_exp[0]
v_E_exp_exp = E_exp_exp[1]



W_gauss_tanh = np.loadtxt(dir + 'galpak_SPT0348_W_run1_gauss_tanh_true_Vrot.dat',unpack=True)
dx_W_gauss_tanh = W_gauss_tanh[0]
v_W_gauss_tanh = W_gauss_tanh[1]
#
# W_gauss_exp = np.loadtxt(dir + 'galpak_SPT0348_W_run1_gauss_exp_true_Vrot.dat',unpack=True)
# dx_W_gauss_exp = W_gauss_exp[0]
# v_W_gauss_exp = W_gauss_exp[1]
# Exponential flux profile
W_exp_arctan = np.loadtxt(dir + 'galpak_SPT0348_W_run1_exp_arctan_true_Vrot.dat',unpack=True)
dx_W_exp_arctan = W_exp_arctan[0]
v_W_exp_arctan = W_exp_arctan[1]
#
W_exp_tanh = np.loadtxt(dir + 'galpak_SPT0348_W_run1_exp_tanh_true_Vrot.dat',unpack=True)
dx_W_exp_tanh = W_exp_tanh[0]
v_W_exp_tanh = W_exp_tanh[1]
#
W_exp_exp = np.loadtxt(dir + 'galpak_SPT0348_W_run1_exp_exp_true_Vrot.dat',unpack=True)
dx_W_exp_exp = W_exp_exp[0]
v_W_exp_exp = W_exp_exp[1]

W_gauss_arctan = np.loadtxt(dir + 'galpak_SPT0348_W_run1_gauss_arctan_true_Vrot.dat',unpack=True) #x,y,z,flux,radius,incl,pa,rV,Vmax,vel disp,reduced chi
dx_W_gauss_arctan = W_gauss_arctan[0] #arcsec
v_W_gauss_arctan = W_gauss_arctan[1] #km/s

W_exp_arctan = np.loadtxt(dir + 'galpak_SPT0348_W_run1_exp_arctan_true_Vrot.dat',unpack=True) #x,y,z,flux,radius,incl,pa,rV,Vmax,vel disp,reduced chi
dx_W_exp_arctan = W_exp_arctan[0] #arcsec
v_W_exp_arctan = W_exp_arctan[1] #km/s

# inclination of SPT0348-E and -W in degrees
i_E = 63.0
i_W = 18.0

# import the PV fits data
filename='/Users/kevinlacaille/Documents/University/PhD/Research/SPT0348/Data/ALMA/Cycle_4/GalPaK/galpak_august20_2018/../../Band7/spt0348_band7_spw1_clean1000_contsub_2sig.image.pv_E_vel.fits'
hdul=fits.open(filename)
hdr=hdul[0].header
wcs=WCS(hdr)
data=hdul[0].data
hdul.close()

# guesstimate where the centre of the PV plot is on y axis
wcs.wcs.crpix[1]=40
wcs.wcs.crval[1]=0
temp=wcs.wcs.cdelt[1]
wcs.wcs.cdelt[1]=temp/1000



# Draw the figure. 
fig=pl.figure()
ax=pl.subplot(111,projection=wcs)
imgplot = ax.imshow(data*1000,cmap='magma_r',aspect='auto',origin='lower') #*1000 to turn to mJy/px
cbar = fig.colorbar(imgplot)
cbar.set_label(r'mJy pixel$^{-1}$')
ax.plot(arcsec_to_px(dx_E_gauss_arctan),vel_to_px(v_E_gauss_arctan),'o-',label='gauss arctan') 
ax.plot(arcsec_to_px(dx_E_exp_arctan),vel_to_px(v_E_exp_arctan),'o-',label='exp arctan', alpha=0.5) 
ax.set_xlabel('Offset (arcsec)')
ax.set_ylabel(r'Velocity (km s$^{-1}$)')

ax.set_ylim(vel_to_px(-500),vel_to_px(500))
ax.set_xlim(arcsec_to_px(-1.5),arcsec_to_px(1.5))
pl.text(arcsec_to_px(-1.35),vel_to_px(380),'SPT0348-E',size=16,color='k')
ax.legend()

pl.savefig('../Figures/GalPaK/PV_E_gauss_arctan.pdf',bbox_inches='tight')
pl.close()









# import the PV fits data
filename='/Users/kevinlacaille/Documents/University/PhD/Research/SPT0348/Data/ALMA/Cycle_4/GalPaK/galpak_august20_2018/../../Band7/spt0348_band7_spw1_clean1000_contsub_2sig.image.pv_W_vel.fits'
hdul=fits.open(filename)
hdr=hdul[0].header
wcs=WCS(hdr)
data=hdul[0].data
hdul.close()

# guesstimate where the centre of the PV plot is on y axis
wcs.wcs.crpix[1]=40
wcs.wcs.crpix[0]=17
wcs.wcs.crval[1]=0
temp=wcs.wcs.cdelt[1]
wcs.wcs.cdelt[1]=temp/1000

# Draw the figure. 
fig=pl.figure()
ax=pl.subplot(111,projection=wcs)
ax.imshow(data*1000,cmap='magma_r',aspect='auto',origin='lower')
ax.plot(arcsec_to_px(dx_W_gauss_arctan),vel_to_px(v_W_gauss_arctan),'o-',label='gauss arctan') 
ax.plot(arcsec_to_px(dx_W_exp_arctan),vel_to_px(v_W_exp_arctan),'o-',label='exp arctan', alpha=0.5) 
cbar = fig.colorbar(imgplot)
cbar.set_label(r'mJy pixel$^{-1}$')
ax.set_xlabel('Offset (arcsec)')
ax.set_ylabel(r'Velocity (km s$^{-1}$)')
ax.legend()
ax.set_ylim(vel_to_px(500),vel_to_px(-750))
ax.set_xlim(arcsec_to_px(-1.5),arcsec_to_px(1.5))
pl.text(arcsec_to_px(-1.35),vel_to_px(-550),'SPT0348-W',size=16)

pl.savefig('../Figures/GalPaK/PV_W_gauss_arctan.pdf',bbox_inches='tight')
pl.close()











################################
# Differences in rotation curves
################################
fig,ax1 = pl.subplots()

ax1.plot(dx_E_gauss_arctan,v_E_gauss_arctan,'--',label='gauss arctan',c='k')
ax1.plot(dx_E_exp_arctan,v_E_exp_arctan,'--',label='exp arctan',c='c')
ax1.set_xlabel('Offset (arcsec)')
ax1.set_ylabel(r'Velocity (km s$^{-1}$)')
ax1.legend(ncol=1,bbox_to_anchor=(-0.07, 1))

ax2 = ax1.twinx()

ax2.plot(dx_E_gauss_arctan,v_E_gauss_arctan - v_E_gauss_tanh,'-',label='gauss arctan - gauss tanh')
ax2.plot(dx_E_gauss_arctan,v_E_gauss_arctan - v_E_gauss_exp,'-',label='gauss arctan - gauss exp')
ax2.plot(dx_E_exp_arctan,v_E_gauss_arctan - v_E_exp_arctan,'-',label='gauss arctan - exp arctan')
ax2.plot(dx_E_exp_arctan,v_E_gauss_arctan - v_E_exp_tanh,'-',label='gauss arctan - exp tanh')
ax2.plot(dx_E_exp_arctan,v_E_gauss_arctan - v_E_exp_exp,'-',label='gauss arctan - exp exp')

ax2.text(-1.2,12,'SPT0348-E',size=16)
ax2.set_ylabel(r'$\Delta$ Velocity (km s$^{-1}$)')
ax2.legend(ncol=1,bbox_to_anchor=(1.55, 1))

pl.savefig('../Figures/GalPaK/Vrot_diff_E.pdf',bbox_inches='tight')
pl.close()





fig,ax1 = pl.subplots()

ax1.plot(dx_W_gauss_arctan,v_W_gauss_arctan,'--',label='gauss arctan',c='k')
ax1.plot(dx_W_exp_arctan,v_W_exp_arctan,'--',label='exp arctan',c='c')
ax1.set_xlabel('Offset (arcsec)')
ax1.set_ylabel(r'Velocity (km s$^{-1}$)')
ax1.legend(ncol=1,bbox_to_anchor=(-0.07, 1))

ax2 = ax1.twinx()

ax2.plot(dx_W_gauss_arctan,v_W_gauss_arctan - v_W_gauss_tanh,'-',label='gauss arctan - gauss tanh')
#ax2.plot(dx_W_gauss_arctan,v_W_gauss_arctan - v_W_gauss_exp,'-',label='gauss arctan - gauss exp')
ax2.plot(-1.5,10,'w.')
ax2.plot(dx_W_exp_arctan,v_W_gauss_arctan - v_W_exp_arctan,'-',label='gauss arctan - exp arctan')
ax2.plot(dx_W_exp_arctan,v_W_gauss_arctan - v_W_exp_tanh,'-',label='gauss arctan - exp tanh')
ax2.plot(dx_W_exp_arctan,v_W_gauss_arctan - v_W_exp_exp,'-',label='gauss arctan - exp exp')

ax2.text(-2,27,'SPT0348-W',size=16)
ax2.set_ylabel(r'$\Delta$ Velocity (km s$^{-1}$)')
ax2.legend(ncol=1,bbox_to_anchor=(1.55, 1))

pl.savefig('../Figures/GalPaK/Vrot_diff_W.pdf',bbox_inches='tight')
pl.close()



