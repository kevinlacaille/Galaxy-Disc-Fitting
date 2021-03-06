import numpy as np
import matplotlib.pyplot as pl

dir = '/Users/kevinlacaille/Documents/University/PhD/Research/SPT0348/Data/ALMA/Cycle_4/GalPaK/galpak_august20_2018/'

#0.1"/px * 6.875 kpc/" * 1000pc/kpc
# SPT0348-E data
# flux profile & rotation curve profile
# Gaussian flux profile
E_gauss_arctan_chain = np.loadtxt(dir + 'galpak_SPT0348_E_run1_gauss_arctan_chain.dat',unpack=True) #x,y,z,flux,radius,incl,pa,rV,Vmax,vel disp,reduced chi
rhalf_E_gauss_arctan = E_gauss_arctan_chain[4] * 0.1*6.875
incl_E_gauss_arctan = E_gauss_arctan_chain[5]
rV_E_gauss_arctan = E_gauss_arctan_chain[7] * 0.1*6.875*1000 #0.1"/px * 6.875 kpc/" * 1000pc/kpc
Vmax_E_gauss_arctan = E_gauss_arctan_chain[8]
sigma_E_gauss_arctan = E_gauss_arctan_chain[-2]
E_gauss_arctan_params = np.loadtxt(dir + 'galpak_SPT0348_E_run1_gauss_arctan_galaxy_parameters.dat',unpack=True) #x,y,z,flux,radius,incl,pa
E_gauss_arctan_Vrot = np.loadtxt(dir + 'galpak_SPT0348_E_run1_gauss_arctan_true_Vrot.dat',unpack=True) #dx,V,flux
#
E_gauss_tanh_chain = np.loadtxt(dir + 'galpak_SPT0348_E_run1_gauss_tanh_chain.dat',unpack=True)
rhalf_E_gauss_tanh = E_gauss_tanh_chain[4] * 0.1*6.875
incl_E_gauss_tanh = E_gauss_tanh_chain[5]
rV_E_gauss_tanh = E_gauss_tanh_chain[7] * 0.1*6.875*1000
Vmax_E_gauss_tanh = E_gauss_tanh_chain[8]
sigma_E_gauss_tanh = E_gauss_tanh_chain[-2]
E_gauss_tanh_params = np.loadtxt(dir + 'galpak_SPT0348_E_run1_gauss_tanh_galaxy_parameters.dat',unpack=True)
E_gauss_tanh_Vrot = np.loadtxt(dir + 'galpak_SPT0348_E_run1_gauss_tanh_true_Vrot.dat',unpack=True)
#
E_gauss_exp_chain = np.loadtxt(dir + 'galpak_SPT0348_E_run1_gauss_exp_chain.dat',unpack=True)
rhalf_E_gauss_exp = E_gauss_exp_chain[4] * 0.1*6.875
incl_E_gauss_exp = E_gauss_exp_chain[5]
rV_E_gauss_exp = E_gauss_exp_chain[7] * 0.1*6.875*1000
Vmax_E_gauss_exp = E_gauss_exp_chain[8]
sigma_E_gauss_exp = E_gauss_exp_chain[-2]
E_gauss_exp_params = np.loadtxt(dir + 'galpak_SPT0348_E_run1_gauss_exp_galaxy_parameters.dat',unpack=True)
E_gauss_exp_Vrot = np.loadtxt(dir + 'galpak_SPT0348_E_run1_gauss_exp_true_Vrot.dat',unpack=True)
# Exponential flux profile
E_exp_arctan_chain = np.loadtxt(dir + 'galpak_SPT0348_E_run1_exp_arctan_chain.dat',unpack=True)
rhalf_E_exp_arctan = E_exp_arctan_chain[4] * 0.1*6.875
incl_E_exp_arctan = E_exp_arctan_chain[5]
rV_E_exp_arctan = E_exp_arctan_chain[7] * 0.1*6.875*1000
Vmax_E_exp_arctan = E_exp_arctan_chain[8]
sigma_E_exp_arctan = E_exp_arctan_chain[-2]
E_exp_arctan_params = np.loadtxt(dir + 'galpak_SPT0348_E_run1_exp_arctan_galaxy_parameters.dat',unpack=True)
E_exp_arctan_Vrot = np.loadtxt(dir + 'galpak_SPT0348_E_run1_exp_arctan_true_Vrot.dat',unpack=True)
#
E_exp_tanh_chain = np.loadtxt(dir + 'galpak_SPT0348_E_run1_exp_tanh_chain.dat',unpack=True)
rhalf_E_exp_tanh = E_exp_tanh_chain[4] * 0.1*6.875
incl_E_exp_tanh = E_exp_tanh_chain[5]
rV_E_exp_tanh = E_exp_tanh_chain[7] * 0.1*6.875*1000
Vmax_E_exp_tanh = E_exp_tanh_chain[8]
sigma_E_exp_tanh = E_exp_tanh_chain[-2]
E_exp_tanh_params = np.loadtxt(dir + 'galpak_SPT0348_E_run1_exp_tanh_galaxy_parameters.dat',unpack=True)
E_exp_tanh_Vrot = np.loadtxt(dir + 'galpak_SPT0348_E_run1_exp_tanh_true_Vrot.dat',unpack=True)
#
E_exp_exp_chain = np.loadtxt(dir + 'galpak_SPT0348_E_run1_exp_exp_chain.dat',unpack=True)
rhalf_E_exp_exp = E_exp_exp_chain[4] * 0.1*6.875
incl_E_exp_exp = E_exp_exp_chain[5]
rV_E_exp_exp = E_exp_exp_chain[7] * 0.1*6.875*1000
Vmax_E_exp_exp = E_exp_exp_chain[8]
sigma_E_exp_exp = E_exp_exp_chain[-2]
E_exp_exp_params = np.loadtxt(dir + 'galpak_SPT0348_E_run1_exp_exp_galaxy_parameters.dat',unpack=True)
E_exp_exp_Vrot = np.loadtxt(dir + 'galpak_SPT0348_E_run1_exp_exp_true_Vrot.dat',unpack=True)












#################
# Incl. R_1/2 rV
#################

fig = pl.figure(figsize=(10,10))
pl.subplots_adjust(hspace=0,wspace=0)
pl.rc('font',size=15)
pl.rc('mathtext', default='regular')
#subplot(rows, columns, index)
ax1 = fig.add_subplot(6,3,1)
ax2 = fig.add_subplot(6,3,2)
ax3 = fig.add_subplot(6,3,3)

ax4 = fig.add_subplot(6,3,4)
ax5 = fig.add_subplot(6,3,5)
ax6 = fig.add_subplot(6,3,6)

ax7 = fig.add_subplot(6,3,7)
ax8 = fig.add_subplot(6,3,8)
ax9 = fig.add_subplot(6,3,9)

ax10 = fig.add_subplot(6,3,10)
ax11 = fig.add_subplot(6,3,11)
ax12 = fig.add_subplot(6,3,12)

ax13 = fig.add_subplot(6,3,13)
ax14 = fig.add_subplot(6,3,14)
ax15 = fig.add_subplot(6,3,15)

ax16 = fig.add_subplot(6,3,16)
ax17 = fig.add_subplot(6,3,17)
ax18 = fig.add_subplot(6,3,18)



bbox=dict(boxstyle="round", ec='k',fc='none',facecolor='none', edgecolor='blue')

########
# Incl.
########
ax1.hist(incl_E_gauss_arctan,20,alpha=0.4,hatch='//',label='gauss arctan')
ax1.axvline(np.median(incl_E_gauss_arctan),ls='--',c='k')
#ax4.text(62,3500,'gauss tanh',size=16)#,bbox=bbox)
ax1.set_xticklabels('')
ax1.set_yticklabels('')
#ax1.text(67,4900,'SPT0348-E',size=20)

ax4.hist(incl_E_gauss_tanh,20,alpha=0.4,hatch='//',label='gauss tanh')
ax4.axvline(np.median(incl_E_gauss_tanh),ls='--',c='k')
#ax4.text(62,3500,'gauss tanh',size=16)#,bbox=bbox)
ax4.set_xticklabels('')
ax4.set_yticklabels('')

ax7.hist(incl_E_gauss_exp,20,alpha=0.4,hatch='//',label='gauss exp')
ax7.axvline(np.median(incl_E_gauss_exp),ls='--',c='k')
#ax7.text(62,3500,'gauss exp',size=16)#,bbox=bbox)
ax7.set_xticklabels('')
ax7.set_yticklabels('')

ax10.hist(incl_E_exp_arctan,20,alpha=0.4,hatch='//',label='exp arctan')
ax10.axvline(np.median(incl_E_exp_arctan),ls='--',c='k')
#ax10.text(62,3500,'exp arctan', size=16)#,bbox=bbox)
ax10.set_xticklabels('')
ax10.set_yticklabels('')

ax13.hist(incl_E_exp_tanh,20,alpha=0.4,hatch='//',label='exp tanh')
ax13.axvline(np.median(incl_E_exp_tanh),ls='--',c='k')
#ax13.text(62,3500,'exp tanh', size=16)#,bbox=bbox)
ax13.set_xticklabels('')
ax13.set_yticklabels('')

ax16.hist(incl_E_exp_exp,20,alpha=0.4,hatch='//',label='exp exp')
ax16.axvline(np.median(incl_E_exp_exp),ls='--',c='k')
#ax16.text(62,3500,'exp exp', size=16)#,bbox=bbox)
ax16.set_xlabel('inclination (deg)')
ax16.set_ylabel('N')
#ax13.set_xticklabels(['170','172','174','176','178'])



########
# R_1/2
########
ax2.hist(rhalf_E_gauss_arctan,20,alpha=0.4,hatch='//',label='gauss arctan')
ax2.axvline(np.median(rhalf_E_gauss_arctan),ls='--',c='k')
ax2.set_yticklabels('')
ax2.set_xticklabels('')
ax2.text(5.17, 7000, 'SPT0348-E', size=20)

ax5.hist(rhalf_E_gauss_tanh,20,alpha=0.4,hatch='//',label='gauss tanh')
ax5.axvline(np.median(rhalf_E_gauss_tanh),ls='--',c='k')
ax5.set_yticklabels('')
ax5.set_xticklabels('')

ax8.hist(rhalf_E_gauss_exp,20,alpha=0.4,hatch='//',label='gauss exp')
ax8.axvline(np.median(rhalf_E_gauss_exp),ls='--',c='k')
ax8.set_yticklabels('')
ax8.set_xticklabels('')

ax11.hist(rhalf_E_exp_arctan,20,alpha=0.4,hatch='//',label='exp arctan')
ax11.axvline(np.median(rhalf_E_exp_arctan),ls='--',c='k')
ax11.set_yticklabels('')
ax11.set_xticklabels('')

ax14.hist(rhalf_E_exp_tanh,20,alpha=0.4,hatch='//',label='exp tanh')
ax14.axvline(np.median(rhalf_E_exp_tanh),ls='--',c='k')
ax14.set_yticklabels('')
ax14.set_xticklabels('')

ax17.hist(rhalf_E_exp_exp,20,alpha=0.4,hatch='//',label='exp exp')
ax17.axvline(np.median(rhalf_E_exp_exp),ls='--',c='k')
ax17.set_yticklabels('')
ax17.set_xlabel(r'$R_{1/2}$ (kpc)')
#ax17.set_xticklabels(['5.20','5.25','5.30'])


#####
# rV
#####
ax3.hist(rV_E_gauss_arctan,20,alpha=0.4,hatch='//',label='gauss arctan')
ax3.axvline(np.median(rV_E_gauss_arctan),ls='--',c='k')
ax3.set_yticklabels('')
ax3.set_xticklabels('')
ax3.set_xscale('log')
ax3.text(1.5e-2,3500,'gauss arctan',size=16)#,bbox=bbox)

ax6.hist(rV_E_gauss_tanh,20,alpha=0.4,hatch='//',label='gauss tanh')
ax6.axvline(np.median(rV_E_gauss_tanh),ls='--',c='k')
ax6.set_yticklabels('')
ax6.set_xticklabels('')
ax6.set_xscale('log')
ax6.text(1.5e-2,3500,'gauss tanh',size=16)#,bbox=bbox)

ax9.hist(rV_E_gauss_exp,20,alpha=0.4,hatch='//',label='gauss exp')
ax9.axvline(np.median(rV_E_gauss_exp),ls='--',c='k')
ax9.set_yticklabels('')
ax9.set_xticklabels('')
ax9.set_xscale('log')
ax9.text(1.5e-2,3500,'gauss exp',size=16)#,bbox=bbox)

ax12.hist(rV_E_exp_arctan,20,alpha=0.4,hatch='//',label='exp arctan')
ax12.axvline(np.median(rV_E_exp_arctan),ls='--',c='k')
ax12.set_yticklabels('')
ax12.set_xticklabels('')
ax12.set_xscale('log')
ax12.text(1.5e-2,3500,'exp arctan',size=16)#,bbox=bbox)

ax15.hist(rV_E_exp_tanh,20,alpha=0.4,hatch='//',label='exp tanh')
ax15.axvline(np.median(rV_E_exp_tanh),ls='--',c='k')
ax15.set_yticklabels('')
ax15.set_xticklabels('')
ax15.set_xscale('log')
ax15.text(1.5e-2,3500,'exp tanh',size=16)#,bbox=bbox)

ax18.hist(rV_E_exp_exp,20,alpha=0.4,hatch='//',label='exp exp')
ax18.axvline(np.median(rV_E_exp_exp),ls='--',c='k')
ax18.set_yticklabels('')
ax18.set_xlabel(r'turnover radius (pc)')
ax18.set_xscale('log')
ax18.text(1.5e-2,3500,'exp exp',size=16)#,bbox=bbox)




# incl xlim
ax1.set_xlim(62,67)
ax4.set_xlim(62,67)
ax7.set_xlim(62,67)
ax10.set_xlim(62,67)
ax13.set_xlim(62,67)
ax16.set_xlim(62,67)

# R_1/2 xlims
r_12_xlo = 5.1
r_12_xhi = 5.5
beam_majoraxis = 0.749*6.875 /2 #0.749"/beam * 6.875kpc/" = kpc/beam, /2 for radius
ax2.set_xlim(r_12_xlo, r_12_xhi) 
ax5.set_xlim(r_12_xlo, r_12_xhi)
ax8.set_xlim(r_12_xlo, r_12_xhi)
ax11.set_xlim(r_12_xlo, r_12_xhi)
ax14.set_xlim(r_12_xlo, r_12_xhi)
ax17.set_xlim(r_12_xlo, r_12_xhi)
ax2_2 = ax2.twiny()
ax2_2.set_xlim(r_12_xlo/beam_majoraxis,r_12_xhi/beam_majoraxis) # kpc / kpc/beam = beam
ax2_2.set_xlim(r_12_xlo/beam_majoraxis,r_12_xhi/beam_majoraxis)
ax2_2.set_xlabel('beams')

# rV xlims
rV_xlo = 1e-2
rV_xhi = 30.0
ax3.set_xlim(rV_xlo, rV_xhi)
ax6.set_xlim(rV_xlo, rV_xhi)
ax9.set_xlim(rV_xlo, rV_xhi)
ax12.set_xlim(rV_xlo, rV_xhi)
ax15.set_xlim(rV_xlo, rV_xhi)
ax18.set_xlim(rV_xlo, rV_xhi)
ax3_2 = ax3.twiny()
ax3_2.set_xscale('log')
ax3_2.set_xlim(rV_xlo/beam_majoraxis*1e-3,rV_xhi/beam_majoraxis*1e-3)
ax3_2.set_xlim(rV_xlo/beam_majoraxis*1e-3,rV_xhi/beam_majoraxis*1e-3)
ax3_2.set_xlabel('beams')

# All ylims to be the same
ax1.set_ylim(0,4500)
ax2.set_ylim(0,4500)
ax3.set_ylim(0,4500)
ax4.set_ylim(0,4500)
ax5.set_ylim(0,4500)
ax6.set_ylim(0,4500)
ax7.set_ylim(0,4500)
ax8.set_ylim(0,4500)
ax9.set_ylim(0,4500)
ax10.set_ylim(0,4500)
ax11.set_ylim(0,4500)
ax12.set_ylim(0,4500)
ax13.set_ylim(0,4500)
ax14.set_ylim(0,4500)
ax15.set_ylim(0,4500)
ax16.set_ylim(0,4500)
ax17.set_ylim(0,4500)
ax18.set_ylim(0,4500)

pl.draw()
pl.savefig('../Figures/GalPaK/SPT0348-E_incl_R12_rV_panels.pdf',bbox_inches='tight')
pl.close()





















'''





# Inclination angle
#ispace_E = np.linspace(62,64.5,20)
pl.hist(incl_E_gauss_arctan,20,alpha=0.4,hatch='//',label='gauss arctan',color='navy')
pl.hist(incl_E_gauss_tanh,20,alpha=0.4,hatch='o',label='gauss tanh',color='sienna')
pl.hist(incl_E_gauss_exp,20,alpha=0.4,hatch="\\",label='gauss exp',color='plum')
pl.hist(incl_E_exp_arctan,20,alpha=0.4,hatch="-",label='exp arctan',color='skyblue')
pl.hist(incl_E_exp_tanh,20,alpha=0.4,hatch="+",label='exp tanh',color='tan')
pl.hist(incl_E_exp_exp,20,alpha=0.4,hatch="*",label='exp exp',color='palegreen')
pl.xlabel('inclination (deg)')
pl.ylabel('N')
pl.legend(ncol=2,loc=2,fontsize=12)
pl.savefig('../Figures/GalPaK/SPT0348-E_incl.pdf',bbox_inches='tight')
pl.close()


# Turnover radius
fig, ax1 = pl.subplots()
#rVspace_E = np.logspace(np.log10(0.006),np.log10(21),20)
ax1.hist(rV_E_gauss_arctan,20,alpha=0.4,hatch='//',label='gauss arctan',color='navy')
ax1.hist(rV_E_gauss_tanh,20,alpha=0.4,hatch='o',label='gauss tanh',color='sienna')
ax1.hist(rV_E_gauss_exp,20,alpha=0.4,hatch="\\",label='gauss exp',color='plum')
ax1.hist(rV_E_exp_arctan,20,alpha=0.4,hatch="-",label='exp arctan',color='skyblue')
ax1.hist(rV_E_exp_tanh,20,alpha=0.4,hatch="+",label='exp tanh',color='tan')
ax1.hist(rV_E_exp_exp,20,alpha=0.4,hatch="*",label='exp exp',color='palegreen')

ax2 = ax1.twiny()
ax1.set_xlim(1e-2, 30)
ax1.set_xlabel('turnover radius (pc)')
ax2.set_xlim(1e-2/0.749,30.0/0.749)
ax2.set_xlabel('beams')

ax1.set_xscale('log') 
ax2.set_xscale('log') 
ax1.set_ylabel('N')
ax1.legend(ncol=2)
pl.savefig('../Figures/GalPaK/SPT0348-E_rV.pdf',bbox_inches='tight')
pl.close()

# Vmax
#Vspace_E = np.linspace(170,185,20)
pl.hist(Vmax_E_gauss_arctan,20,alpha=0.4,hatch='//',label='gauss arctan',color='navy')
pl.hist(Vmax_E_gauss_tanh,20,alpha=0.4,hatch='o',label='gauss tanh',color='sienna')
pl.hist(Vmax_E_gauss_exp,20,alpha=0.4,hatch="\\",label='gauss exp',color='plum')
pl.hist(Vmax_E_exp_arctan,20,alpha=0.4,hatch="-",label='exp arctan',color='skyblue')
pl.hist(Vmax_E_exp_tanh,20,alpha=0.4,hatch="+",label='exp tanh',color='tan')
pl.hist(Vmax_E_exp_exp,20,alpha=0.4,hatch="*",label='exp exp',color='palegreen')
#pl.xscale('log') 
pl.xlabel(r'V$_{max}$ (km s$^{-1}$)')
pl.ylabel('N')
pl.legend(loc=2)
#pl.xlim(165,185)
pl.savefig('../Figures/GalPaK/SPT0348-E_Vmax.pdf',bbox_inches='tight')
pl.close()

# sigma
#Vspace_E = np.linspace(170,185,20)
pl.hist(sigma_E_gauss_arctan,20,alpha=0.4,hatch='//',label='gauss arctan',color='navy')
pl.hist(sigma_E_gauss_tanh,20,alpha=0.4,hatch='o',label='gauss tanh',color='sienna')
pl.hist(sigma_E_gauss_exp,20,alpha=0.4,hatch="\\",label='gauss exp',color='plum')
pl.hist(sigma_E_exp_arctan,20,alpha=0.4,hatch="-",label='exp arctan',color='skyblue')
pl.hist(sigma_E_exp_tanh,20,alpha=0.4,hatch="+",label='exp tanh',color='tan')
pl.hist(sigma_E_exp_exp,20,alpha=0.4,hatch="*",label='exp exp',color='palegreen')
#pl.xscale('log') 
pl.xlabel(r'$\sigma_{max}$ (km s$^{-1}$)')
pl.ylabel('N')
pl.legend(loc=2)
#pl.xlim(165,185)
pl.savefig('../Figures/GalPaK/SPT0348-E_sigma.pdf',bbox_inches='tight')
pl.close()


# R_1/2
#Vspace_E = np.linspace(170,185,20)
fig, ax1 = pl.subplots()
ax1.hist(rhalf_E_gauss_arctan,20,alpha=0.4,hatch='//',label='gauss arctan',color='navy')
ax1.hist(rhalf_E_gauss_tanh,20,alpha=0.4,hatch='o',label='gauss tanh',color='sienna')
ax1.hist(rhalf_E_gauss_exp,20,alpha=0.4,hatch="\\",label='gauss exp',color='plum')
ax1.hist(rhalf_E_exp_arctan,20,alpha=0.4,hatch="-",label='exp arctan',color='skyblue')
ax1.hist(rhalf_E_exp_tanh,20,alpha=0.4,hatch="+",label='exp tanh',color='tan')
ax1.hist(rhalf_E_exp_exp,20,alpha=0.4,hatch="*",label='exp exp',color='palegreen')

ax2 = ax1.twiny()
ax1.set_xlim(5.1,5.5)
ax1.set_xlabel(r'R$_{1/2}$ (kpc)')
ax2.set_xlim(5.1/0.749,5.5/0.749)
ax2.set_xlabel('beams')

#ax1.set_xscale('log') 

#pl.xlabel(r'R$_{1/2}$ (kpc)')
ax1.set_ylabel('N')
ax1.legend(loc=1,ncol=2)
#pl.xlim(165,185)
pl.savefig('../Figures/GalPaK/SPT0348-E_Rhalf.pdf',bbox_inches='tight')
pl.close()















#0.1"/px * 6.875 kpc/" * 1000pc/kpc
# SPT0348-W data
# flux profile & rotation curve profile
# Gaussian flux profile
W_gauss_arctan_chain = np.loadtxt(dir + 'galpak_SPT0348_W_run1_gauss_arctan_chain.dat',unpack=True) #x,y,z,flux,radius,incl,pa,rV,Vmax,vel disp,reduced chi
incl_W_gauss_arctan = W_gauss_arctan_chain[5]
rV_W_gauss_arctan = W_gauss_arctan_chain[7] #* 0.1*6.875*1000 #0.1"/px * 6.875 kpc/" * 1000pc/kpc
Vmax_W_gauss_arctan = W_gauss_arctan_chain[8]
W_gauss_arctan_params = np.loadtxt(dir + 'galpak_SPT0348_W_run1_gauss_arctan_galaxy_parameters.dat',unpack=True) #x,y,z,flux,radius,incl,pa
W_gauss_arctan_Vrot = np.loadtxt(dir + 'galpak_SPT0348_W_run1_gauss_arctan_true_Vrot.dat',unpack=True) #dx,V,flux
#
W_gauss_tanh_chain = np.loadtxt(dir + 'galpak_SPT0348_W_run1_gauss_tanh_chain.dat',unpack=True)
incl_W_gauss_tanh = W_gauss_tanh_chain[5]
rV_W_gauss_tanh = W_gauss_tanh_chain[7] #* 0.1*6.875*1000
Vmax_W_gauss_tanh = W_gauss_tanh_chain[8]
W_gauss_tanh_params = np.loadtxt(dir + 'galpak_SPT0348_W_run1_gauss_tanh_galaxy_parameters.dat',unpack=True)
W_gauss_tanh_Vrot = np.loadtxt(dir + 'galpak_SPT0348_W_run1_gauss_tanh_true_Vrot.dat',unpack=True)
#

# Exponential flux profile
W_exp_arctan_chain = np.loadtxt(dir + 'galpak_SPT0348_W_run1_exp_arctan_chain.dat',unpack=True)
incl_W_exp_arctan = W_exp_arctan_chain[5]
rV_W_exp_arctan = W_exp_arctan_chain[7] #* 0.1*6.875*1000
Vmax_W_exp_arctan = W_exp_arctan_chain[8]
W_exp_arctan_params = np.loadtxt(dir + 'galpak_SPT0348_W_run1_exp_arctan_galaxy_parameters.dat',unpack=True)
W_exp_arctan_Vrot = np.loadtxt(dir + 'galpak_SPT0348_W_run1_exp_arctan_true_Vrot.dat',unpack=True)
#
W_exp_tanh_chain = np.loadtxt(dir + 'galpak_SPT0348_W_run1_exp_tanh_chain.dat',unpack=True)
incl_W_exp_tanh = W_exp_tanh_chain[5]
rV_W_exp_tanh = W_exp_tanh_chain[7] #* 0.1*6.875*1000
Vmax_W_exp_tanh = W_exp_tanh_chain[8]
W_exp_tanh_params = np.loadtxt(dir + 'galpak_SPT0348_W_run1_exp_tanh_galaxy_parameters.dat',unpack=True)
W_exp_tanh_Vrot = np.loadtxt(dir + 'galpak_SPT0348_W_run1_exp_tanh_true_Vrot.dat',unpack=True)
#
W_exp_exp_chain = np.loadtxt(dir + 'galpak_SPT0348_W_run1_exp_exp_chain.dat',unpack=True)
incl_W_exp_exp = W_exp_exp_chain[5]
rV_W_exp_exp = W_exp_exp_chain[7] #* 0.1*6.875*1000
Vmax_W_exp_exp = W_exp_exp_chain[8]
W_exp_exp_params = np.loadtxt(dir + 'galpak_SPT0348_W_run1_exp_exp_galaxy_parameters.dat',unpack=True)
W_exp_exp_Vrot = np.loadtxt(dir + 'galpak_SPT0348_W_run1_exp_exp_true_Vrot.dat',unpack=True)



# Inclination angle
#ispace_E = np.linspace(62,64.5,20)
pl.hist(incl_W_gauss_arctan,20,alpha=0.4,hatch='//',label='gauss arctan',color='navy')
pl.hist(incl_W_gauss_tanh,20,alpha=0.4,hatch='o',label='gauss tanh',color='sienna')
#pl.hist(incl_W_gauss_exp,20,alpha=0.4,hatch="\\",label='gauss exp',color='plum')
pl.hist(incl_W_exp_arctan,20,alpha=0.4,hatch="-",label='exp arctan',color='skyblue')
pl.hist(incl_W_exp_tanh,20,alpha=0.4,hatch="+",label='exp tanh',color='tan')
pl.hist(incl_W_exp_exp,20,alpha=0.4,hatch="*",label='exp exp',color='palegreen')
pl.xlabel('inclination (deg)')
pl.ylabel('N')
pl.legend(ncol=2,loc=1,fontsize=12)
pl.savefig('../Figures/GalPaK/SPT0348-W_incl.pdf',bbox_inches='tight')
pl.close()

# Turnover radius
rVspace_W = np.logspace(np.log10(0.01),np.log10(1050),20)
pl.hist(rV_W_gauss_arctan,rVspace_W,alpha=0.4,hatch='//',label='gauss arctan',color='navy')
pl.hist(rV_W_gauss_tanh,rVspace_W,alpha=0.4,hatch='o',label='gauss tanh',color='sienna')
#pl.hist(rV_W_gauss_exp,rVspace_W,alpha=0.4,hatch="\\",label='gauss exp',color='plum')
pl.hist(rV_W_exp_arctan,rVspace_W,alpha=0.4,hatch="-",label='exp arctan',color='skyblue')
pl.hist(rV_W_exp_tanh,rVspace_W,alpha=0.4,hatch="+",label='exp tanh',color='tan')
pl.hist(rV_W_exp_exp,rVspace_W,alpha=0.4,hatch="*",label='exp exp',color='palegreen')
pl.xscale('log')
pl.yscale('log') 
pl.xlabel('turnover radius (pc)')
pl.ylabel('N')
pl.legend()
pl.savefig('../Figures/GalPaK/SPT0348-W_rV.pdf',bbox_inches='tight')
pl.close()

# Vmax
Vspace_W = np.linspace(950,1390,20)
pl.hist(Vmax_W_gauss_arctan,Vspace_W,alpha=0.4,hatch='//',label='gauss arctan',color='navy')
pl.hist(Vmax_W_gauss_tanh,Vspace_W,alpha=0.4,hatch='o',label='gauss tanh',color='sienna')
#pl.hist(Vmax_W_gauss_exp,Vspace_W,alpha=0.4,hatch="\\",label='gauss exp',color='plum')
pl.hist(Vmax_W_exp_arctan,Vspace_W,alpha=0.4,hatch="-",label='exp arctan',color='skyblue')
pl.hist(Vmax_W_exp_tanh,Vspace_W,alpha=0.4,hatch="+",label='exp tanh',color='tan')
pl.hist(Vmax_W_exp_exp,Vspace_W,alpha=0.4,hatch="*",label='exp exp',color='palegreen')
#pl.xscale('log') 
pl.xlabel(r'V$_{max}$ (km s$^{-1}$)')
pl.ylabel('N')
pl.legend(loc=2)
#pl.xlim(165,185)
pl.savefig('../Figures/GalPaK/SPT0348-W_Vmax.pdf',bbox_inches='tight')
pl.close()


'''