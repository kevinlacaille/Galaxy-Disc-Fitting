
imagename='spt0348_band7_spw1_clean1000_contsub_2sig.image'
outfile='spt0348_band7_spw1_clean1000_contsub_2sig.image.pv_W_vel'
#mask='NGC5257_12CO21_combine_noise40.mask'
mode='length'
center=["03h48m42.090s","-62d20m51.195s"]
length = {"value": 3.4, "unit": "arcsec"}
pa={"value": -56.5, "unit": "deg"}

impv(imagename=imagename,outfile=outfile,mode=mode,center=center,length=length,pa=pa)#,mask=mask)
exportfits(imagename=outfile,fitsimage='spt0348_band7_spw1_clean1000_contsub_2sig.image.pv_W_vel.fits',velocity=True, dropdeg=True)

