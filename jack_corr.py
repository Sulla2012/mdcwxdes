import h5py
import treecorr
import astropy
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import pickle as pk

from astropy.cosmology import Planck15 as cosmo
from astropy import constants as const
from astropy import units as u

h100 = cosmo.H0/(1.00E+02*u.km*u.s**(-1)*u.Mpc**(-1))
h  = const.h.value





print('load gold')
gold = fits.open('des_gold_dec_cut.fits')
gold_ras, gold_decs, gold_zs = np.array(gold[1].data['RADeg']), np.array(gold[1].data['decDeg']), np.array(gold[1].data['z'])
print('load gold rands')
gold_rand = fits.open('des_rand_dec_cut.fits')
gold_rand_ras, gold_rand_decs = np.array(gold_rand[1].data['RADeg']), np.array(gold_rand[1].data['decDeg'])
print('here')
mdcw = fits.open('MaDCoWS_DES_mass_v1.0.fits')
mdcw_ras, mdcw_decs, mdcw_zs = np.array(mdcw[1].data['RADeg']), np.array(mdcw[1].data['decDeg']), np.array(mdcw[1].data['redshift'])
dec_cut = np.where((mdcw_decs >-30))[0]
mdcw_ras, mdcw_decs, mdcw_zs = mdcw_ras[dec_cut], mdcw_decs[dec_cut], mdcw_zs[dec_cut]

mdcw_rand = fits.open('MaDCoWS_randoms_DES_mask.fits')
mdcw_rand_ras, mdcw_rand_decs = np.array(mdcw_rand[1].data['RADeg']), np.array(mdcw_rand[1].data['decDeg'])
dec_cut = np.where((mdcw_rand_decs >-30))[0]
mdcw_ras_rand, mdcw_decs_rand = mdcw_rand_ras[dec_cut], mdcw_rand_decs[dec_cut]

#Rands cat the same for all bins
print(' rand cats')
gold_rand_cat = treecorr.Catalog(ra = gold_rand_ras, dec = gold_rand_decs, ra_units = 'deg', dec_units = 'deg', npatch = 50)
mdcw_rand_cat = treecorr.Catalog(ra = mdcw_rand_ras, dec = mdcw_rand_decs, ra_units = 'deg', dec_units = 'deg',
        patch_centers=gold_rand_cat.patch_centers)


#For each mdcws cluster, we compute xi with that cluster removed
xi_dir = {}
zs = np.arange(0.75, 1.5, 0.025)
for i in range(len(mdcw_ras)):
    #Figure out what z bin it's in
    j = (mdcw_zs[i]-0.75)//0.025
    print(mdcw_zs[i], j)
    gold_flags = np.where((gold_zs >= zs[j]) & (gold_zs < zs[j+1]))[0]
    gold_cat = treecorr.Catalog(ra = gold_ras[gold_flags], dec = gold_decs[gold_flags], ra_units = 'deg', dec_units = 'deg',
            patch_centers=gold_rand_cat.patch_centers)


