import treecorr
import astropy
from astropy.io import fits
import numpy as np
import dill as dl
import time
import os.path

from astropy.cosmology import Planck15 as cosmo
from astropy import constants as const
from astropy import units as u

# Constants
# --------------------------------------------------------

h100 = cosmo.H0/(1.00E+02*u.km*u.s**(-1)*u.Mpc**(-1))
h  = const.h.value


#load gold rands
gold_rand = fits.open('des_rand_dec_cut.fits')
gold_rand_ras, gold_rand_decs = np.array(gold_rand[1].data['RADeg']), np.array(gold_rand[1].data['decDeg'])

#load mdcw rands

mdcw_rand = fits.open('MaDCoWS_randoms_DES_mask.fits')
mdcw_rand_ras, mdcw_rand_decs = np.array(mdcw_rand[1].data['RADeg']), np.array(mdcw_rand[1].data['decDeg'])
dec_cut = np.where((mdcw_rand_decs >-30))[0]
mdcw_ras_rand, mdcw_decs_rand = mdcw_rand_ras[dec_cut], mdcw_rand_decs[dec_cut]


#make randoms catalogs
print('Making randoms catalogs')
gold_rand_cat = treecorr.Catalog(ra = gold_rand_ras, dec = gold_rand_decs, ra_units = 'deg', dec_units = 'deg')
mdcw_rand_cat = treecorr.Catalog(ra = mdcw_rand_ras, dec = mdcw_rand_decs, ra_units = 'deg', dec_units = 'deg')

zs = np.arange(0.75, 1.5, 0.025)

print('starting correlation')
#Calculate rr for each z bin
for i in range(len(zs)-1):
    print(zs[i])
    fname = '/global/cscratch1/sd/jorlo/mdcwxdes/rr_{}_{}.dl'.format(zs[i], zs[i+1])
    if os.path.isfile(fname):
        print('rr already exists. If you were trying to remake this rr, please delete it')
        continue

    toc = time.time()
    
    #Calculate angular separation for z bin
    ang_dia_dist = cosmo.angular_diameter_distance((zs[i+1]+zs[i])/2)
    ang_dia_dist *= u.radian**(-1)
    #min_sep = 0.2h^-1 Mpc/ang_dia_dist = radians

    min_sep = 0.2*h100**(-1)*u.Mpc
    max_sep = 60*h100**(-1)*u.Mpc

    min_sep = (min_sep/ang_dia_dist).to(u.arcmin).value
    max_sep = (max_sep/ang_dia_dist).to(u.arcmin).value

    #Make rr object and process
    r1r2 = treecorr.NNCorrelation(min_sep=min_sep, max_sep=max_sep, nbins=25.,
                            sep_units='arcmin')
    r1r2.process(gold_rand_cat, mdcw_rand_cat)

    dl.dump(r1r2, open(fname, "wb"))


    tic = time.time()

    print('Took ', tic-toc, ' to calculate rr')




