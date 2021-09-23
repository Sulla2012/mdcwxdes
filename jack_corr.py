import h5py
import treecorr
import astropy
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import dill as pk
import time

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
gold_rand_cat = treecorr.Catalog(ra = gold_rand_ras, dec = gold_rand_decs, ra_units = 'deg', dec_units = 'deg')
mdcw_rand_cat = treecorr.Catalog(ra = mdcw_rand_ras, dec = mdcw_rand_decs, ra_units = 'deg', dec_units = 'deg')

print(len(mdcw_zs))
#For each mdcws cluster, we compute xi with that cluster removed

try: 
    xi_dir = pk.load(open('./jack_xi_dir.dl', 'rb'))
    w_dir = pk.load(open('./jack_w_dir.dl', 'rb'))
except:
    xi_dir = {}
    w_dir = {}

zs = np.arange(0.75, 1.525, 0.025)

for i in range(len(mdcw_ras)):
    
    #Skip completed xi's
    if i in xi_dir.keys(): 
        print('Skipped as already calculated')
        continue

    #Skip clusters w/o z
    if np.isnan(mdcw_zs[i]): continue
    
    #Skip one cluster which is above highest z bin
    if mdcw_zs[i] > 1.5: continue

    #Figure out what z bin it's in
    j = (mdcw_zs[i]-0.75)//0.025
    j = int(j)
    print(zs[j],mdcw_zs[i], zs[j+1])
    gold_flags = np.where((gold_zs >= zs[j]) & (gold_zs < zs[j+1]))[0]
    gold_cat = treecorr.Catalog(ra = gold_ras[gold_flags], dec = gold_decs[gold_flags], ra_units = 'deg', dec_units = 'deg')
   
    mdcw_cat = treecorr.Catalog( ra = [mdcw_ras[i]], dec = [mdcw_decs[i]], ra_units = 'deg', dec_units = 'deg')

    ang_dia_dist = cosmo.angular_diameter_distance((zs[j+1]+zs[j])/2)
    ang_dia_dist *= u.radian**(-1)
    #min_sep = 0.2h^-1 Mpc/ang_dia_dist = radians

    min_sep = 0.2*h100**(-1)*u.Mpc
    max_sep = 60*h100**(-1)*u.Mpc

    min_sep = (min_sep/ang_dia_dist).to(u.arcmin).value
    max_sep = (max_sep/ang_dia_dist).to(u.arcmin).value
    
    path = '/global/cscratch1/sd/jorlo/mdcwxdes/'
    toc = time.time()
    #Calculate d1d2 with just the one mdcw point
    d1d2_jack = treecorr.NNCorrelation(min_sep=min_sep, max_sep=max_sep, nbins=25., sep_units='arcmin')
    d1d2_jack.process(gold_cat, mdcw_cat)

    #D1R2 is gold_data x mdcw_rands so no need to replicate that here: we're not actually changing anything
    
    #Calculate d2r1 with just the one point
    d2r1_jack = treecorr.NNCorrelation(min_sep=min_sep, max_sep=max_sep, nbins=25., sep_units='arcmin')
    d2r1_jack.process(mdcw_cat, gold_rand_cat)
    tic = time.time()

    print('Took ', int(tic-toc), ' seconds to calculate jackknife dd and dr')
    #And randoms are again the same

    #Load full data
    d1d2 = pk.load(open(path+'d1d2_{}_{}.dl'.format(zs[j], zs[j+1]), 'rb'))
    d1r2 = pk.load(open(path+'d1r2_{}_{}.dl'.format(zs[j], zs[j+1]), 'rb'))
    d2r1 = pk.load(open(path+'d2r1_{}_{}.dl'.format(zs[j], zs[j+1]), 'rb'))
    r1r2 = pk.load(open(path+'rr_{}_{}.dl'.format(zs[j], zs[j+1]), 'rb'))

    #Subtract pairs due to this one mdcw cluster
    d1d2.npairs = d1d2.npairs - d1d2_jack.npairs
    d2r1.npairs = d2r1.npairs - d2r1_jack.npairs

    xi = 0
    N_cl = 0
    sigmabar_g = 0
    for k, z in enumerate(zs):
        #Somewhat inelligant way of stopping at the last bin
        if k == len(zs)-1: continue
        gold_flags = np.where((gold_zs >= zs[k]) & (gold_zs < zs[k+1]))[0]
        mdcw_flags = np.where((mdcw_zs >= zs[k]) & (mdcw_zs < zs[k+1]))[0]
        if k == j:
            xi_curr, varxi = d1d2.calculateXi(r1r2, d1r2, d2r1)
            xi += xi_curr * (len(mdcw_flags)-1)
            N_cl += len(mdcw_flags)-1
            sigmabar_g += len(gold_flags)/4143.2*(len(mdcw_flags)-1)
        else:
            N_cl += len(mdcw_flags)
            sigmabar_g += len(gold_flags)/4143.2*len(mdcw_flags)
            try:
                #Load the existing dd, dr, and rr per redshift bin. Note for some redshift bins
                #There are no clusters, so there is no xi. This it the reason for the try
                #except statement
                d1d2_curr = pk.load(open(path+'d1d2_{}_{}.dl'.format(zs[k], zs[k+1]), 'rb'))
                d1r2_curr = pk.load(open(path+'d1r2_{}_{}.dl'.format(zs[k], zs[k+1]), 'rb'))
                d2r1_curr = pk.load(open(path+'d2r1_{}_{}.dl'.format(zs[k], zs[k+1]), 'rb'))
                r1r2_curr = pk.load(open(path+'rr_{}_{}.dl'.format(zs[k], zs[k+1]), 'rb'))
            except:
                print('Skipping bin ', zs[k], zs[k+1])
                continue

            xi_curr, varxi = d1d2_curr.calculateXi(r1r2_curr, d1r2_curr, d2r1_curr)

            xi += xi_curr * len(mdcw_flags)
    print(N_cl)
    w_dir[i] = xi[:]
    print('w: ', w_dir[i])
    xi_new = xi / N_cl
    sigmabar_g /= N_cl
    print(sigmabar_g)
    xi_new *= sigmabar_g
    print('xi: ', xi_new)
    print('w: ' ,w_dir[i])
    xi_dir[i] = xi_new

    pk.dump(xi_dir, open('jack_xi_dir.dl', 'wb'))
    pk.dump(w_dir, open('jack_w_dir.dl', 'wb'))

