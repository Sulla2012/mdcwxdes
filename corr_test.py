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

# Constants
# --------------------------------------------------------

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



zs = np.arange(0.75, 1.5, 0.025)
xi_dir = {}

for i in range(len(zs)-1):
    print(zs[i])
    gold_flags = np.where((gold_zs >= zs[i]) & (gold_zs < zs[i+1]))[0]
    gold_cat = treecorr.Catalog(ra = gold_ras[gold_flags], dec = gold_decs[gold_flags], ra_units = 'deg', dec_units = 'deg', 
            patch_centers=gold_rand_cat.patch_centers)

    mdcw_flags = np.where((mdcw_zs >= zs[i]) & (mdcw_zs < zs[i+1]))[0]
    if len(mdcw_flags) == 0:
        continue
    mdcw_cat = treecorr.Catalog(ra = mdcw_ras[mdcw_flags], dec = mdcw_decs[mdcw_flags], ra_units = 'deg', dec_units = 'deg', 
            patch_centers=gold_rand_cat.patch_centers)

    ang_dia_dist = cosmo.angular_diameter_distance((zs[i+1]+zs[i])/2)
    ang_dia_dist *= u.radian**(-1)
    #min_sep = 0.2h^-1 Mpc/ang_dia_dist = radians
    
    min_sep = 0.2*h100**(-1)*u.Mpc
    max_sep = 60*h100**(-1)*u.Mpc
    
    min_sep = (min_sep/ang_dia_dist).to(u.arcmin).value
    max_sep = (max_sep/ang_dia_dist).to(u.arcmin).value

    print(min_sep, max_sep)
    
    
    #According to the treecorr docs, jackkifing the randoms is not really neccesarry:
    #https://rmjarvis.github.io/TreeCorr/_build/html/cov.html
    #But we can add it later if needed, which is described above
    d1d2 = treecorr.NNCorrelation(min_sep=min_sep, max_sep=max_sep, nbins=25.,
                            sep_units='arcmin',var_method='jackknife')
    
    print('dd')
    d1d2.process(gold_cat, mdcw_cat)
    
    d1r2 = treecorr.NNCorrelation(min_sep=min_sep, max_sep=max_sep, nbins=25.,
                            sep_units='arcmin')
    print('d1r2')
    d1r2.process(gold_cat, mdcw_rand_cat)
    
    d2r1 = treecorr.NNCorrelation(min_sep=min_sep, max_sep=max_sep, nbins=25.,
                            sep_units='arcmin')
    print('d2r1')
    d2r1.process(mdcw_cat, gold_rand_cat)
    
    r1r2 = treecorr.NNCorrelation(min_sep=min_sep, max_sep=max_sep, nbins=25.,
                            sep_units='arcmin')
    print('r1r2')
    r1r2.process(gold_rand_cat, mdcw_rand_cat)
    
    print('dd ', d1d2.npairs)
    print('d1r2 ', d1r2.npairs)
    print('d2r1 ', d2r1.npairs)
    print('rr ', r1r2.npairs)

    xi_dir[zs[i]] = {}
    
    xi, varxi = d1d2.calculateXi(r1r2, d1r2, d2r1)

    cov = d1d2.cov
    
    xi_dir[zs[i]]['xi'] = xi
    xi_dir[zs[i]]['varxi'] = varxi
    xi_dir[zs[i]]['cluster_num'] = len(mdcw_flags)
    xi_dir[zs[i]]['galaxy_num'] = len(gold_flags)
    xi_dir[zs[i]]['cov'] = cov
    
    
pk.dump(xi_dir, open('xi_dir.pk', "wb"))

cluster_sum = 0
xi = 0
var_xi = 0

for key in xi_dir.keys():
    xi += xi_dir[key]['cluster_num']*xi_dir[key]['xi']
    cluster_sum += xi_dir[key]['cluster_num']
    var_xi += xi_dir[key]['varxi']**2
    
xi /= cluster_sum

r = np.exp(d1d2.meanlogr)
sig = np.sqrt(var_xi)

plt.plot(r, xi, color='blue')
plt.plot(r, -xi, color='blue', ls=':')
plt.errorbar(r[xi>0], xi[xi>0], yerr=sig[xi>0], color='blue', lw=0.1, ls='')
plt.errorbar(r[xi<0], -xi[xi<0], yerr=sig[xi<0], color='blue', lw=0.1, ls='')
leg = plt.errorbar(-r, xi, yerr=sig, color='blue')

plt.xscale('log')
plt.yscale('log', nonposy='clip')
plt.xlabel(r'$\theta$ (arcmin)')

plt.legend([leg], [r'$w(\theta)$'], loc='lower left')
plt.xlim([1,100])
plt.savefig('mdcwxdes_cor.pdf')
plt.show()

d1d2.write('d1d2.out')
d1r2.write('d1r2.out')
d2r1.write('d2r1.out')
r1r2.write('r1r2.out')

