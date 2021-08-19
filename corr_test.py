import h5py
import treecorr
import astropy
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

#filename = "/global/cfs/cdirs/des/data_actxdes/des_data/save_gold.h5"

#full_cat = h5py.File(filename,'r')

#des_ras, des_decs = np.array(full_cat['gold_ra']), np.array(full_cat['gold_dec'])
#dec_cut = np.where((des_decs >-30))[0]
#des_ras, des_decs = des_ras[dec_cut], des_decs[dec_cut]

#des_rand_ras, des_rand_decs = np.array(full_cat['rand_ra']), np.array(full_cat['rand_dec'])
#dec_cut = np.where((des_rand_decs >-30))[0]
#des_rand_ras, des_rand_decs = des_rand_ras[dec_cut], des_rand_decs[dec_cut]

gold = fits.open('des_gold_dec_cut.fits')
gold_ras, gold_decs, gold_zs = np.array(gold[1].data['RADeg']), np.array(gold[1].data['decDeg']), np.array(gold[1].data['z'])

gold_rand = fits.open('des_rand_dec_cut.fits')
gold_rand_ras, gold_rand_decs, gold_rand_zs = np.array(gold_rand[1].data['RADeg']), np.array(gold_rand[1].data['decDeg']), np.array(gold_rand[1].data['z'])

#cat_lim = int(1e8)
#flags = np.random.randint(len(des_ras), size = cat_lim)
gold_cat = treecorr.Catalog(ra = gold_ras, dec = gold_decs, ra_units = 'deg', dec_units = 'deg')
#flags = np.random.randint(len(des_rand_ras), size = cat_lim)
gold_rand_cat = treecorr.Catalog(ra = gold_rand_ras, dec = gold_rand_decs, ra_units = 'deg', dec_units = 'deg')

mdcw = fits.open('MaDCoWS_DES_mass_v1.0.fits')
mdcw_ras, mdcw_decs = np.array(mdcw[1].data['RADeg']), np.array(mdcw[1].data['decDeg'])
dec_cut = np.where((mdcw_decs >-30))[0]
mdcw_ras, mdcw_decs = mdcw_ras[dec_cut], mdcw_decs[dec_cut]
mdcw_cat = treecorr.Catalog(ra = mdcw_ras, dec = mdcw_decs, ra_units = 'deg', dec_units = 'deg')

mdcw_rand = fits.open('MaDCoWS_randoms_DES_mask.fits')
mdcw_rand_ras, mdcw_rand_decs = np.array(mdcw_rand[1].data['RADeg']), np.array(mdcw_rand[1].data['decDeg'])
dec_cut = np.where((mdcw_rand_decs >-30))[0]
mdcw_ras, mdcw_decs = mdcw_rand_ras[dec_cut], mdcw_rand_decs[dec_cut]
mdcw_rand_cat = treecorr.Catalog(ra = mdcw_rand_ras, dec = mdcw_rand_decs, ra_units = 'deg', dec_units = 'deg')

d1d2 = treecorr.NNCorrelation(min_sep=1., max_sep=100., nbins=25.,
                            sep_units='arcmin')

d1d2.process(gold_cat, mdcw_cat)

d1r2 = treecorr.NNCorrelation(min_sep=1., max_sep=100., nbins=25.,
                            sep_units='arcmin')

d1r2.process(gold_cat, mdcw_rand_cat)

d2r1 = treecorr.NNCorrelation(min_sep=1., max_sep=100., nbins=25.,
                            sep_units='arcmin')

d2r1.process(mdcw_cat, gold_rand_cat)

r1r2 = treecorr.NNCorrelation(min_sep=1., max_sep=100., nbins=25.,
                            sep_units='arcmin')

r1r2.process(gold_rand_cat, mdcw_rand_cat)

print('dd ', d1d2.npairs)
print('d1r2 ', d1r2.npairs)
print('d2r1 ', d2r1.npairs)
print('rr ', r1r2.npairs)

xi, varxi = d1d2.calculateXi(r1r2, d1r2, d2r1)

r = np.exp(d1d2.meanlogr)
sig = np.sqrt(varxi)

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

