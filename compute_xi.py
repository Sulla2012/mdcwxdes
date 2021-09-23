import dill as pk
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
area_of_des = 4143.2

#Load xis
xi_dir = pk.load(open('./jack_xi_dir.dl', 'rb'))

#Load mdcws and gold catalogs
print('load gold')
gold = fits.open('des_gold_dec_cut.fits')
gold_ras, gold_decs, gold_zs = np.array(gold[1].data['RADeg']), np.array(gold[1].data['decDeg']), np.array(gold[1].data['z'])
print('load mdcw')
mdcw = fits.open('MaDCoWS_DES_mass_v1.0.fits')
mdcw_ras, mdcw_decs, mdcw_zs = np.array(mdcw[1].data['RADeg']), np.array(mdcw[1].data['decDeg']), np.array(mdcw[1].data['redshift'])
dec_cut = np.where((mdcw_decs >-30))[0]
mdcw_ras, mdcw_decs, mdcw_zs = mdcw_ras[dec_cut], mdcw_decs[dec_cut], mdcw_zs[dec_cut]

#First we need to compute the total number of clusters, which we did not divide by in jack_corr, but which we should (see eqn 25)
zs = np.arange(0.75, 1.525, 0.025)
N_cl = 0
for k in range(len(zs)-1):
    mdcw_flags = np.where((mdcw_zs >= zs[k]) & (mdcw_zs < zs[k+1]))[0]
    N_cl += len(mdcw_flags)

#Now compute sigma^bar_g, eqn 27
sigmabar_g = 0

for k in range(len(zs)-1):
    gold_flags = np.where((gold_zs >= zs[k]) & (gold_zs < zs[k+1]))[0]
    mdcw_flags = np.where((mdcw_zs >= zs[k]) & (mdcw_zs < zs[k+1]))[0]
    sigmabar_g += len(mdcw_flags) * len(gold_flags)/area_of_des
sigmabar_g /= N_cl    


temp = np.zeros((len(xi_dir.keys()), len(xi_dir[2])))
for i, key in enumerate(xi_dir.keys()):
    #Note in jack_corr we already multiplied w by N_cl,i, so here we divide by N_cl and multiply by sigma^bar_g
    temp[i] = xi_dir[key]/N_cl*sigmabar_g

print(temp.shape)
print(np.mean(temp, axis = 0))
print(np.std(temp, axis = 0))
xi = np.mean(temp, axis = 0)
r = np.logspace(np.log10(0.2), np.log10(60), 25)
plt.plot(r, xi, color='blue')
plt.plot(r, -xi, color='blue', ls=':')

plt.xscale('log')
plt.yscale('log', nonposy='clip')
plt.xlabel(r'$\theta$ (arcmin)')



plt.xlim([1,100])

plt.xscale('log')
plt.yscale('log')
plt.savefig('corrected_xi.pdf')

plt.close()

plt.figure()
for key in xi_dir.keys():
    plt.plot(r, xi_dir[key], alpha = 0.25, c = 'b')
plt.xscale('log')
plt.yscale('log', nonposy='clip')
plt.xlabel(r'$\theta$ (arcmin)')



plt.xlim([1,100])
plt.savefig('allxis.pdf')
