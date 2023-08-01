import matplotlib.pyplot as plt
import numpy as np

path = 'effspher-py/resultats/'
deff, ttab = np.loadtxt(path + 'eff.csv').T
moy_deff = 1.001 * np.mean(ttab)
idx = np.abs(ttab - moy_deff).argmin()
deff_seuil = deff[idx]
plt.figure()
plt.title(r"$t_{eff}$ en fonction de $\delta_{eff}$")
plt.plot(deff, ttab)
plt.plot(deff_seuil, ttab[idx], 'x', label=f'({deff_seuil:.2f}, {ttab[idx]: .2f})')
plt.hlines(moy_deff, deff[0], deff[-1], colors='red', linestyles='--')
plt.xlabel(r'$\delta_{eff}$')
plt.ylabel('Temps (Gyr)')
plt.yscale('log')
plt.legend()
plt.show()
