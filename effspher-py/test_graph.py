import matplotlib.pyplot as plt
import scipy.integrate as itg
import numpy as np
import os
from main import *
from Sol_lin import *
from Sol_nonlin import *


os.chdir("resultats")
# tests sur l'équation linéaire
plt.figure()
plt.title("Tests sur l'équation linéaire")
euler = Euler_lin(1e-3, 1e4)

w = rk4_lin(1e-3, 1e4)
z = rk2_lin(1e-3, 1e4)
c_i = [1e-3, H(ti) * 1e-3]
y = itg.odeint(equation, c_i, ttab)
plt.plot(ttab, euler[:, 0], label='Euler')
plt.plot(ttab, y[:, 0], label='Scipy')
plt.plot(ttab, z[:, 0], label='Runge-Kutta 2')
plt.plot(ttab, w[:, 0], label='Runge-Kutta 4')
plt.plot(ttab, ttab ** (2 / 3) * (1 / 4.4) * 1e-11, label=r'$t^{2/3}$')
plt.xlabel('Temps (en secondes)')
plt.ylabel('Valeur de la surdensité')
plt.legend()

# Temps d'effondrement en fonction de la surdensité initiale
t_eff = np.genfromtxt("temps_eff.txt")
for i in range(0, len(t_eff)):
    t_eff[i] = float(t_eff[i]) / milliard_annee
surdens_init = np.genfromtxt("surdensité_initiale.txt")
plt.figure()
plt.plot(surdens_init, t_eff, '+')
plt.ylabel("Temps d'effondrement (en milliards d'années)")
plt.xlabel("Surdensité initiale")
plt.yscale('log')
plt.title("Temps d'effondrement en fonction de la surdensité initiale")

plt.figure()
plt.title("Temps d'effondrement calculé avec RK4")
x, teff = rk4(4 * surd_mini, 3 * 1e4)
plt.plot(ttab1, x[:, 0], label="Surdensité initiale =" + str(round(4 * surd_mini, 2)))
plt.xlabel("Temps (en milliards d'années")
plt.ylabel("Valeur de la surdensité")
plt.xscale('log')
plt.yscale('log')
plt.legend()

# Temps d'effondrement en fonction de la surdensité d'effondrement

t_deltaeff = np.genfromtxt("temps_deltaeff.txt")
for i in range(0, len(t_deltaeff)):
    t_deltaeff[i] = float(t_deltaeff[i]) / milliard_annee
surdens_eff = np.genfromtxt("Delta_eff.txt")
effmax = t_deltaeff.argmax()
tmax = t_deltaeff[effmax]
plt.figure()
plt.plot(surdens_eff, t_deltaeff)
plt.title("Temps d'effondrement en fonction de la surdensité d'effondrement")
plt.plot(surdens_eff[effmax], tmax, "x", label="Temps d'effondrement constant à partir"
                                               " de $\delta_{eff}$ =" + str(round(surdens_eff[effmax])))
plt.xlabel("Valeur de la surdensité d'effondrement")
plt.ylabel("Temps d'effondrement (milliards d'années)")
plt.xscale('log')
plt.legend()

# Temps d'effondrement en fonction du pas de temps
temps_fctN = np.genfromtxt("temps_fctN.txt")
N = np.genfromtxt("pas_temps.txt")
pas_temps = np.zeros(len(N))
for i in range(0, len(N)):
    pas_temps[i] = float((tf - ti) / N[i])

# plt.figure()
# plt.title("Influence du pas de temps sur le temps d'effondrement")
# plt.xlabel("Valeur du pas de temps (milliards d'années)")
# plt.ylabel("Temps d'effondrement (milliards d'années)")
# plt.xscale('log')
# plt.plot(pas_temps / milliard_annee, temps_fctN, "+")
# plt.show()
