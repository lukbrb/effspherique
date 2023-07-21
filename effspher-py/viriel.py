import numpy as np
from cosmofunc import rho_m, G, a, ttab1, ttab, H, dt, milliard_annee
import nonlinear
import matplotlib.pyplot as plt

x, teff = nonlinear.rk4(4 * nonlinear.surd_mini, 3 * 1e4)
Masse = 1e16  # N'influe pas sur le résultat final de la surdensité
delta = x[:, 0]
dv = x[:, 1]


# Fonction pour déterminer t lorsque Rviriel est atteint

def t_viriel(borne_min, borne_max, tolerance, r, r_vir):
    for i in r[borne_min:borne_max]:
        inter = abs(i - r_vir)
        t = ttab[borne_min]
        while inter > tolerance:
            milieu = (borne_min + borne_max) / 2
            if i > r_vir:
                borne_min = milieu
            else:
                borne_max = milieu
            inter = borne_max - borne_min
            t += 1
        return t


# Fonction qui fait évoluer la surdensité après la surdensité viriel

# TODO: Vectoriser cette fonction
def surdensite(index):
    for i in range(0, index):
        delta[i] = delta[i]
    for i in range(index, len(delta)):
        delta[i] = rho_vir / rho_m(ttab[i]) - 1
    return delta


# Vitesse calculée analytiquement
# TODO: Probablement vectorisable
def vitesse(r_final):
    r = r_final
    _v = np.zeros(len(r))
    v_test = np.zeros(len(r))
    for i in range(0, len(r)):
        t = ttab[i]
        _v[i] = H(t) * r[i] * (1 - (dv[i] * r[i] ** 3) / (9 * Masse * G * t))
        v_test[i] = H(t) * r[i] - (r[i] * dv[i]) / (3 * (1 + delta[i]))
    return _v, v_test


# def Vitesse_alter(R):  # VITESSE CALCULEE DIRECTEMENT A PARTIR DE R ET dt
#     V = np.zeros(len(R))
#     for i in range(0, len(R) - 1):
#         dR = R[i + 1] - R[i]
#         V[i] = dR / dt
#     return V
#
def Vitesse_alter(r):  # VITESSE CALCULEE DIRECTEMENT A PARTIR DE R ET dt
    dr = np.diff(r)
    _v = dr / dt
    return np.pad(_v, (0, 1), mode='constant')  # numpy pad pour avoir len(v) == len(r)


# Formule pour l'énergie cinétique
# def Ecin(Rcourt, vitesse):
#     R = Rcourt
#     Ecin = np.zeros(len(R))
#     for i in range(0, len(R)):
#         v = vitesse[i]
#         Ecin[i] = (3 * Masse * v ** 2) / 10
#     return Ecin

def Ecin(_v):
    return (3 * Masse * _v ** 2) / 10


# Calcul du rayon viriel en prenant la moitié de son rayon max (méthode 1)

R = ((3 * Masse) / (4 * np.pi * rho_m(ttab) * (1 + delta))) ** (1 / 3)
Rmax, imax = R.max(), R.argmax()
tmax = ttab[imax]
Rvir = Rmax * 0.5

# On ajuste les valeurs de la surdensité une fois R_vir atteint

tvir1 = t_viriel(377, 378, 1e-6, R, Rvir)
rho_vir = (3 * Masse) / (4 * np.pi * Rvir ** 3)  # Densité volumique de virialisation
delta_vir = rho_vir / rho_m(tvir1) - 1  # Surdensité viriel
delta[delta > delta_vir] = delta_vir  # Obsolète
indice = delta.argmax()
print("Indice de la virialisation :", indice)
print("Temps viriel calculé avec la fonction", tvir1)

# Rayon et surdensité finaux
delta_final = surdensite(indice)
R_final = ((3 * Masse) / (4 * np.pi * rho_m(ttab) * (1 + delta_final))) ** (1 / 3)

# Calcul des vitesses via méthodes numériques et formule analytique
v, vtest = vitesse(R_final)
Vitesse2 = Vitesse_alter(R_final)

# Calcul des énergies - On regarde où est-ce qu'elles se croisent pour déterminer le rayon viriel (méthode 2)
Epp = -(3 * G * Masse ** 2) / (5 * R_final)
Ec = Ecin(v)
Ectest = Ecin(Vitesse2)
Epvir = (3 * G * Masse ** 2) / (5 * Rvir)

surd_viriel_finale = 1 + delta_vir * (a(teff) / a(tvir1)) ** 3
print("Surdensité viriel finale :", surd_viriel_finale)
print(delta_final[410] + 1)  # Ici l'indice n'est valable que pour une certaine valeur de N (pour N=5000)

if __name__ == '__main__':
    plt.figure()
    plt.plot(ttab, 2 * Ectest, label="Energie cinétique en fonction de t")
    plt.plot(ttab, abs(Epp), label="Energie potentielle en fonction de t")
    plt.plot(tvir1, Epvir, "x", label='Energie potentielle au rayon viriel')

    plt.xscale('log')
    plt.yscale('log')
    plt.legend()

    plt.figure()
    # plt.plot(tmax,Rmax,'x',label="Rayon max")
    plt.plot(tvir1 / milliard_annee, Rvir, 'x', label="Rayon Viriel")
    plt.plot(tvir1 / milliard_annee, delta_vir, "+", label="Surdensité viriel")
    plt.plot(ttab1, R_final, label="Rayon")
    plt.plot(ttab1, delta_final, label="Surdensité")
    plt.plot(tmax / milliard_annee, Rmax, 'x', label="Rayon max")
    plt.xscale('log')
    plt.yscale('log')
    plt.legend()

    plt.figure()

    plt.plot(ttab, abs(v), label="Vitesse théorique")
    plt.plot(ttab, abs(Vitesse2), label="Vitesse numérique")
    plt.xscale('log')
    plt.yscale('log')
    plt.legend()
    plt.show()
