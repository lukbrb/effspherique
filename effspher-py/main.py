import numpy as np

# Déclaration des variables
Mp = 3.08567758 * 1e22  # Valeur d'un mégaparsec
H0 = 7 * 1e4 / Mp  # Constante de Hubble Actuelle
G = 6.6742 * 1e-11  # Constante de gravitation
rho_crit = (3 * H0 ** 2) / (8 * np.pi * G)  # Densité critique
sig = 1  # Densité de l'univers
milliard_annee = 3600 * 24 * 365 * 1e9

# CONDITIONS INITIALES
# CONDITIONS INITIALES
# On veut ti très petit devant 1/H0, or 1/H0= 4.4e17s (=13 952 308 472 ans)
# On peut donc prendre ti=300000 ans par exemple
# Grille de temps

ti = 300000 * 365 * 24 * 3600
tf = 1 / H0
N = 5000
ttab = np.linspace(ti, tf, N)
dt = (tf - ti) / N  # DEFINITION DU PAS DE TEMPS
ttab1 = ttab / milliard_annee


def H(t):
    return 2 / (3 * t)


def a(t):
    return (1.5 * H0 * t) ** (2 / 3)


def rho_m(t):
    return (sig * rho_crit) / (a(t) ** 3)
