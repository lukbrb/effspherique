from Sol_nonlin import rk4
from main import H0



"""On cherche ici la surdensité minimale telle que le temps d'effondrement soit égal à celui de l'âge de l'univers, 
dans un univers de Einstein - De Sitter"""

t_eff = 2 / (3 * H0) # L'âge de l'univers

def cond_init(d_i_min, d_i_max, tolerance, delta_eff, t_eff):
    # Chercher le temps d'effondrement pour d_i_min:
    sol1, temps_eff_max = rk4(d_i_min, delta_eff)
    if temps_eff_max < t_eff:
        print("Surdensité minimum trop grande")
    # Chercher le temps d'effondrement pour d_i_min:
    sol2, temps_eff_min = rk4(d_i_max, delta_eff)
    if t_eff < temps_eff_min:
        print("Surdensité maximum trop petite")

    inter = abs(d_i_max - d_i_min)
    while inter > tolerance:
        milieu = (d_i_min + d_i_max) / 2
        sol, temps_moy = rk4(milieu, delta_eff)
        if temps_moy > t_eff:
            d_i_min = milieu
        else:
            d_i_max = milieu
        inter = d_i_max - d_i_min
    surd_init = (d_i_max + d_i_min) / 2
    return surd_init


sur_init = cond_init(1e-4, 1e-2, 1e-7, 1e4, t_eff)
print("Surdensité initiale:", sur_init)
