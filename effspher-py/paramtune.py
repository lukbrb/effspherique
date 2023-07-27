"""
Sont regroupées dans ce fichier les différentes fonctions permettant de
modifier les paramètres du modèle, et d'étudier leur impact sur la solution.
"""
import numpy as np

from solvers import rk4
from cosmofunc import age_univ, surd_mini, ti, tf, H, eq_diff, milliard_annee

tf /= milliard_annee
ti /= milliard_annee
age_univ /= milliard_annee


def cond_init(d_i_min, d_i_max, tolerance, delta_eff, t_eff):
    """
    On cherche ici la surdensité minimale telle que le temps d'effondrement soit égal à celui de l'âge de l'univers,
    dans un univers d'Einstein - De Sitter.
    """
    # Chercher le temps d'effondrement pour d_i_min
    init = (d_i_min, d_i_min * H(ti), ti)
    sol1, ttab = rk4(init, eq_diff, tf, max_density=delta_eff)
    temps_eff_max = ttab[-1]
    if temps_eff_max < t_eff:
        print("Surdensité minimum trop grande")
    # Chercher le temps d'effondrement pour d_i_max
    init = (d_i_max, d_i_max * H(ti), ti)
    sol2, ttab = rk4(init, eq_diff, tf, max_density=delta_eff)
    temps_eff_min = ttab[-1]
    if t_eff < temps_eff_min:
        print("Surdensité maximum trop petite")

    inter = abs(d_i_max - d_i_min)
    while inter > tolerance:
        milieu = (d_i_min + d_i_max) / 2
        init = (milieu, milieu * H(ti), ti)
        sol, temps_moy = rk4(init, eq_diff, tf, max_density=delta_eff)
        if temps_moy[-1] > t_eff:
            d_i_min = milieu
        else:
            d_i_max = milieu
        inter = d_i_max - d_i_min

    return (d_i_max + d_i_min) / 2


def iterate_on_param(kind: int, iterateur: np.ndarray, save_res: bool = True):
    choix = {0: 'init', 1: 'eff', 2: 'dt'}
    if kind not in choix.keys():
        raise ValueError("'kind' doit être égal soit à :"
                         "- 0 pour la surdensité initiale"
                         "- 1 pour la surdensité d'effondrement"
                         "- 2 pour le pas de temps.")

    init = (4 * surd_mini, 4 * surd_mini * H(ti), ti)
    results = list()

    if kind == 2:
        for e in iterateur:
            delta, ttab = rk4(init, eq_diff, t_max=tf, dt=e)
            results.append([e, ttab[-1]])
    elif kind == 1:
        for e in iterateur:
            delta, ttab = rk4(init, eq_diff, t_max=tf, max_density=e)
            results.append([e, ttab[-1]])
    else:
        for e in iterateur:
            init = (e, e * H(ti), ti)
            delta, ttab = rk4(init, eq_diff, t_max=tf, max_density=1e4)
            results.append([e, ttab[-1]])

    results = np.array(results)
    if save_res:
        np.savetxt(fname=f"resultats/{choix[kind]}.txt", X=results)

    return results


# ITERATIONS SUR LES SURDENSITÉS INITIALES
# On part de la surdensité qui s'effondre à l'âge de l'univers


if __name__ == '__main__':
    sur_init = cond_init(1e-4, 1e-2, 1e-7, 1e4, age_univ)
    print("Surdensité initiale donnant t_eff égale à l'âge de l'univers :", sur_init)

    surd_init = np.linspace(0.001777656936645508, 15 * 0.001777656936645508, 100)
    iterate_on_param(kind=0, iterateur=surd_init)
    surd_eff = np.linspace(1e3, 1e6, 1000)
    iterate_on_param(kind=1, iterateur=surd_eff)
    N = list(range(5000, 10000))
    iterate_on_param(kind=2, iterateur=N)
