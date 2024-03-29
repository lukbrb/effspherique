"""
Sont regroupées dans ce fichier les différentes fonctions permettant de
modifier les paramètres du modèle, et d'étudier leur impact sur la solution.
"""
import numpy as np

from solvers import integrateur
from cosmofunc import age_univ, surd_mini, ti, tf, H, eq_diff


def cond_init(d_i_min, d_i_max, tolerance, delta_eff, t_eff):
    """
    On cherche ici la surdensité minimale telle que le temps d'effondrement soit égal à celui de l'âge de l'univers,
    dans un univers d'Einstein - De Sitter.
    """
    # Chercher le temps d'effondrement pour d_i_min
    sol1, _, temps_eff_max = integrateur(d_i_min, ti, eq_diff, tf, max_density=delta_eff)

    if temps_eff_max < t_eff:
        print("Surdensité minimum trop grande")
    # Chercher le temps d'effondrement pour d_i_max
    sol2, _, temps_eff_min = integrateur(d_i_max, ti, eq_diff, tf, max_density=delta_eff)

    if t_eff < temps_eff_min:
        print("Surdensité maximum trop petite")

    inter = abs(d_i_max - d_i_min)
    while inter > tolerance:
        milieu = (d_i_min + d_i_max) / 2
        sol, _, temps_moy = integrateur(milieu, ti, eq_diff, tf, max_density=delta_eff)
        if temps_moy > t_eff:
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
            delta, _, teff = integrateur(4 * surd_mini, ti, eq_diff, t_max=tf, dt=e)
            results.append([e, teff])
    elif kind == 1:
        for e in iterateur:
            delta, _, teff = integrateur(4 * surd_mini, ti, eq_diff, t_max=tf, max_density=e)
            results.append([e, teff])
    else:
        for e in iterateur:
            delta, _, teff = integrateur(e, ti, eq_diff, t_max=tf, max_density=1e4)
            results.append([e, teff])

    results = np.array(results)
    if save_res:
        np.savetxt(fname=f"resultats/{choix[kind]}.txt", X=results)

    return results


if __name__ == '__main__':
    sur_init = cond_init(1e-4, 1e-2, 1e-7, 1e4, age_univ)
    print("Surdensité initiale donnant t_eff égale à l'âge de l'univers :", sur_init)

    surd_init = np.linspace(0.001777656936645508, 15 * 0.001777656936645508, 100)
    iterate_on_param(kind=0, iterateur=surd_init)
    surd_eff = np.linspace(1e3, 1e6, 1000)
    iterate_on_param(kind=1, iterateur=surd_eff)
    N = np.array(list(range(5000, 10000)))
    iterate_on_param(kind=2, iterateur=N)
