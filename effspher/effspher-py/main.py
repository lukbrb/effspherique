import numpy as np

from cosmofunc import age_univ, surd_mini, ti, tf, eq_diff, eq_diff_lin
from paramtune import cond_init, iterate_on_param
from viriel import Surdensite

PARAMS = False
VIRIEL = True
if PARAMS:
    print("Détermination de \delta_i tel que teff = âge de l'univers")
    sur_init = cond_init(1e-4, 1e-2, 1e-7, 1e4, age_univ)
    print("Surdensité initiale donnant t_eff égal à l'âge de l'univers :", sur_init)

    print("Détermination de l'influence des paramètres :")
    print("1 - Itération sur la surdensité initiale...")
    surd_init = np.linspace(0.001777656936645508, 15 * 0.001777656936645508, 100)
    iterate_on_param(kind=0, iterateur=surd_init)

    print("2 - Itération sur la surdensité d'effondrement...")
    surd_eff = np.linspace(1e3, 1e6, 1000)
    iterate_on_param(kind=1, iterateur=surd_eff)

    print("3 - Itération sur le pas de temps...")
    N = np.array(list(range(5000, 10000)))
    iterate_on_param(kind=2, iterateur=N)

    print("Détermination de la surdensité viriel:")

if VIRIEL:
    surdensite_nonlin = Surdensite(4 * surd_mini, ti, eq_diff, tf, dt=1e-5, max_density=4 * 1e7)
    delta_ta_calc = surdensite_nonlin.surd_ta()
    delta_ta_th = (9/16) * np.pi**2 - 1
    surd_vir_calc = surdensite_nonlin.surd_finale()
    surd_vir_th = 18 * np.pi**2

    print("Régime non-linéaire")
    print('-' * 10)
    print("Surdensité volte-face:", delta_ta_calc)
    print("Surdensité volte-face théorique:", delta_ta_th)
    print(f"Différence relative: {(abs(delta_ta_calc - delta_ta_th))/delta_ta_th : %}")
    print('-' * 10)
    print("Surdensité viriel finale:", surd_vir_calc)
    print("Surdensité viriel théorique:", surd_vir_th)
    print(f"Différence relative: {(abs(surd_vir_calc - surd_vir_th))/surd_vir_th : %}")
    print('-' * 10)
    # =================== Linéaire ===================
    surdensite_lin = Surdensite(4 * surd_mini, ti, eq_diff_lin, tf, dt=1e-5, max_density=1000)
    delta_ta_calc_lin = surdensite_lin.surd_ta(energie=False)
    delta_ta_th_lin = (3 / 20) * (6 * np.pi)**(2/3)
    surd_eff_calc_lin = surdensite_lin.surd_eff(equation=eq_diff_lin)
    surd_eff_th_lin = (3 / 5) * (1.5 * np.pi)**(2/3)

    print("Régime linéaire")
    print('-' * 10)
    print("Surdensité volte-face:", delta_ta_calc_lin)
    print("Surdensité volte-face théorique:", delta_ta_th_lin)
    print(f"Différence relative: {(abs(delta_ta_calc_lin - delta_ta_th_lin)) / delta_ta_th_lin : %}")
    print('-' * 10)
    print("Surdensité effondrement:", surd_eff_calc_lin)
    print("Surdensité effondrement théorique:", surd_eff_th_lin)
    print(f"Différence relative: {(abs(surd_eff_calc_lin - surd_eff_th_lin)) / surd_eff_th_lin : %}")


# if PLOT:
#     plt.figure()
#     plt.title("Évolution du rayon")
#     plt.plot(ttab, R, 'r', label="Rayon")
#     plt.plot(tmax, Rta, 'x', label="Rayon max")
#     plt.plot(tvir, Rvir, 'x', label="Rayon viriel")
#     plt.xlabel('Temps (Gyr)')
#     plt.ylabel('R(t)')
#     # plt.yscale('log')
#     plt.legend()
#
#     drdt = vitesse(R)
#     plt.figure()
#     plt.title("Profil des vitesses")
#     plt.plot(ttab, drdt, '--k', label="Vitesse")
#     plt.vlines(tmax, drdt.min(), drdt.max(), colors='red', label="Vitesse au rayon max")
#     plt.vlines(tvir, drdt.min(), drdt.max(), colors='orange', label="Vitesse au rayon viriel")
#     plt.xlabel('Temps (Gyr)')
#     plt.ylabel('R(t)')
#     # plt.yscale('log')
#     plt.legend()
#
#     plt.figure()
#     plt.title("Évolution de la surdensité")
#     plt.plot(ttab, delta, '--k')
#     plt.xlabel('Temps (Gyr)')
#     plt.ylabel('Surdensité $\delta$')
#     plt.yscale('log')
#     plt.show()