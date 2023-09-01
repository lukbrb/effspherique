module cosmofunc
    implicit none
    private

    public :: H, Hvec, G, eq_diff_lin, eq_diff, di, ti, t_max, age_univers, rho_m, masse, PI, exp_param
    
    real, parameter :: Mp = 3.08567758 * 1E22  ! Valeur d un mégaparsec
    real, parameter :: H0 = 7 * 1e4 / Mp
    real, parameter :: G = 6.6742 * 1E-11
    real, parameter :: Gyr = 3600. * 24. * 365. * 1E09
    real, parameter :: ti = (300000. * 365. * 24. * 3600.) / Gyr
    real, parameter :: di = 0.0017104343414306644
    real, parameter :: t_max = (1 / H0) / Gyr
    real, parameter :: age_univers =  (2 / (3 * H0)) / Gyr
    real, parameter :: masse = 1E16
    real, parameter :: PI = 4.0 * DATAN(1.d0)
    real, parameter :: rho_crit = (3 * H0 ** 2) / (8 * PI * G)  ! Densité critique
    real, parameter :: sig = 1.  ! densité de matière

contains

    real function H(t) result(hubble_cst)
        real, intent(in) :: t
        hubble_cst = 2 / (3 * t)
    end function H

    real function eq_diff_lin(d, p, t) result(val_eq)
        real, intent(in) :: d, p, t
        val_eq = -2 * H(t) * p + 1.5 * (H(t) ** 2) * d
    end function eq_diff_lin

    real function eq_diff(d, p, t) result(val_eq)
        real, intent(in) :: d, p, t
        val_eq = -2 * H(t) * p + 1.5 * (H(t) ** 2) * d * (1. + d) + (4. / 3.) * ((p ** 2) / (1. + d))
    end function eq_diff

! --------------- Vectorized functions (until I learn how to handle this better) -------------------
    function Hvec(t) result(hubble_cst)
        real, dimension(:), intent(in) :: t
        real, dimension(size(t)) :: hubble_cst
        hubble_cst = 2 / (3 * t)
    end function Hvec

    function exp_param(t) result(a)
        real, dimension(:), intent(in) :: t
        real, dimension(size(t)) :: a
        a = (1.5 * H0 * t) ** (2. / 3.)
    end function exp_param

    function rho_m(t) result(densite_moy)
        real, dimension(:), intent(in) :: t
        real, dimension(size(t)) :: densite_moy
        densite_moy = (sig * rho_crit) / (exp_param(t) ** 3)
    end function rho_m


    function eq_diff_linvec(d, p, t) result(val_eq)
        real, dimension(:), intent(in) :: d, p, t
        real, dimension(size(t)) :: val_eq
        val_eq = -2 * Hvec(t) * p + 1.5 * (Hvec(t) ** 2) * d
    end function eq_diff_linvec

    function eq_diffvec(d, p, t) result(val_eq)
        real, dimension(:), intent(in) :: d, p, t
        real, dimension(size(t)) :: val_eq
    val_eq = -2 * Hvec(t) * p + 1.5 * (Hvec(t) ** 2) * d * (1. + d) + (4. / 3.) * ((p ** 2) / (1. + d))
    end function eq_diffvec
    
end module cosmofunc