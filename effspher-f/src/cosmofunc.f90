module cosmofunc
    implicit none
    private

    public :: H, G, eq_diff_lin, eq_diff, di, ti, t_max, age_univers, rho_m, masse, PI, exp_param
    
    real(kind = 8), parameter :: Mp = 3.08567758 * 1E22  ! Valeur d un mégaparsec
    real(kind = 8), parameter :: H0 = 7 * 1e4 / Mp
    real(kind = 8), parameter :: G = 6.6742 * 1E-11
    real(kind = 8), parameter :: Gyr = 3600. * 24. * 365. * 1E09
    real(kind = 8), parameter :: ti = (300000. * 365. * 24. * 3600.) / Gyr
    real(kind = 8), parameter :: di = 0.0017104343414306644
    real(kind = 8), parameter :: t_max = (1 / H0) / Gyr
    real(kind = 8), parameter :: age_univers =  (2 / (3 * H0)) / Gyr
    real(kind = 8), parameter :: masse = 1E16
    real(kind = 8), parameter :: PI = 4.0 * DATAN(1.d0)
    real(kind = 8), parameter :: rho_crit = (3 * H0 ** 2) / (8 * PI * G)  ! Densité critique
    real(kind = 8), parameter :: sig = 1.  ! densité de matière

    interface H
        procedure H_scalar, H_vec
    end interface H

    interface eq_diff_lin
        procedure eq_diff_lin_scalar, eq_diff_lin_vec
    end interface eq_diff_lin

    interface eq_diff
        procedure eq_diff_scalar, eq_diff_vec
    end interface eq_diff


contains

    real(kind = 8) function H_scalar(t) result(hubble_cst)
        real(kind = 8), intent(in) :: t
        hubble_cst = 2 / (3 * t)
    end function H_scalar

    real(kind = 8) function eq_diff_lin_scalar(d, p, t) result(val_eq)
        real(kind = 8), intent(in) :: d, p, t
        val_eq = -2 * H(t) * p + 1.5 * (H(t) ** 2) * d
    end function eq_diff_lin_scalar

    real(kind = 8) function eq_diff_scalar(d, p, t) result(val_eq)
        real(kind = 8), intent(in) :: d, p, t
        val_eq = -2 * H(t) * p + 1.5 * (H(t) ** 2) * d * (1. + d) + (4. / 3.) * ((p ** 2) / (1. + d))
    end function eq_diff_scalar

! --------------- Vectorized functions (until I learn how to handle this better) -------------------
    function H_vec(t) result(hubble_cst)
        real(kind = 8), dimension(:), intent(in) :: t
        real(kind = 8), dimension(size(t)) :: hubble_cst
        hubble_cst = 2 / (3 * t)
    end function H_vec

    function exp_param(t) result(a)
        real(kind = 8), dimension(:), intent(in) :: t
        real(kind = 8), dimension(size(t)) :: a
        a = (1.5 * H0 * t) ** (2. / 3.)
    end function exp_param

    function rho_m(t) result(densite_moy)
        real(kind = 8), dimension(:), intent(in) :: t
        real(kind = 8), dimension(size(t)) :: densite_moy
        densite_moy = (sig * rho_crit) / (exp_param(t) ** 3)
    end function rho_m


    function eq_diff_lin_vec(d, p, t) result(val_eq)
        real(kind = 8), dimension(:), intent(in) :: d, p, t
        real(kind = 8), dimension(size(t)) :: val_eq
        val_eq = -2 * H(t) * p + 1.5 * (H(t) ** 2) * d
    end function eq_diff_lin_vec

    function eq_diff_vec(d, p, t) result(val_eq)
        real(kind = 8), dimension(:), intent(in) :: d, p, t
        real(kind = 8), dimension(size(t)) :: val_eq
    val_eq = -2 * H(t) * p + 1.5 * (H(t) ** 2) * d * (1. + d) + (4. / 3.) * ((p ** 2) / (1. + d))
    end function eq_diff_vec
    
end module cosmofunc