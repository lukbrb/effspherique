module viriel
    use cosmofunc, only: Hvec, G, di, ti, rho_m, masse, PI, exp_param
    use io, only: read_data
    use solvers, only: rk4
    implicit none

    private
    character(256), parameter :: filename = 'results/rk4-results.dat'
    public :: surd_vir, surd_finale, surd_ta
contains

    function vitesse(t, r, ddelta, delta) result(v)
        real, dimension(:), allocatable, intent(in) :: t, r, ddelta, delta
        real, dimension(:), allocatable :: v
        allocate(v(size(t)))
        v = Hvec(t) * r - (r * ddelta) / (3 * (1 + delta))
    end function vitesse

    ! function energie_cin(t, r, ddelta, delta) result(ec)
    !     real, dimension(:), allocatable, intent(in) :: t, r, ddelta, delta
    !     real, dimension(size(r)) :: v, ec
    !     v = vitesse(t, r, ddelta, delta)
    !     ec = (3 * masse * v ** 2) / 10
    ! end function energie_cin

    ! function energie_pot(r) result(ep)
    !     real, intent(in) :: r
    !     real :: ep
    !     ep = - (3 * G * masse ** 2) / (5 * r)
    ! end function energie_pot

    function find_rta(t, r, ddelta, delta, energie) result(rta_and_tta)
        logical, intent(in) :: energie
        real, allocatable, dimension(:), intent(in) :: t, r, ddelta, delta
        real, dimension(2) :: rta_and_tta
        integer, dimension(1) :: imax
        if (energie) then
            imax = minloc(abs(vitesse(t, r, ddelta, delta)))
        else
            imax = maxloc(r)
        end if
        rta_and_tta = [r(imax(1)), t(imax(1))]
    end function find_rta

    function surd_ta() result(sta)
        real, dimension(:), allocatable:: t, ddelta, delta, r
        real, dimension(2) :: rta_et_tta, sol
        real :: sta
        call read_data(t, ddelta, delta)
        allocate(r(size(t)))
        r = ((3 * masse) / (4 * PI * rho_m(t) * (1 + delta))) ** (1. / 3.)
        rta_et_tta = find_rta(t, r, ddelta, delta, energie=.true.)
        deallocate(ddelta)
        deallocate(delta)
        sol = rk4(4. * di, ti, t_max=rta_et_tta(2), dt=1E-5, max_density=1E4)
        sta = sol(2)
    end function surd_ta


    subroutine find_rvir(t_and_r, t, r, ddelta, delta)
        real, dimension(:), allocatable, intent(in) :: t, r, ddelta, delta
        real, dimension(:), allocatable :: vvir
        real, dimension(2) :: t_and_r
        real :: rvir, tvir
        integer, dimension(1) :: index_vir

        allocate(vvir(size(t)))
        vvir = -sqrt(G * masse / r)
        index_vir = minloc(abs(vitesse(t, r, ddelta, delta) - vvir))
        deallocate(vvir)

        tvir = t(index_vir(1))
        rvir = r(index_vir(1))
        t_and_r = [tvir, rvir]        
    end subroutine find_rvir

    function surd_vir() result(svir)
        real, dimension(:), allocatable:: t, r, ddelta, delta
        real, dimension(2) :: t_and_r, sol
        real :: svir
        call read_data(t, ddelta, delta)
        ! write(*, '(5F10.4)') delta
        allocate(r(size(t)))
        r = ((3 * masse) / (4 * PI * rho_m(t) * (1 + delta))) ** (1. / 3.)
        call find_rvir(t_and_r, t, r, ddelta, delta)
        deallocate(ddelta)
        deallocate(delta)
        sol = rk4(4. * di, ti, t_max=t_and_r(1), dt=1E-5, max_density=1E4)
        svir = sol(2)
    end function surd_vir

    function surd_finale() result(dfinal)
        real, dimension(:), allocatable:: t, ddelta, delta
        real, dimension(2) :: tvir_dvir
        real, dimension(1) :: tfin, tvir, dfinal
        call read_data(t, ddelta, delta)
        deallocate(ddelta)
        deallocate(delta)
        tfin(1) = t(size(t))
        tvir_dvir = surd_vir()
        tvir(1) = tvir_dvir(1)
        dfinal = 1 + tvir_dvir(2) * (exp_param(tfin) / exp_param(tvir)) ** 3
    end function surd_finale

end module viriel