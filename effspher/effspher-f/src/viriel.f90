module viriel
    use cosmofunc, only: Hvec, G, di, ti, rho_m, masse, PI, exp_param
    use io, only: read_data
    use solvers, only: rk4
    implicit none

    private
    character(256), parameter :: filename = 'results/rk4-results.dat'
    public :: surd_vir, surd_finale, surd_ta
contains

    function radius(t, delta) result(r)
        real, dimension(:), allocatable, intent(in) :: t, delta
        real, dimension(size(t)) :: r

        r = ((3 * masse) / (4 * PI * rho_m(t) * (1 + delta))) ** (1. / 3.)
    end function radius


    function vitesse(t, ddelta, delta) result(v)
        real, dimension(:), allocatable, intent(in) :: t, ddelta, delta
        real, dimension(size(t)) :: r, v

        r = radius(t, delta)
        v = Hvec(t) * r - (r * ddelta) / (3. * (1. + delta))
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

    function find_rta(t, ddelta, delta, use_energie) result(rta_and_tta)
        logical, intent(in) :: use_energie
        real, allocatable, dimension(:), intent(in) :: t, ddelta, delta
        real, dimension(size(t)) :: r
        real, dimension(2) :: rta_and_tta
        integer, dimension(1) :: imax

        r = radius(t, delta)
        if (use_energie) then
            imax = minloc(abs(vitesse(t, ddelta, delta)))
        else
            imax = maxloc(r)
        end if
        rta_and_tta = [r(imax(1)), t(imax(1))]
    end function find_rta


    function surd_ta() result(sol)
        real, dimension(:), allocatable:: t, ddelta, delta, r
        real, dimension(2) :: rta_and_tta, sol

        call read_data(t, ddelta, delta)
        allocate(r(size(t)))

        r = radius(t, delta)
        rta_and_tta = find_rta(t, ddelta, delta, use_energie=.true.)
        sol = rk4(4. * di, ti, t_max=rta_and_tta(2), dt=1E-5, max_density=1E4)

        deallocate(t)
        deallocate(ddelta)
        deallocate(delta)
        deallocate(r)
    end function surd_ta


    function find_rvir(t, ddelta, delta) result(t_and_r)
        real, dimension(:), allocatable, intent(in) :: t, ddelta, delta
        real, dimension(size(t)) :: r, vvir
        real, dimension(2) :: t_and_r
        real :: rvir, tvir
        integer, dimension(1) :: index_vir

        r = radius(t, delta)
        vvir = -sqrt(G * masse / r)
        index_vir = minloc(abs(vitesse(t, ddelta, delta) - vvir))

        tvir = t(index_vir(1))
        rvir = r(index_vir(1))
        t_and_r = [tvir, rvir]        
    end function find_rvir


    function surd_vir() result(sol)
        real, dimension(:), allocatable:: t, r, ddelta, delta
        real, dimension(2) :: t_and_r, sol

        call read_data(t, ddelta, delta)
        allocate(r(size(t)))

        r = radius(t, delta)
        t_and_r = find_rvir(t, ddelta, delta)
        sol = rk4(4. * di, ti, t_max=t_and_r(1), dt=1E-5, max_density=1E4)
        
        deallocate(t)
        deallocate(ddelta)
        deallocate(delta)
        deallocate(r)
    end function surd_vir


    function surd_finale() result(dfinal)
        real, dimension(:), allocatable:: t, ddelta, delta
        real, dimension(2) :: tvir_dvir
        real, dimension(1) :: tfin, tvir, dfinal

        call read_data(t, ddelta, delta)

        tfin(1) = t(size(t))
        tvir_dvir = surd_vir()
        tvir(1) = tvir_dvir(1)
        dfinal = 1 + tvir_dvir(2) * (exp_param(tfin) / exp_param(tvir)) ** 3

        deallocate(t)
        deallocate(ddelta)
        deallocate(delta)
    end function surd_finale


end module viriel