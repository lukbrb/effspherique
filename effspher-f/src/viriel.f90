module viriel
    use cosmofunc, only: Hvec, G, di, ti, rho_m, masse, PI, exp_param

    use solvers, only: rk4
    implicit none

    private
    character(256), parameter :: filename = 'results/rk4-results.dat'
    public :: surd_vir, surd_finale
contains

    subroutine read_t_and_delta(t, ddelta, delta)
        real, dimension(:), allocatable, intent(out) :: t, ddelta, delta
        character(len=512) :: msg
        integer :: num_rows, i
        integer :: io, status
        logical :: exists
        
        inquire(file=filename, exist=exists)
        if (exists) then
            open(unit=io, file=filename, status='old', action='read', iostat=status, iomsg=msg)
            if (status /= 0) then
                write(*,*) "Erreur lors de l'ouverture du fichier ", trim(filename), ':', trim(msg)
                stop
            end if
        end if
        
        ! Compter le nombre de lignes dans le fichier
        num_rows = 0
        do
            read(io, *, iostat=status)
            if (status /= 0) exit
            num_rows = num_rows + 1
        end do
    
        ! Allouer de la mémoire pour les tableaux
        allocate(t(num_rows))
        allocate(ddelta(num_rows))
        allocate(delta(num_rows))
    
        ! Revenir au début du fichier
        rewind(io)
        ! Lire les données dans les tableaux
        do i = 1, num_rows
            read(io, *) t(i), ddelta(i), delta(i)
        end do

        close(io)
        
        end subroutine read_t_and_delta

    function vitesse(t, r, ddelta, delta) result(v)
        real, dimension(:), allocatable, intent(in) :: t, r, ddelta, delta
        real, dimension(:), allocatable:: v
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

    subroutine find_rvir(t_and_r, t, r, ddelta, delta)
        real, dimension(:), allocatable, intent(in) :: t, r, ddelta, delta
        real, dimension(:), allocatable :: vvir, vdiff
        real, dimension(2) :: t_and_r
        real :: rvir, tvir
        integer, dimension(1) :: index_vir

        allocate(vvir(size(t)))
        allocate(vdiff(size(t)))
        vvir = -sqrt(G * masse / r)
        vdiff = vitesse(t, r, ddelta, delta) - vvir
        index_vir = minloc(abs(vdiff))
        deallocate(vvir)
        deallocate(vdiff)
        tvir = t(index_vir(1))
        rvir = r(index_vir(1))
        t_and_r = [tvir, rvir]        
    end subroutine find_rvir

    function surd_vir() result(sol)
        real, dimension(:), allocatable:: t, r, ddelta, delta
        real :: tvir
        real, dimension(2) :: t_and_r, sol
        call read_t_and_delta(t, ddelta, delta)
        ! write(*, '(5F10.4)') delta
        allocate(r(size(t)))
        r = ((3 * masse) / (4 * PI * rho_m(t) * (1 + delta))) ** (1. / 3.)
        call find_rvir(t_and_r, t, r, ddelta, delta)
        tvir = t_and_r(1)
        print *, 'tvir=', tvir
        sol = rk4(4. * di, ti, t_max=tvir, dt=1E-5, max_density=1E4)
    end function surd_vir

    function surd_finale() result(dfinal)
        real, dimension(:), allocatable:: t, ddelta, delta
        real, dimension(2) :: tvir_dvir
        real, dimension(1) :: tfin, tvir, dfinal
        call read_t_and_delta(t, ddelta, delta)
        deallocate(ddelta)
        deallocate(delta)
        tfin(1) = t(size(t))
        tvir_dvir = surd_vir()
        tvir(1) = tvir_dvir(1)
        dfinal = 1 + tvir_dvir(2) * (exp_param(tfin) / exp_param(tvir)) ** 3
    end function surd_finale

end module viriel