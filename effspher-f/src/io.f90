module io
    implicit none

    private
    character(256), parameter :: filename = 'results/rk4-results.dat'
    public :: read_data

contains

    subroutine read_data(t, ddelta, delta)
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

        end subroutine read_data

end module io