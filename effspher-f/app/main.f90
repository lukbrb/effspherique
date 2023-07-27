program main
  use solvers, only: rk2, rk4
  use cosmofunc, only: di, ti, t_max
  implicit none

  real :: dt = 1e-5
  real :: max_density = 1e4
  real, dimension(2) :: result2, result4

  result2 = rk2(di, ti, t_max, dt, max_density)
  result4 = rk4(di, ti, t_max, dt, max_density)
  print *, 'Résultat RK2:', result2
  print *, 'Résultat RK4:', result4
  ! call say_hello()

end program main
