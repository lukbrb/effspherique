program main
  use solvers, only: rk2, rk4, euler
  use cosmofunc, only: di, ti, t_max
  use paramtune, only: cond_init, iterate_on_param
  implicit none

  real :: dt = 1e-5
  real :: max_density = 1e4
  real, dimension(2) :: resulte, result2, result4

  resulte = euler(4. * di, ti, t_max, dt, max_density)
  result2 = rk2(4. * di, ti, t_max, dt, max_density)
  result4 = rk4(4. * di, ti, t_max, dt, max_density)
  print *, 'Résultat Euler explicit :', resulte
  print *, 'Résultat RK2 :', result2
  print *, 'Résultat RK4 :', result4

  ! call say_hello()

end program main
