program main
  use effspher_f, only: say_hello
  use solvers, only: rk4, H
  implicit none

  real :: di, pi, ti, t_max, dt, max_density, Gyr, Mp, H0, result

  Mp = 3.08567758 * 1E22  ! Valeur d un mégaparsec
  H0 = 7 * 1e4 / Mp
  Gyr = 3600 * 24 * 365 * 1E09
  ti = 300000 * 365 * 24 * 3600
  ti = ti / Gyr
  di = 0.0017104343414306644
  pi = H(ti) * di
  t_max = (1 / H0) / Gyr
  dt = 1e-5
  max_density = 1e4


  result = rk4(di, pi, ti, t_max, dt, max_density)
  print *, 'Résultat RK4:', result
  ! call say_hello()

end program main
