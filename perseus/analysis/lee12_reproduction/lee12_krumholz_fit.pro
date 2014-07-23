function func2, X, P
  Z = 1.0
  a = 0.2
  f_diss = 0.1
  phi_mol = 10.0
  c = 3.0 * 10.0^(10.0)
  mu_H = 2.3 * 10.0^(-24.0)
  sigma_d_Solar = 10.0^(-21.0)
  R_d_Solar = 10.0^(-16.5)
  E_0_Solar = 7.5 * 10.0^(-4.0)
  sigma_d = sigma_d_Solar * Z
  R_d = R_d_Solar * Z
  chi = ((f_diss * sigma_d_Solar * c * E_0_Solar) * (1.0 + (3.1 * Z^(0.365)))) / (31.0 * P[0] * R_d_Solar)

  print, 'chi'
  print, chi

  psi = chi * (2.5 + chi) / (2.5 + (chi * exp(1.0)))

  print, 'psi'
  print, psi
  
  tau_c = (3.0 * X * 2.0 * 10.0^(33.0) * sigma_d) / (4.0 * (3.1 * 10.0^(18.0))^(2.0) * mu_H)
  
  print, 'tau_c'
  print, tau_c
  
  f_H2_sub1 = (3.0 * psi) / (4.0 * tau_c)
  f_H2_sub2 = (4.0 * a * psi * phi_mol) / ((4.0 * tau_c) + (3.0 * (phi_mol - 1.0) * psi))

  print, 'fh2_sub1, fh2_sub2'
  print, f_H2_sub1, f_H2_sub2

  f_H2 = 1.0 - (f_H2_sub1 / (1.0 + f_H2_sub2))
  f_HI = 1.0 - f_H2

  print, 'fhi, fh2'
  print, f_HI, f_H2

  return, f_H2 / f_HI
end
