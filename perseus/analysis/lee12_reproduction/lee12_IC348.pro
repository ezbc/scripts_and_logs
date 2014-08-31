pro IC348_all_params 

; .sav files are copied from /d/leffe2/lee/perseus_cloud/Min/R_H2/R_H2_102311/errors/. 
restore, /ver, 'R_H2_plot_IC348_all.sav'
restore, /ver, 'R_H2_plot_sigma_IC348_all.sav'

; 1000 Gaussian random numbers 
ran_num1 = randomn(seed, 1000)
ran_num2 = randomn(seed, 1000)

tot_sden_arr = dblarr(n_elements(tot_sden),1000)
R_H2_arr = dblarr(n_elements(R_H2),1000)
for i = 0, 999 do begin
  tot_sden_arr[*,i] = (ran_num1[i] * tot_sden_sig) + tot_sden
  R_H2_arr[*,i] = (ran_num2[i] * R_H2_sig) + R_H2
endfor

; Version 1 
; Z varies from 0.5 to 1.5 with an initial guess of 1.0. 
; phi_CNM varies from 1.0 to 20.0 with an initial guess of 8.0. 
; phi_mol = 10.0 is fixed.
; IDL> print, median(Z_arr, /even), stddev(Z_arr)
;      1.0255398      0.11705563
; IDL> print, median(phi_CNM_arr, /even), stddev(phi_CNM_arr)
;      8.0685647      0.98124434
; IDL> print, median(bestnorm_arr, /even), stddev(bestnorm_arr)
;      18.265689       12.427484
; IC348_params_ver1 and IC348_all_params_ver1 are consistent within 1-sigma uncertainties. 
;Z_arr = dblarr(1000)
;phi_CNM_arr = dblarr(1000)
;bestnorm_arr = dblarr(1000)
;err = dblarr(n_elements(tot_sden)) + 1.0
;parinfo = replicate({value:0.D, limited:[0,0], limits:[0.D,0]}, 2)
;parinfo[0].limited[0] = 1
;parinfo[0].limits[0] = 0.5D
;parinfo[0].limited[1] = 1
;parinfo[0].limits[1] = 1.5D
;parinfo[1].limited[0] = 1
;parinfo[1].limits[0] = 1.0D
;parinfo[1].limited[1] = 1
;parinfo[1].limits[1] = 20.0D
;parinfo[*].value = [1.0D, 8.0D]
;for i = 0, 999 do begin
;  tot_sden_fit = reform(tot_sden_arr[*,i])
;  R_H2_fit = reform(R_H2_arr[*,i])
;  result = mpfitfun('func1', tot_sden_fit, R_H2_fit, err, parinfo=parinfo, bestnorm=bestnorm)
;  Z_arr[i] = result[0]
;  phi_CNM_arr[i] = result[1]
;  bestnorm_arr[i] = bestnorm
;endfor

;save, Z_arr, phi_CNM_arr, bestnorm_arr, filename='IC348_all_params_ver1.sav', /verbose

; Version 2
; phi_CNM varies from 1.0 to 20.0 with an initial guess of 8.0.
; Z = 1.0 is fixed.
; phi_mol = 10.0 is fixed.
; IDL> print, median(phi_CNM_arr, /even), stddev(phi_CNM_arr)
;      8.3736688      0.87665885
; IDL> print, median(bestnorm_arr, /even), stddev(bestnorm_arr)
;      18.462623       12.493571
; IC348_params_ver2 and IC348_all_params_ver2 are consistent within 1-sigma uncertainties. 
; Version 1 and Version 2 are consistent within 1-sigma uncertainties. 
phi_CNM_arr = dblarr(1000)
bestnorm_arr = dblarr(1000)
err = dblarr(n_elements(tot_sden)) + 1.0
parinfo = replicate({value:0.D, limited:[0,0], limits:[0.D,0]}, 1)
parinfo[0].limited[0] = 1
parinfo[0].limits[0] = 1.0D
parinfo[0].limited[1] = 1
parinfo[0].limits[1] = 20.0D
parinfo[*].value = [8.0D]

fname='/d/bip3/ezbc/perseus/data/lee12_reproduction/IC348_h_sd.txt' 
OPENW,1,fname 
PRINTF,1, tot_sden,FORMAT='(F7.2,1X)'
CLOSE,1

fname='/d/bip3/ezbc/perseus/data/lee12_reproduction/IC348_RH2.txt' 
OPENW,1,fname 
PRINTF,1, R_H2,FORMAT='(F7.2,1X)'
CLOSE,1

for i = 0, 999 do begin
  tot_sden_fit = reform(tot_sden_arr[*,i])
  R_H2_fit = reform(R_H2_arr[*,i])
  ;result = mpfitfun('func2', tot_sden_fit, R_H2_fit, err, parinfo=parinfo, bestnorm=bestnorm)
  ;result = 0
  ;bestnorm = 0
  phi_CNM_arr[i] = result[0]
  bestnorm_arr[i] = bestnorm
endfor

save, phi_CNM_arr, bestnorm_arr, filename='IC348_all_params_ver2.sav', /verbose

stop
end
