pro getPlanckDataTaurus

; Script for extracting planck data and building a header for the image

; First download the planck map
; Centered on RA = 4 33m, Dec = 26d 47m
; Centered on RA = 68.3 deg, Dec = 26.8 deg
; Choose 10 deg area, thus 600 pix * 1/60 deg/pix = 10 deg
; Starting 5 deg east and north of center

img = getimgfromhealedit([68,26],[1000,1000],[1./60,1./60],/dust,$
                         column=3,hdr=hdr)
writefits,'taurus.planck.ebv.fits',img,hdr

img = getimgfromhealedit([68,26],[1000,1000],[1./60,1./60],/dust,$
                         column=4,hdr=hdr)
writefits,'taurus.planck.ebverr.fits',img,hdr

end
