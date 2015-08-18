
import os
os.chdir('/d/bip3/ezbc/shield/')

image = '749237_cube_regrid_casa'
importfits(fitsimage='749237_cube_regrid.fits',
           imagename=image + '.image',
           overwrite=True)

# define spectral axis as velocity
import os
os.system('rm -rf ' + image + '_vel.image')
ia.open(image + '.image')
ia.regrid(outfile=image + '_vel.image',
          asvelocity=True,
          overwrite=True)
ia.close()

# write with velocity axis
exportfits(fitsimage=image + '_vel.fits',
           imagename=image + '_vel.image',
           velocity=True,
           dropstokes=True,
           overwrite=True)

# Write without velocity axis
exportfits(fitsimage=image + '_freq.fits',
           imagename=image + '_vel.image',
           velocity=False,
           dropstokes=True,
           overwrite=True)

