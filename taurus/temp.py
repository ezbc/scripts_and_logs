



# import external modules
import pyfits as pf
import numpy as np
from mycoords import make_velocity_axis
from scipy import ndimage
from skimage.morphology import watershed
from skimage.feature import peak_local_max
import matplotlib.pyplot as plt
from sklearn.feature_extraction import image as sklearn_image
from sklearn.cluster import spectral_clustering


# define directory locations
output_dir = '/d/bip3/ezbc/taurus/data/python_output/nhi_av/'
figure_dir = '/d/bip3/ezbc/taurus/figures/'
av_dir = '/d/bip3/ezbc/taurus/data/av/'
hi_dir = '/d/bip3/ezbc/taurus/data/galfa/'

# Load av fits file
av_image, av_header = pf.getdata(av_dir + 'taurus_av_k09_regrid.fits',
        header=True)

image = av_image[110:250,110:280]

mask = image > 8

label_im, nb_labels = ndimage.label(mask)

# Select larger regions
sizes = ndimage.sum(mask, label_im, range(nb_labels + 1))
mask_size = sizes < 1
remove_pixel = mask_size[label_im]
label_im[remove_pixel] = 0
labels = np.unique(label_im)
label_clean = np.searchsorted(labels, label_im)
label_im, nb_labels = ndimage.label(mask)

fig, axes = plt.subplots(ncols=2, figsize=(12, 5))
ax0, ax1 = axes

ax0.imshow(label_im, cmap=plt.cm.spectral, origin='lower')

ax1.imshow(label_clean, cmap=plt.cm.spectral, origin='lower')

for ax in axes:
    ax.axis('off')

fig.subplots_adjust(hspace=0.01, wspace=0.01, top=1, bottom=0, left=0,
                    right=1)
fig.show()


