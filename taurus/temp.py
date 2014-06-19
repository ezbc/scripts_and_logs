import numpy as np
import matplotlib.pyplot as plt

def polyfit2d(x, y, z, degree=1):

    ''' 2D polynomial fit to z at the given x and y coordinates. z(x,y) =
    p[0] * x**deg + p[1] * y**deg + ... + p[deg].
    '''

    from scipy.optimize import curve_fit

    popt, pcov = curve_fit(poly2d_fit, (x, y), z)

    return popt, pcov

def poly2d((x, y), coeffs=(0, 1, 1), degree=1):

    if len(coeffs) != 2 * degree + 1:
        raise ValueError('Coeffslength must be 2*degree + 1')
    else:
    	pass

    z = coeffs[0] * np.ones(x.shape)

    for i in xrange(1, degree + 1, 2):
        z += coeffs[i] * x**i + coeffs[i + 1] * y**i

    return z.ravel()

def poly2d_fit((x, y), a=0, b=1, c=1):

    return poly2d((x,y), coeffs=(a, b, c), degree=1)

shape = (10, 10)

x = np.linspace(0, shape[0], shape[0])
y = np.linspace(0, shape[1], shape[1])

x, y = np.meshgrid(x,y)

z = poly2d((x, y), coeffs=(1,15,5))

z_noisy = z + 15 * np.random.normal(size=z.shape)

print z

fig = plt.figure()
ax = fig.add_subplot(111)
ax.imshow(z_noisy.reshape(10, 10))
fig.show()

popt, pcov = polyfit2d(x, y, z_noisy)

z_fit = poly2d_fit((x, y), *popt)

print z_fit

fig = plt.figure()
ax = fig.add_subplot(111)
ax.imshow(z_fit.reshape(10, 10))
fig.show()


