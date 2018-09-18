import numpy as np
import pylab as pl

from GRF_routines import addGRF


"""
   N.B.: the grains are generated randomly => a single grain is not necessarily 
      representative of the fractal parameters. Only a statistically significant
      sample can provide a reliable view of the underlying properties. 
"""

# Tunable parameters
# - Main parameters
N = 100              # size of the grain image
sigma_dens = 1.0     # Width of the density distribution
beta = 3.0           # Slope of the power spectrum (=probability density function=PDF) for k>Kmin

# - Secondary parameters
Kmin = 3             # Minimum wavenumber for PDF slope = beta
beta_in = 4          # Slope of PDF between k=0 and Kmin
Npad = 30            # Padding around the map to get rid of the peak near k=0
Adens = 1.           # Mean density - should not impact dust grain

doplot = 1           # Visual check of grain calculation: 0= desactivated, 1= only final plot, 2= animation of search for the main group + final plot


# while loop to make sure the peak is sufficiently close to the center of the cube
tol = 0.4            # Tolerance for location of the density maximum - between 0 (very tolerant) and 0.49999 (very tough)
mmax = [0,0]
ii = 0
while ((mmax[0]<N*tol) | (mmax[0]>N*(1-tol)) | (mmax[1]<N*tol) | (mmax[1]>N*(1-tol))):

   # Build a fractal cube using Gaussian Random Field
   cube = np.zeros((N,N))
   cube = addGRF(cube, sigma_dens, Adens, beta, Npad, Kmin, beta_in)

   # Search the location of the max value
   mmax = np.nonzero(cube==np.amax(cube))
   print mmax
   ii += 1
   if (ii>1000):
      print "Failed generating a GRF with this tolerance value:", tol
      print "Try a lower value. It should be between 0 (very tolerant) and 0.49999 (very tough)."
      exit()


# Roll the cube to put the max at the center
cube = np.roll(cube, N/2-mmax[0][0], axis=0)
cube = np.roll(cube, N/2-mmax[1][0], axis=1)

# Apodize the cube
I, J = np.indices(np.shape(cube))
Rad = (I-N/2)**2+(J-N/2)**2
sig = N*0.1
cube *= np.exp(-Rad/(2*sig**2))

# Compute the location of the mass center
cdm = [np.sum(J*cube)/np.sum(cube),np.sum(I*cube)/np.sum(cube)]
print cdm

# Roll the cube to put the mass center at the center
cube = np.roll(cube, int(N/2-cdm[0]), axis=0)
cube = np.roll(cube, int(N/2-cdm[1]), axis=1)

# Threshold
grain = np.zeros(np.shape(cube))
m = np.nonzero(cube>np.percentile(cube,70))
grain[m] = 1.

# Identify the pixels attached to the main region
def group(grain, x0,y0):
   # Warning: no precaution for edges
   goon = 1
   progress = np.zeros(np.shape(grain))
   progress[x0,y0] = 2
   if (doplot==2): pl.ion()
   while (goon==1):
      goon = 0
      m = np.nonzero(progress==2)   # search last included friends
      for i in range(len(m[0])):
         i1, j1 = m[0][i], m[1][i]
         # Search only touching pixels
         if (0):
            searchind = [[i1+1,j1], [i1+1,j1-1], [i1,j1-1], [i1-1,j1-1], \
                         [i1-1,j1], [i1-1,j1+1], [i1,j1+1], [i1+1,j1+1]]
         # Search touching pixels + 1
         if (1):
            searchind = [[i1+1,j1], [i1+1,j1-1], [i1,j1-1], [i1-1,j1-1], \
                         [i1-1,j1], [i1-1,j1+1], [i1,j1+1], [i1+1,j1+1], \
                         [i1+2,j1], [i1+2,j1-1], [i1+2,j1-2], [i1+1,j1-2], \
                         [i1,j1-2], [i1-1,j1-2], [i1-2,j1-2], [i1-2,j1-1], \
                         [i1-2,j1], [i1-2,j1+1], [i1-2,j1+2], [i1-1,j1+2], \
                         [i1,j1+2], [i1+1,j1+2], [i1+2,j1+2], [i1+2,j1+1]]
         for inds in searchind:
            inds = np.asarray(inds)
            if ((np.all(inds>=0)) & (np.all(inds<N))):
               if (grain[inds[0],inds[1]]==1):
                  if (progress[inds[0],inds[1]]==0): # not found before
                     progress[inds[0],inds[1]] = 3   # 3 for freshly identified
                     goon = 1
      if (doplot==2): pl.imshow(progress, origin='low', cmap='bone', interpolation='nearest')
      if (doplot==2): pl.pause(0.0001)
      progress[m] = 1
      progress[np.nonzero(progress==3)] = 2

   progress[np.nonzero(progress>1)] = 1
   pl.ioff()
   return progress

mmax = np.nonzero(cube==np.amax(cube))
grain2 = group(grain, mmax[0][0], mmax[1][0])

# Write down the grain image in ASCII format (very inefficient format, but easy to control)
np.savetxt("Grain_N%i_S%ip%i_B%ip%i.txt" % (N, int(sigma_dens), \
                                               int((sigma_dens-int(sigma_dens))*10), \
                                               int(beta), \
                                               int((beta-int(beta))*10)), \
                                               grain2, fmt='%i')


if (doplot>=1): 
   pl.figure(figsize=(20,10))
   pl.subplot(131)
   pl.title('Original gaussian density distribution')
   pl.plot(cdm[0], cdm[1], '+k')
   pl.imshow(np.log10(cube), origin='low', interpolation='nearest')
   pl.colorbar(shrink=0.4)

   pl.subplot(132)
   pl.title('Threshold cut of the original distribution')
   pl.imshow(grain, origin='low', cmap='bone', interpolation='nearest')
   pl.colorbar(shrink=0.4)

   pl.subplot(133)
   pl.title('Final Grain - only the central group')
   pl.imshow(grain2, origin='low', cmap='bone', interpolation='nearest')
   pl.colorbar(shrink=0.4)

   pl.show()

