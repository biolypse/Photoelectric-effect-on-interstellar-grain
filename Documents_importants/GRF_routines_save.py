#!/usr/bin/python

from numpy import *
from matplotlib.pylab import *
import numpy as np



def addGRF(cube, sigma_dens, Adens, beta, Npad, Kmin, beta_in, smooth=None):

   ### Get shape and size of the cube
   dims = shape(cube)
   nAxes = len(dims)   # number of axes


   ### Generate a GRF cube

   # Create a padded cube in Fourier space
   Cubedims = zeros(nAxes, int32)
   for i in range(nAxes):
      Cubedims[i] = dims[i] + Npad

   FourierCube = ones(Cubedims, complex)

   # Create cubes of Fourier coordinates
   centre = Cubedims/2
   Coords = indices(Cubedims, complex)
   for i in range(nAxes): Coords[i] -= centre[i]
   Kradius = sqrt( sum(Coords**2,axis=0) )

   # Fill the Fourier cube with the requiered power law(s)
   FourierCube = Kradius**(-beta*0.5)
   m_in = nonzero(Kradius<Kmin)
   FourierCube[m_in] = (Kmin**(-beta*0.5) / Kmin**4) * Kradius[m_in]**4

   # Fix NaNs
   FourierCube[nonzero(FourierCube!=FourierCube)] = 0.

   # Randomize the phases of the Fourier cube
   Phase = ones(Cubedims, complex)
   Phase = np.random.seed()   # necessary to have different random number in parallel threads (if forked at the same time, they get the same seed unless np.random.seed() is called)
   Phase = uniform(size=Cubedims)*2*pi*1j
   Phase *= exp(Phase)
   FourierCube *= Phase

   # Multiply by a smoothing Gaussian
   if (smooth):
      Ksmooth = float(dims[0])/smooth
      FourierCube *= exp(-(Kradius)**2/(2*Ksmooth**2))

   # Compute the density cube by inverse Fourier transform the Fourier cube
   RealCube = fftn(ifftshift(FourierCube))

   Density = ones(dims, float32)
   Indices = []
   for i in range(nAxes):
      Indices.append(slice(Npad/2,Cubedims[i]-Npad/2-(Npad%2)))
   Indices = tuple(Indices)

   Density = exp( sigma_dens*real(RealCube[Indices]) )
   Density /= mean(Density.flatten())

   return Density * Adens


def addGRF2(cube, sigma_dens, Adens, beta, Npad, Kmin, beta_in):

   ### Get shape and size of the cube
   dims = shape(cube)
   nAxes = len(dims)   # number of axes


   ### Generate a GRF cube

   # Create a padded cube in Fourier space
   Cubedims = zeros(nAxes, int32)
   for i in range(nAxes):
      Cubedims[i] = dims[i] + Npad

   FourierCube1 = ones(Cubedims, complex)
   FourierCube2 = ones(Cubedims, complex)

   # Create cubes of Fourier coordinates
   centre = Cubedims/2
   Coords = indices(Cubedims, complex)
   for i in range(nAxes): Coords[i] -= centre[i]
   Kradius = sqrt( sum(Coords**2,axis=0) )

   # Fill the Fourier cube with the requiered power law(s)
   FourierCube1 = exp( -beta*0.5 * log(Kradius) )
   FourierCube2 = exp( -beta*0.5 * log(Kradius) )
   m = nonzero(Kradius>10)
   FourierCube2[m] = exp( -beta*0.6 * log(Kradius[m]) )
   loglog(Kradius.flatten(), FourierCube1.flatten(), '.b')
   loglog(Kradius.flatten(), FourierCube2.flatten(), '.r')
   show()
   m_in = nonzero(Kradius<Kmin)
#   FourierCube[m_in] = (Kmin**(-beta*0.5) / Kmin**4) * Kradius[m_in]**4
   FourierCube1[m_in] = 0.
   FourierCube2[m_in] = 0.

   # Fix NaNs
   FourierCube1[nonzero(FourierCube1!=FourierCube1)] = 0.
   FourierCube2[nonzero(FourierCube2!=FourierCube2)] = 0.

   # Randomize the phases of the Fourier cube
   Phase = ones(Cubedims, complex)
   Phase = random(Cubedims)*2*pi*1j
   Phase *= exp(Phase)
   FourierCube1 *= Phase
   FourierCube2 *= Phase

   # Compute the density cube by inverse Fourier transform the Fourier cube
   RealCube1 = fftn(ifftshift(FourierCube1))
   RealCube2 = fftn(ifftshift(FourierCube2))
   Density1 = ones(dims, float32)
   Density2 = ones(dims, float32)
   Indices = []
   for i in range(nAxes):
      Indices.append(slice(Npad/2,Cubedims[i]-Npad/2))
   Indices = tuple(Indices)

   Density1 = Adens * exp( sigma_dens*real(RealCube1[Indices]) )
   Density1 /= mean(Density1.flatten())
   Density2 = Adens * exp( sigma_dens*real(RealCube2[Indices]) )
   Density2 /= mean(Density2.flatten())

   return Density1, Density2




def addGRF_to_CRTfile(modelname, sigma_dens, Adens, beta, Npad, Kmin, beta_in):

   ### Load the Galaxy cube
   filename = '%s.crt' % (modelname)

   fp = file(filename, 'rb')

   nX, nY, nZ = fromfile(fp, int32, 3)
   print nX, nY, nZ

   cube = fromfile(fp, float32, nX*nY*nZ)
   Cube = cube.reshape((nZ, nY, nX))


   ### Generate a GRF cube

   # Create a padded cube in Fourier space
   Ncube1 = nZ + Npad   # 300 for N=256 (700 <=> ~26 Gbytes of RAM memory)
   Ncube2 = nY + Npad   # 300 for N=256 (700 <=> ~26 Gbytes of RAM memory)
   Ncube3 = nX + Npad   # 300 for N=256 (700 <=> ~26 Gbytes of RAM memory)
   Cubedims = asarray([Ncube1,Ncube2,Ncube3], int32)
   FourierCube = ones((Ncube1,Ncube2,Ncube3), complex)

   # Create cubes of Fourier coordinates
   centre=[Ncube1/2,Ncube2/2,Ncube3/2]
   X, Y, Z   = indices(Cubedims, complex)
   X        -= centre[0]
   Y        -= centre[1]
   Z        -= centre[2]
   Kradius = sqrt(X*X+Y*Y+Z*Z)

   # Fill the Fourier cube with the requiered power law(s)
   FourierCube = Kradius**(-beta*0.5)
   m_in = nonzero(Kradius<Kmin)
   FourierCube[m_in] = (Kmin**(-beta*0.5) / Kmin**4) * Kradius[m_in]**4

   # Fix NaNs
   FourierCube[nonzero(FourierCube!=FourierCube)] = 0.

   # Randomize the phases of the Fourier cube
   Phase = ones((Ncube1,Ncube2,Ncube3), complex)
   Phase = random(Cubedims)*2*pi*1j
   Phase *= exp(Phase)
   FourierCube *= Phase

   # Compute the density cube by inverse Fourier transform the Fourier cube
   RealCube = fftn(ifftshift(FourierCube))
   Density = ones((nZ,nY,nX), float32)

   npad = 0    # can be usefull to take only the inner part of the RealCube
   Density[npad:nZ-npad,npad:nY-npad,npad:nX-npad] = \
      exp( \
         sigma_dens*real(RealCube[Npad/2:Ncube1-Npad/2, \
                                    Npad/2:Ncube2-Npad/2, \
                                    Npad/2:Ncube3-Npad/2]) \
         + log(Adens) \
         )[npad:nZ-npad,npad:nY-npad,npad:nX-npad]

   Density /= mean(Density.flatten())

   # Use the GRF density cube to scale the input density cube
   Galaxy = Cube * Density

   # Plot coldens check
   if (0):
      coldens = ones((nY,nX))
      coldens = sum(Density, axis=0)
      coldens = coldens.flatten()
      hist(log10(coldens))
      show()

   # Drop Galaxy in a binary file for CRT
   if (1):
      GalaxyCube = real((Galaxy)).flatten()

      fp = file('%s.crt' % modelname, 'wb')
      dims = asarray([nX,nY,nZ], int32)
      dims.tofile(fp)
      GalaxyCube.tofile(fp)
      fp.close()




