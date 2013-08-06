#
# Still in development 
#
from dolfin import *
import matplotlib.pyplot as plt
import sys
import numpy as np
sys.path.append("/home/huskeypm/sources/modified-pb/example")
import poissonboltzmann as pb
parms = pb.parms

from scipy.interpolate import griddata

# for interpolating solutions onto 2d mesh 
def interp2d(mesh,x):
      maxr = parms.domRad    
      minr = -parms.domRad    
      res = 10000
      (gx,gy) = np.mgrid[0:0:1j,minr:maxr:res*1j]
      interp = griddata(mesh.coordinates(),x.vector(),(gx,gy))
      interp[np.isnan(interp)]=0
      interp = np.reshape(interp,res)
      gy = np.reshape(gy,res)

      return gy,interp


# for validating PB equation against Fig 3 of Bourbatacha
def test1():
  plt.figure()
  parms.mode='linear'
  # Mesh is made according to Bourbatache
  # R = 10^-8 [M] --> 100 nm --> 1000A
  # separation is 
  # 2e-8 m --> 20 nm, which implies R + 10 nm per box half-width --> 1100A
  fileIn= "/net/home/huskeypm/Sources/smolhomog/example/validation/cylinder.xml"
  fileIn= "/home/huskeypm/sources/smolhomog/example/validation/cylinder.xml"
  parms.molRad=1000  # A 
  parms.domRad=1100. # A 
  mesh = Mesh(fileIn)
  
  # DH   
  #parms.ionC = 0.1 # M 
  #parms.z = -2
  #boundaryPotential="DH"
  #(V,x)= SolvePoissonBoltzmann(mesh,boundaryPotential=boundaryPotential)
  #gy,interp=interp2d(mesh,x)
  #plt.plot(gy,interp,"b",label="DH [I] %3.1f M" % parms.ionC)
  
  # VALIDATED Potential looks validated against Fig 3 of Bourbatche 2012 paper
  sigma = -0.01 # C/m^2
  cb = 50 * 1e-3 # mol/m^3 --> M
  parms.ionC = cb
  boundaryPotential = pb.Grahame(sigma,cb)
  (V,x)= pb.SolvePoissonBoltzmann(mesh,boundaryPotential=boundaryPotential)
  gy,interp=interp2d(mesh,x)
  plt.plot(gy,interp,"r",label="psi0=%4.1f mV [I] %3.1f M" % (boundaryPotential,parms.ionC))

  cb = 200 * 1e-3 # mol/m^3 --> M
  parms.ionC = cb
  boundaryPotential = pb.Grahame(sigma,cb)
  (V,x)= pb.SolvePoissonBoltzmann(mesh,boundaryPotential=boundaryPotential)
  gy,interp=interp2d(mesh,x)
  plt.plot(gy,interp,"g",label="psi0=%4.1f mV [I] %3.1f M" % (boundaryPotential,parms.ionC))
   
  plt.legend(loc=0)
  plt.xlim([parms.molRad,parms.domRad])
  #plt.ylim([0,-45])
  plt.xlabel("1e-8 m")
  plt.gca().invert_yaxis()
  plt.gcf().savefig("bourbfig3.png")


  #Fig. 9. Relative homogenized diffusion coefficient 
  #report Dhom/D for ep=0 -- 1.
  # for sig = 0.07, 0.05,  0.001
  #might also compare w Fig 7 of Murad, but should probably suffice to show opposite
  #tends noted for D- vs D+ [actually, they aren't opposites, so look closer]  


# for validating against fig 9 of bourbatache 
def test2():
  sys.path.append("/home/huskeypm/sources/smolhomog/example/noobstacle/")
  import testing as test
  
  # all cases 
  sigmas = np.array([0.07, 0.05,  0.001 ]) 
  nSigma = np.shape(sigmas)[0]
  molRads = np.linspace(100,900,10)
  nMolRads = np.shape(molRads)[0]  
  outs = np.zeros([nSigma,nMolRads])
  cb = 50 * 1e-3 # mol/m^3 --> M  
  vFracs = np.zeros(nMolRads)  
  for i,molRad in enumerate(molRads): 
        print molRad
        for j,sigma in enumerate(sigmas):
            boundaryPotential = pb.Grahame(sigma,cb)
            #(V,potential)= pb.SolvePoissonBoltzmann(mesh,boundaryPotential=boundaryPotential)
            #
            #diff test.doit(mode="hack2",fileIn=fileIn,potential=potential)
            #outs[j,i] = Deff
            # vFracs[i] = V
            
            
  ## SINGLE CASE   
  # add routine to do a number of mesh sizes   
  # now need 1e-8 m cell (100 A)   
  fileIn= "/home/huskeypm/sources/smolhomog/example/validation/cylinder2.xml"
  
        
        
  parms.molRad= 900  # A 
  parms.domRad=1000. # A 
  #mesh = Mesh(fileIn)

  mesh = Mesh(fileIn)
  # VALIDATED Potential looks validated against Fig 3 of Bourbatche 2012 paper

  i=0; 
  j=0
  sigma = -0.01 # C/m^2
  sigma = sigmas[j] 
  cb = 50 * 1e-3 # mol/m^3 --> M
  parms.ionC = cb
  boundaryPotential = pb.Grahame(sigma,cb)
  #(V,potential)= pb.SolvePoissonBoltzmann(mesh,boundaryPotential=boundaryPotential)
  potential = Function(FunctionSpace(mesh,"CG",1))
  potential.vector()[:] = 0.  
  results = test.doit(mode="hack2",discontinuous=False,dim=2,fileIn=fileIn,potential=potential)
  vFracs[i] = results.phi
  outs[j,i] = results.Ds[0]


  plt.figure()
  plt.plot(vFracs, outs[j,:],label="$\sigma$=%5.1f" % sigmas[j] )
  plt.legend(loc=0)
  plt.gcf().savefig("fig9.png")           
  
    
    





      

#test1()
test2()
