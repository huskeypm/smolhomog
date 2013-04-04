import sys
import numpy as np
import matplotlib.pylab as plt
sys.path.append("/home/huskeypm/sources/smolfin/")
sys.path.append("/home/huskeypm/sources/homogenization/")

from dolfin import *
#import imp
#pb  = imp.load_source('poissonboltzmann.py', '/home/huskeypm/sources/smolfin/')
import poissonboltzmann as pb
meshName = "meshes/volFrac_0.10.xml"

# quick test
#mesh = Mesh(meshName)
#import poissonboltzmann as pb
#pb.params.molRad = 12.5
#pb.params.domRad = 200
#(V,x) = pb.SolvePoissonBoltzmann(mesh,meshType="gamer")

import homog
debug = 0
z = -0.09  # warning, values above 0.2 seem to give unphysical values
molRad = 12.5
homog.parms.d = 1. # diff constant 

#import homog_w_electro as hwe
  #Dx = hwe.doit(meshFile="meshes/volFrac_0.10.xml",meshType="gamer")

# molRad - radius of enzeyme 
# z - enzyme charge
# q - substrate charge 
def CalcPotentialtential(mesh,molRad,z,q,meshType="dolfin"):
   pb.params.z = z
   pb.params.q = q
   pb.params.molRad = 12.5
   pb.params.sphericalDomainBoundary=False
   (Vdumm,psi) = pb.SolvePoissonBoltzmann(mesh,meshType=meshType)
   name = "elect.pvd"
   File(name) << psi

   return psi 

  

def call(meshPrefix,q,volFrac,molGamer=0):
#  hwe.params.z = -1
#  hwe.pb.molRad = 12.5
#  hwe.pb.sphericalDomainBoundary=False
#  hwe.params.molRad = 12.5
#  hwe.params.sphericalDomainBoundary=False
#  hwe.params.q=1.


  # mesh 
  meshName = meshPrefix+"_mesh.xml.gz"
  mesh = Mesh(meshName)

  # ESP 
  #psi = CalcPotentialtential(mesh,molRad,z,q,meshType="gamer")
  psi = CalcPotentialtential(mesh,molRad,z,q)


  
  # double check size just in case
  mn = np.min(mesh.coordinates(),axis=0)
  mx = np.max(mesh.coordinates(),axis=0)
  boxVol = np.prod(mx-mn)
  sphereVol = 4/3.*np.pi*(molRad**3)
  print meshName
  print "Vol frac %f vs %f " % (sphereVol/boxVol,volFrac)


  #molDomUnit.smolMode = smolMode
  #solve_homogeneous_unit(molDomUnit,type="field",debug=debug,smolMode=smolMode)
  #homog.SolveHomogSystem(molPrefix=meshPrefix,molGamer=molGamer)  
  #quit()
  results = homog.SolveHomogSystem(molPrefix=meshPrefix,molGamer=molGamer,\
    smolMode="true",smolPsi=psi,smolq = q)
  Dx = results.molDomUnit.problem.d_eff[0]

  #(V,x) = pb.SolvePoissonBoltzmann(mesh)
  # do electro homog
  # store value
  #results[i,j] = np.random.rand() +  i*j

  # moleculardomain.q =q 
  # moleculardomain.psi = psi
  
  #Dx = hwe.doit(meshFile=meshName, meshType="gamer")
  return Dx
  
  
def doit(fileIn):
  ## params 

  volFracs = np.array([0.1,0.2,0.27,0.34,0.5])
  #meshes = np.array([0.1,0.2,0.27,0.34,0.5])
  meshes = volFracs
  # kappa hard coded into PB solver
  qs = np.array([-1,0,1,2])
  #qs = np.array([-1,0])           
  
  
  if(debug==1):
    qs=[0]
    volFracs=[0.1]
    qs=[2]
    volFracs=[0.5]
    meshes = volFracs

  results = np.zeros([ np.shape(volFracs)[0], np.shape(qs)[0] ])
  
  for i, volFrac in enumerate(volFracs):
    for j, q in enumerate(qs):
      #meshName = "meshes/volFrac_%4.2f.xml" % meshes[i]
      #meshName = "meshes/volFrac_%4.2f_mesh.xml.gz" % meshes[i]
      meshPrefix= "meshes/volFrac_%4.2f" % meshes[i]
      #molGamer = 1
      molGamer = 0
      results[i,j] = call(meshPrefix,q,volFrac,molGamer=molGamer)

  if(debug==1):
    return
  
  
  plt.figure()
  col = ["r-","k-","b-","b--"]
  for j, q in enumerate(qs):
    label = "qSubstrate = %d " % q
    plt.plot(volFracs,results[:,j], col[j],label=label)
  
  plt.title("Diffusion hindrance with respect to volume fraction")
  plt.xlabel("Volume fraction")
  plt.ylabel("Dx/Dbulk")
  plt.legend(loc='lower left')
  
  plt.gcf().savefig("volfrac.png")
  








import sys
#
# Revisions
#       10.08.10 inception
#


if __name__ == "__main__":
  import sys
  scriptName= sys.argv[0]
  msg="""
Purpose: 
  Compute diffusion constant for a variety of geometries 
 
Usage:
"""
  msg+="  %s -debug" % (scriptName)
  msg+="""
  
 
Notes:

"""
  remap = "none"

  if len(sys.argv) < 2:
      raise RuntimeError(msg)

  fileIn= sys.argv[1]
  if(len(sys.argv)==3):
    print "arg"

  for i,arg in enumerate(sys.argv):
    if(arg=="-debug"):
      debug=1




  doit(fileIn)


