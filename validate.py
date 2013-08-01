#
# Creates 2D mesh 
#



import numpy as np
from dolfin import * 

def FDsym(x,marg=2., midpoint=5., slope=10.):
  midpointLeft = midpoint-marg
  midpointRight = midpoint+marg
  zL = 1/(np.exp(-(slope*x-slope*midpointLeft)) + 1)
  zR = 1 - 1/(np.exp(-(slope*x-slope*midpointRight)) + 1)
  zP = zL*zR
  return zP 

def FD2d(x,y,amplitude=1., marg=2., midpoint=5., slope=10.):
  zx = FDsym(x,marg=marg,midpoint=midpoint,slope=slope)
  zy = FDsym(y,marg=marg,midpoint=midpoint,slope=slope)

  z = amplitude*zx*zy
  return z



  
def doit():
  import homog 
  
  ## Grid setup 
  mesh = UnitSquare(50,50)
  mesh.coordinates()[:] = mesh.coordinates()[:] * 10.
  V = FunctionSpace(mesh,"CG",1) 
  
  ## subdomains (do nothing)
  subdomains = MeshFunction("uint",mesh,1)
  
  ## PMF setup 
  print "WARNING: need to use dof map here"
  pmf  = Function(V)  
  amplitude = 2. # [kcal/mol]
  marg=3.;
  midpoint=5.;
  slope=10.
  xs=mesh.coordinates()[:,0];ys=mesh.coordinates()[:,1]
  
  vfunc = np.vectorize(FD2d,otypes=[np.float])
  pmfVals = vfunc(xs,ys, -1*amplitude,marg,midpoint,slope)
  pmf.vector()[:] = pmfVals
  
  ## write files, but should be able to do this without writing files 
  molPrefix = "temp"
  meshFileInner = molPrefix+"_mesh.xml.gz"
  subdomainFileInner = molPrefix+"_subdomains.xml.gz"
  potentialFileInner = molPrefix+"_values.xml.gz"
  File(meshFileInner) << mesh
  File(subdomainFileInner) << subdomains
  File(potentialFileInner) << pmf
  
  results = homog.SolveHomogSystem(
    molPrefix=molPrefix,smolMode=True,molGamer=False)


#!/usr/bin/env python
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
 
Usage:
"""
  msg+="  %s -validation" % (scriptName)
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
    if(arg=="-validation"):
      doit()






