#
# Creates 2D mesh 
#



import homog 
import numpy as np
from dolfin import * 

def FDsym(x,marg=2., midpoint=5., slope=10.):
  midpointLeft = midpoint-marg
  midpointRight = midpoint+marg
  yL = 1/(np.exp(-(slope*x-slope*midpointLeft)) + 1)
  yR = 1 - 1/(np.exp(-(slope*x-slope*midpointRight)) + 1)
  yP = yL*yR
  return yP 



## Grid setup 
mesh = UnitSquare(50,50)
mesh.coordinates()[:] = mesh.coordinates()[:] * 10.
xs = mesh.coordinates()[:,0]
V = FunctionSpace(mesh,"CG",1) 

## subdomains (do nothing)
subdomains = MeshFunction("uint",mesh,1)

## PMF setup 
print "WARNING: need to use dof map here"
pmf  = Function(V)  
amplitude = 1. # [kcal/mol]
fd1d = np.sqrt(amplitude) * FDsym(xs,marg=1.,midpoint=5.,slope=10.)
fd2d = np.outer(fd1d,fd1d)
pmf.vector()[:] = fd2d

## write files, but should be able to do this without writing files 
meshFileInner = "temp_mesh.xml.gz"
subdomainFileInner = "temp_subdomains.xml.gz"
potentialFileInner = "temp_values.xml.gz"
File(meshFileInner) << mesh
File(subdomainFileInner) << subdomains
File(potentialFileInner) << pmf






