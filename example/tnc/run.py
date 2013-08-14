
import numpy as np
import matplotlib.pylab as plt
from dolfin import *

mesh = Mesh("p_mesh.xml.gz")

# <codecell>

np.max(mesh.coordinates())



import poissonboltzmann as pb
from dolfin import *
import matplotlib.pylab as plt 
import numpy as np
import homoglight as hl

print "WARNING: bounds are hardcoded"
mins=np.array([-50,-50,-50])
maxs=np.array([50,50,50])

# <codecell>

# Define Dirichlet boundary (x = 0 or x = 1)
tol = 0.6
def isInner(x,y,z,mins=np.array([0,0,0]),maxs=np.array([1,1,1])):
  if(x>mins[0]+tol and x<maxs[0]-tol and\
     y>mins[1]+tol and y<maxs[1]-tol and\
     z>mins[2]+tol and z<maxs[2]-tol
     ):
        return True

  return False            
            
                                                    

class leftBoundary(SubDomain):
  def inside(self,x,on_boundary):
    outer  = isInner(x[0],x[1],x[2],mins,maxs)
    return outer and on_boundary

# Define Dirichlet boundary (x = 0 or x = 1)
class elseBoundary(SubDomain):
  def inside(self,x,on_boundary):
    outer  = isInner(x[0],x[1],x[2],mins,maxs)
    if(outer):
        return False
    return on_boundary


def runCase(ionC=0.15):
    # load mesh 
    meshFile = "p_mesh.xml.gz"
    mesh = Mesh(meshFile)
    
    
    #dims
    mins=np.min(mesh.coordinates(),axis=0)
    maxs=np.max(mesh.coordinates(),axis=0)
    print maxs

    # env/system info
    parms = pb.parms
    sigma = -0.01 # C/m^2
    parms.ionC = ionC # M 
    boundaryPotential = pb.Grahame(sigma,parms.ionC)
    print boundaryPotential
    parms.update()
    
    
    # Mark BCs 
    V = FunctionSpace(mesh,"CG",1)
    dim=3
    subdomains = MeshFunction("uint",mesh,dim-1)
    boundary = leftBoundary(); 
    boundary.min=mins;boundary.max=maxs
    leftmarker = 2
    boundary.mark(subdomains,leftmarker)
    boundary = elseBoundary()
    boundary.min=mins;boundary.max=maxs
    elsemarker = 3
    boundary.mark(subdomains,elsemarker)
    #NOT USING ANY NEUMANN COND ds = Measure("ds")[subdomains]
    
    bcs=[]    
    f = Constant(boundaryPotential)    
    bcs.append(DirichletBC(V, f, subdomains,leftmarker))    
    #import view
    #view.PrintBoundary(mesh,bcs)

    ## Solve PB eqn 
    (V,potential)= pb.PBEngine(mesh,V,subdomains,bcs)
    File("psi.pvd") << potential
    scaledPotential = Function(V)
 
    zLigs = np.array([-2,-1,0,1,2])
    Ds = np.zeros(np.shape(zLigs)[0])
    for i,zLig in enumerate(zLigs):
      parms.zLig = zLig

      #scaledPotential.vector()[:] = parms.Fz_o_RT*potential.vector()[:]  
      scaledPotential.vector()[:] = parms.F_o_RT*potential.vector()[:]
      scaledPotential.vector()[:] = parms.kT * scaledPotential.vector()[:]  # (z=+,psi=+) --> V
      results = hl.runHomog(fileXML=meshFile,psi=scaledPotential,q=parms.zLig,smolMode=True)
      Ds[i] =results.d_eff[0]
        
    return zLigs,Ds    

# <codecell>

def doit():
    zLigs,Ds1p000 = runCase(ionC=1.)
    zLigs,Ds0p150 = runCase(ionC=0.15)

    zLigsLabs=['ATP(-2)','ADP(-1)','(0)','K(+1)',"Ca(2+)"]
    width=0.3
    fig,ax=plt.subplots()
    ax.bar(zLigs,Ds0p150,width,color='r',label="I=0.15 [M]")
    ax.bar(zLigs+width,Ds1p000,width,color='b',label="I=1.00 [M]")
    ax.set_ylim([0,1])
    ax.set_xlabel("Net charge") 
    ax.set_ylabel("D")
    ax.legend(loc=0)
    ax.set_xticks(zLigs+width)
    ax.set_xticklabels(zLigsLabs)
    plt.gcf().savefig("myofibril.png")

# <codecell>


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

# <codecell>


