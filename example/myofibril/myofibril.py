
import sys 
sys.path.append("/net/home/huskeypm/Sources/modified-pb/example/")

import poissonboltzmann as pb
from dolfin import *
import matplotlib.pylab as plt 
import numpy as np
import homoglight as hl

print "WARNING: bounds are hardcoded"
dims = 2


# Define Dirichlet boundary (x = 0 or x = 1)
class BcThang():
  def __init__(self,dims=2,tol=10):
    self.dims = dims
    self.tol = tol
    self.mins=np.zeros(dims)
    self.maxs=np.array([547.8,316.272])

  def isOuter(self,x,y,z=0):
    outerXY =  ((x>self.maxs[0]-self.tol or x<self.mins[0]+self.tol) or\
                (y>self.maxs[1]-self.tol or y<self.mins[1]+self.tol)) 
    
    outerZ = True

    if dims==3:
      outerZ =  (z>self.maxs[2]-self.tol or  z<self.mins[2]+self.tol)

    if(outerZ and outerXY):
      #print "TRUE ", x,y
      #quit()
      return True

    #print "False ", x,y
    return False            
            

bcThang = BcThang()
                                                    

class outerBoundary(SubDomain):
  def inside(self,x,on_boundary):
    outer  = bcThang.isOuter(x[0],x[1]) # ,mins,maxs)

    #if(outer==False and on_boundary==True): 
    #  print "sdfsdf"
    #  quit()
    return outer and on_boundary

# Define Dirichlet boundary (x = 0 or x = 1)
class innerBoundary(SubDomain):
  def inside(self,x,on_boundary):
    outer  = bcThang.isOuter(x[0],x[1]) # ,mins,maxs)
    if(outer):
        return False
    return on_boundary



# sigma [C/m^2]
# ionC  [M] 
#  
# pass in boundary class where you want charge applied
def runCase(ionC=0.15,meshFile="0.xml",sigma=-0.01,singleCase=False, chargeBoundary=innerBoundary):
    # load mesh 
    mesh = Mesh(meshFile)
    
    
    #dims
    bcThang.mins=np.min(mesh.coordinates(),axis=0)
    bcThang.maxs=np.max(mesh.coordinates(),axis=0)

    #bcThang.mins=np.array([240.8,240.272])
    #bcThang.maxs=np.array([250.8,250.272])

    # env/system info
    parms = pb.parms
    parms.ionC = ionC # M 
    boundaryPotential = pb.Grahame(sigma,parms.ionC)
    print boundaryPotential
    parms.update()
    
    
    # Mark BCs 
    V = FunctionSpace(mesh,"CG",1)
    subdomains = MeshFunction("uint",mesh,dims-1)
    boundary = chargeBoundary(); 
    chargeMarker = 2
    boundary.mark(subdomains,chargeMarker)
    #boundary = innerBoundary()
    #innerMarker = 3
    #boundary.mark(subdomains,innerMarker)
    #NOT USING ANY NEUMANN COND ds = Measure("ds")[subdomains]
    
    bcs=[]    
    f = Constant(boundaryPotential)    
    bcs.append(DirichletBC(V, f, subdomains,chargeMarker))    

    ## Solve PB eqn 
    (V,potential)= pb.PBEngine(mesh,V,subdomains,bcs)
    File("psi.pvd") << potential
    scaledPotential = Function(V)
 
    if(singleCase):
      zLigs = np.array([-1]) 
    else: 
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
      runCase(ionC=0.15,meshFile="0.xml",sigma=-0.01,singleCase=True)   
    if(arg=="-run"):
      doit()


