from dolfin import *
import numpy as np
import matplotlib.pylab as plt
import sys
sys.path.append("/home/huskeypm/sources/smolhomog/example/myofibril/")
sys.path.append("/net/home/huskeypm/Sources/smolhomog/example/")


import caseRunner as cr # just for PB stuff


cr.tol = 4.
class outerBoundary(SubDomain):
  
  def inside(self,x,on_boundary):
    outer = x[1] > (cr.maxs[1]-cr.tol)  or x[1] < (cr.mins[1]+cr.tol)
    #print "Using outer"
    #if outer:
    #if(1):
    #  print x[1], cr.mins[1],cr.maxs[1], outer
    return outer and on_boundary

class innerBoundary(SubDomain):
  def inside(self,x,on_boundary):
    inner = x[1] <= (cr.maxs[1]+cr.tol) or x[1] >= (cr.mins[1]-cr.tol)
    return inner and on_boundary


# <codecell>


# <codecell>


# ionC [M]
# sigma [C/m^2]
def runCase(meshFile,ionC=0.15,sigma=-0.01):
    mesh = Mesh(meshFile)
    
    cr.mins=np.min(mesh.coordinates(),axis=0)
    cr.maxs=np.max(mesh.coordinates(),axis=0)
    
    
    V = FunctionSpace(mesh,"CG",1)
    dim=2
    subdomains = MeshFunction("uint",mesh,dim-1)
    boundary = outerBoundary();
    outerMarker = 2
    boundary.mark(subdomains,outerMarker) 
  
    #f = Constant(1.)    
    #bc = DirichletBC(V, f, subdomains,outerMarker)
    #view = Function(V)
    #bc.apply(view.vector())
    #File("view.pvd") << view
    
    #import view
    #view.PrintBoundary(mesh,bcs)
    parms = pb.parms
    parms.ionC = ionC # M 
    boundaryPotential = pb.Grahame(sigma,parms.ionC)
    print boundaryPotential
    parms.update()

    bcs=[]
    f = Constant(boundaryPotential)
    bcs.append(DirichletBC(V, f, subdomains,outerMarker))
    
    ## Solve PB eqn 
    (V,potential)= pb.PBEngine(mesh,V,subdomains,bcs)
    File("psi2.pvd") << potential
    scaledPotential = Function(V)
    print "return for now" 
    return 1,2
    
    zLigs = np.array([-1,0,1])
    Ds = np.zeros(np.shape(zLigs)[0])
    for i,zLig in enumerate(zLigs):
      parms.zLig = zLig

      #scaledPotential.vector()[:] = parms.Fz_o_RT*potential.vector()[:]  
      scaledPotential.vector()[:] = parms.F_o_RT*potential.vector()[:]
      scaledPotential.vector()[:] = parms.kT * scaledPotential.vector()[:]  # (z=+,psi=+) --> V
      results = hl.runHomog(fileXML=meshFile,psi=scaledPotential,\
                            reflectiveBoundary="backfront",q=parms.zLig,smolMode=True)
      #print "y dir ", results.d_eff[1] 
      Ds[i] =results.d_eff[0]
        
    return zLigs, Ds

     

# <codecell>
# overwrite definition 
#cr.outerBoundary = outerBoundary
#cr.innerBoundary = outerBoundary
#cr.chargeBoundary=outerBoundary

def test():
  meshFile="meshShallow.xml"
  mesh = Mesh(meshFile)
  cr.mins=np.min(mesh.coordinates(),axis=0)
  cr.maxs=np.max(mesh.coordinates(),axis=0)

  z,Ds=cr.runCase(meshFile=meshFile,ionC=5.00,sigma=-0.1,chargeBoundary=outerBoundary)     
  print Ds

def run():
  meshFile="meshShallow.xml"
  mesh = Mesh(meshFile)
  cr.mins=np.min(mesh.coordinates(),axis=0)
  cr.maxs=np.max(mesh.coordinates(),axis=0)

  class empty:pass
  shallow = empty()
  shallow.meshFile="meshShallow.xml"
  shallow.zLigs,shallow.Ds=cr.runCase(meshFile=shallow.meshFile,\
    ionC=0.15,sigma=-0.01,chargeBoundary=outerBoundary)
  
  deep = empty()
  deep.meshFile="meshDeep.xml"
  deep.zLigs,deep.Ds=cr.runCase(meshFile=deep.meshFile,\
    ionC=0.15,sigma=-0.01,chargeBoundary=outerBoundary)   
  
  N=np.shape(deep.zLigs)[0]
  ind = np.arange(N)  # the x locations for the groups
  width = 0.35       # the width of the bars
  
  fig, ax = plt.subplots()
  rects1 = ax.bar(ind, shallow.Ds, width, color='r')
  rects2 = ax.bar(ind+width, deep.Ds, width, color='y')
  
  ax.set_ylabel('D')
  ax.set_xlabel('z')
  ax.set_title('Corrugated mesh with $\sigma$=-0.01 [C/m2]')
  ax.set_xticks(ind+width)
  ax.set_xticklabels( deep.zLigs )
  
  ax.legend((rects1[0], rects2[0]), ('Shallow', 'Deep') )
  plt.gcf().savefig("corrugated.png")
  

#!/usr/bin/env python
import sys
#
# Revisions
#       10.08.10 inception
#


def helpmsg():
  scriptName= sys.argv[0]
  msg="""
Purpose: 
 
Usage:
"""
  msg+="  %s -validation" % (scriptName)
  msg+="""
  
 
Notes:

"""
  return msg

if __name__ == "__main__":
  import sys
  msg = helpmsg()
  remap = "none"

  if len(sys.argv) < 2:
      raise RuntimeError(msg)

  fileIn= sys.argv[1]
  if(len(sys.argv)==3):
    print "arg"

  for i,arg in enumerate(sys.argv):
    if(arg=="-validation"):
      run()            
      quit()
    if(arg=="-test"):
      test()            
      quit()
  





  raise RuntimeError("Arguments not understood")




