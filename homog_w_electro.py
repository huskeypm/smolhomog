#
# Revisions
#
"""
export PYTHONPATH=$PYTHONPATH:/home/huskeypm/bin/dolfin/
"""
from dolfin import *
import numpy as np
import sys
sys.path.append("/home/huskeypm/sources/smolfin/")
import poissonboltzmann as pb

boundsMin = np.zeros(3)
boundsMax = np.zeros(3)
meshFile="none"
dim=-1
electro=1.

class params:
  beta = 1/0.59
  q=1.
  z = -3

EPS = 0.001 
def DefineBoundary(x,btype,on_boundary):
  if(not on_boundary):
    return 0 

  # need this, since for some reason 'dim' is not global
  dim = np.shape(x)[0]


  lr = ( np.abs(x[0]-boundsMin[0]) < EPS or np.abs(x[0]-boundsMax[0]) < EPS)
  tb = ( np.abs(x[1]-boundsMin[1]) < EPS or np.abs(x[1]-boundsMax[1]) < EPS)
  if(dim==3):
    fb = ( np.abs(x[2]-boundsMin[2]) < EPS or np.abs(x[2]-boundsMax[2]) < EPS)
  else:
    fb=0
  obs = (not tb and not lr and not fb)

  if(btype=="lr"):
    return lr

  if(btype=="obs"):
    return obs

  if(btype=="tb"):
    return tb  

  if(btype=="fb"):
    return fb  

class XBoundary(SubDomain):
  def inside(self, x, on_boundary):
    #return (( x[0] < 0 or x[0] > 7) and on_boundary)
    return (DefineBoundary(x,"lr",on_boundary))

class YBoundary(SubDomain):
  def inside(self, x, on_boundary):
    return (DefineBoundary(x,"tb", on_boundary))

class ZBoundary(SubDomain):
  def inside(self, x, on_boundary):
    return (DefineBoundary(x,"fb", on_boundary))


class ObsBoundary(SubDomain):
  def inside(self, x, on_boundary):
    return (DefineBoundary(x,"obs", on_boundary))

##
##
##

def doit(meshFile="none",tag="",meshType="gamer"): 


  ## Mesh
  if 0:
    1
  else:
    mesh = Mesh(meshFile)     
    #gamer=0
    #
    for i in np.arange( mesh.ufl_cell().geometric_dimension() ):
      print "Marking direction %d" % i
      boundsMin[i] = np.min(mesh.coordinates()[:,i])
      boundsMax[i] = np.max(mesh.coordinates()[:,i])

  # Function spaces
  V = VectorFunctionSpace(mesh,"CG",1)
  Vscalar = FunctionSpace(mesh,"CG",1)

  ## first do electro
  if electro:
   print "Doing pb"
   pb.params.z = params.z
   pb.params.q = params.q
   (Vdumm,psi) = pb.SolvePoissonBoltzmann(mesh,meshType=meshType)
   name = "elect.pvd" 
   File(name) << psi      

  else:
   print "Neglecting PB"
   psi = Function(Vscalar)
   psi.vector()[:] = 0.

  ## NOw homogenization  
  
  dim = mesh.ufl_cell().geometric_dimension()
  print "Using %d dim mesh" % dim
  print boundsMin
  print boundsMax
  
  # Test boundaries
  test= 0
  if (test==1):
    V = FunctionSpace(mesh,"CG",1)
    x = Function(V)
    #bc_lr = DirichletBC(V,Constant(1),XBoundary())
    #bc_lr.apply(x.vector())
    #bc_obs = DirichletBC(V,Constant(20),ObsBoundary())
    #bc_obs.apply(x.vector())
    #bc_tb = DirichletBC(V,Constant(5),YBoundary())
    #bc_tb.apply(x.vector())
    bc_fb = DirichletBC(V,Constant(7),ZBoundary())
    bc_fb.apply(x.vector())
    File("testbound.pvd") << x
    quit()
    
  
  bcs = [] 
  if 0:
    1
  else:
    bcs.append(DirichletBC(V.sub(0),Constant(0),XBoundary()))
    bcs.append(DirichletBC(V.sub(1),Constant(0),YBoundary()))
    if(dim==3):
      bcs.append(DirichletBC(V.sub(2),Constant(0),ZBoundary()))
    
  if(dim==2):
    Dii  = Constant((1.0,1.0))
  elif(dim==3):
    Dii  = Constant((1.0,1.0,1.0))
  Aij = diag(Dii) 
  #Delta = Identity( mesh.ufl_cell().geometric_dimension()) #
  Delta = Identity( dim )
  
  u = TrialFunction(V)
  v = TestFunction(V)
  if(meshType=="gamer"):
    form = inner(Aij*(grad(u) + Delta), grad(v))*dx(1)
  else:
    form = inner(Aij*(grad(u) + Delta), grad(v))*dx
  a = lhs(form)
  L = rhs(form)
  
  x = Function(V)

  # NOT GOOD FOR MOST LARGE MESHES solve(a == L, x, bcs)
  # Use linear variational solver 
  # Equation is linear, since diffusion coefficient doesn't depend on u,
  print "verify that this is avtually a linear poroblem"
  lvproblem = LinearVariationalProblem(a,L, x, bcs=bcs)        
  solver = LinearVariationalSolver(lvproblem)
  solver.parameters["linear_solver"] = "gmres"
  solver.parameters["preconditioner"] = "ilu"
  solver.solve()

  #plot(u,interactive=True)
  print "WARNING: not sure that this is correct yet for adding electro contribution"
  expnpmf = Function(Vscalar)
  expnpmf.vector()[:] = np.exp(-1*params.beta*params.q*psi.vector()[:])
  File("distro.pvd") << project(expnpmf)

  #xproj = project(x[0]*expnpmf,V=Vscalar)
  xproj = project(x*expnpmf)
  name = "%dd%s.pvd" % (dim,tag)
  File(name) << xproj
  
  omegas = np.zeros(dim)
  for i in range(dim):
    v = [0,0,0]
    v[i] = 1
    if(dim==2):
      grad_Xi_component = inner(grad(x[i]),Constant((v[0],v[1]))) + Constant(1)
    elif(dim==3):
      grad_Xi_component = inner(grad(x[i]),Constant((v[0],v[1],v[2]))) + Constant(1)
    grad_Xi_component = grad_Xi_component*expnpmf
    print "WARNING: not sure that this is correct yet for adding electro contribution"
    outname = "diff%d.pvd" % i
    File(outname)<<project(grad_Xi_component,V=Vscalar)
  
    if(meshType=="gamer"): 
      form = grad_Xi_component * dx(1)
    else:
      form = grad_Xi_component * dx
    integrand = assemble(form)
    omegas[i] = integrand
  
  if(meshType=="gamer"): 
    vol = assemble( Constant(1)*dx(1),mesh=mesh ) 
    surf = assemble( Constant(1)*ds,mesh=mesh ) 
  else:
    vol = assemble( Constant(1)*dx,mesh=mesh ) 
    surf = assemble( Constant(1)*ds,mesh=mesh ) 
  
  print "D*"
  Ds = omegas/vol
  print Ds
  print "Dx"
  print Ds[0]

  print "Vol: %f SA: %f\n" % (vol, surf)
  return Ds[0]
  
  # for debug
  # Only relevant for cubes
  #print "Vol bounds  volfrac  2/3p 2p/2-p"
  #v = np.prod( (boundsMax-boundsMin)[0:dim])
  #vf = vol /v
  #print "%f %f %f %f" % (v,vf,2/3.*vf, 2*vf/(3-vf))

if __name__ == "__main__":
  msg="""
Purpose: 

Usage:
  homog3d.py 

  assumes sphere2d.xml
Notes:

"""
  tag = "_welectro"

  error="""
Results don't agree w homog.py, so don't use
 
"""

  raise RuntimeError(error) 

  import sys
  if len(sys.argv) < 2:
      raise RuntimeError(msg)

  fileIn= sys.argv[1]
  if(len(sys.argv)==3):
    print "arg"

  for i,arg in enumerate(sys.argv):
    if(arg=="-mesh"):
      meshFile=sys.argv[i+1]
    if(arg=="-negative"):
      q = -1 * q
      tag = "_negative"
    if(arg=="-noelectro"):
      electro=False
      tag = "_noelectro"


  meshFile = "sphere2d.xml"







  doit(meshFile=meshFile,tag=tag)
