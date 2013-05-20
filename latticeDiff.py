#
# Time dependent solver, single species 
# 
from dolfin import *
import numpy as np


class LeftBoundary(SubDomain):
  def inside(self,x,on_boundary):
    edge = (np.abs(x[0]- -8) < DOLFIN_EPS) 
    #print x[0], edge, on_boundary
    return on_boundary and edge

class RightBoundary(SubDomain):
  def inside(self,x,on_boundary):
    edge = (np.abs(x[0]- 8) < DOLFIN_EPS) 
    #print x[0], edge, on_boundary
    return on_boundary and edge


# Nonlinear equation 
class MyEquation(NonlinearProblem):
    def __init__(self, a, L,bcs):
        NonlinearProblem.__init__(self)
        self.L = L
        self.a = a
        self.bcs = bcs
        self.reset_sparsity = True
    def F(self, b, x):
        assemble(self.L, tensor=b)
        for bc in self.bcs:
          bc.apply(b,x)
    def J(self, A, x):
        assemble(self.a, tensor=A, reset_sparsity=self.reset_sparsity)
        for bc in self.bcs:
          bc.apply(A)
        self.reset_sparsity = False


def tsolve(Diff=1.,fileName="m25.xml",outName="output.pvd",mode="pointsource"):
  D = Constant(Diff) 
  # Create mesh and define function space
  mesh = Mesh(fileName)     
  V = FunctionSpace(mesh, "Lagrange", 1)
  
  
  # Define trial and test functions
  du    = TrialFunction(V)
  q  = TestFunction(V)
  # Define functions
  u   = Function(V)  # current solution
  u0  = Function(V)  # solution from previous converged step
  #init_cond = InitialConditions()

  ## mark boundaries 
  subdomains = MeshFunction("uint",mesh,1)
  boundary = LeftBoundary()
  lMarker = 2
  boundary.mark(subdomains,lMarker)
  boundary = RightBoundary()
  rMarker = 3
  boundary.mark(subdomains,rMarker)

  ## decide on source of probability density  
  bcs=[]
  if(mode=="pointsource"): 
    f = Expression("sqrt(x[0]*x[0] + x[1]*x[1])")
    limit = 0.3
    vmin = 1
    vmax = 10
    init_cond = conditional(ge(f, limit), vmin, vmax) 
    z = project(init_cond,V)
    File("test.pvd") << z
    #u.interpolate(init_cond)
    #u0.interpolate(init_cond)
    u.vector()[:] = z.vector()[:]
    u.vector()[np.where(u.vector()) < vmin] = vmin
    u.vector()[np.where(u.vector()) > vmax] = vmax
    u0.vector()[:] = u.vector()[:]
  
    print "total conc", assemble(u*dx)

  else: 
    bc = DirichletBC(V,Constant(1.0),subdomains,lMarker)
    bcs.append(bc)
    f = Function(V)
    bc.apply(f.vector())
    File("bc.pvd") << f

  ds = Measure("ds")[subdomains]

  
  # weak form 
  dt=15.0   
  T=200 
  RHS = -inner(D*grad(u), grad(q))*dx
  L = u*q*dx - u0*q*dx - dt * RHS
  
  # Compute directional derivative about u in the direction of du (Jacobian)
  a = derivative(L, u, du)
  
  
  problem = MyEquation(a,L,bcs)
  solver = NewtonSolver("gmres")         
  file = File(outName, "compressed")
  
  
  t = 0.0
  concs=[]
  ts = []
  while (t < T):
      # report on prev iter
      uds = assemble(u*ds(rMarker),mesh=mesh)
      area = assemble(Constant(1.)*ds(rMarker),mesh=mesh)
      conc = uds/area
      ts.append(t)
  
      # advance 
      t += dt
      u0.vector()[:] = u.vector()
      solver.solve(problem, u.vector())
      file << (u,t) 
  
      concs.append(conc)


  ts = np.asarray(ts)
  concs = np.asarray(concs)

  return (ts,concs)
  

def doit():
  import matplotlib.pylab as plt
  plt.figure()
  #
  mode = "pointsource"
  mode = "bc"
  #(ts,concs) = tsolve(Diff=1.,fileName="m25.xml",outName="o25_1.pvd",mode=mode) 
  #plt.plot(ts,concs,"k-",label="Diff=1., m25")
  #
  (ts,concs) = tsolve(Diff=1.0,fileName="m15.xml",outName="o15.pvd",mode=mode) 
  plt.plot(ts,concs,"b-",label="Diff=1.0, m15",lw=1)
  #
  (ts,concs) = tsolve(Diff=1.0,fileName="m25.xml",outName="o25.pvd",mode=mode) 
  plt.plot(ts,concs,"b-",label="Diff=1.0, m25",lw=2)
  #
  (ts,concs) = tsolve(Diff=1.0,fileName="m50.xml",outName="o50.pvd",mode=mode) 
  plt.plot(ts,concs,"b-.",label="Diff=1.0, m50",lw=3)
  #
  (ts,concs) = tsolve(Diff=1.0,fileName="m75.xml",outName="o75.pvd",mode=mode) 
  plt.plot(ts,concs,"b-.",label="Diff=1.0, m75",lw=4)
  #
  (ts,concs) = tsolve(Diff=1.0,fileName="m85.xml",outName="o85.pvd",mode=mode) 
  plt.plot(ts,concs,"b--",label="Diff=1.0, m85",lw=5)
    
  
  plt.legend(loc=4)
  title = "Diffusion profile: " + mode
  plt.title(title)
  plt.ylabel("Conc (at outer bound)") 
  plt.xlabel("time") 
  figname = mode+".png"
  plt.gcf().savefig(figname)


import sys

if __name__ == "__main__":
  import sys
  scriptName= sys.argv[0]
  msg="""
Purpose: 
 
Usage:
"""
  msg+="  %s validation/intact" % (scriptName)
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
    if(arg=="-arg1"):
      arg1=sys.argv[i+1] 




  doit()


