#
# Time dependent solver, single species 
# 
from dolfin import *
import numpy as np
import matplotlib.pylab as plt


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


def tsolve(Diff=1.,fileName="m25.xml",outName="output.pvd",mode="pointsource",pmf=1.):
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

  
  ## weak form 
  # weak form of smol is
  # Int[ u*v - u0*v - -dt(e^-*grad(e^+ u)*grad(v))] 
  # let u' = e^+ * u
  # then 
  # Int[ e^-*(u'*v - u0'*v - -dt(grad(u')*grad(v))] 
  dt=15.0   
  T=100 
  expnpmf=Function(V)
  expnpmf.vector()[:] = np.exp(-pmf/0.6)
  RHS = -inner(D*expnpmf*grad(u), grad(q))*dx
  L = (expnpmf*u*q*dx - expnpmf*u0*q*dx - dt*RHS)
  
  # Compute directional derivative about u in the direction of du (Jacobian)
  a = derivative(L, u, du)
  
  
  problem = MyEquation(a,L,bcs)
  solver = NewtonSolver("gmres")         
  file = File(outName, "compressed")
  
  
  t = 0.0
  concs=[]
  ts = []
  while (t < T):
      # advance 
      t0=t
      t += dt
      u0.vector()[:] = u.vector()
      solver.solve(problem, u.vector())

      # remap to u from u' = e^+ * u (see above) 
      up = Function(V)
      u0p = Function(V)
      u0p.vector()[:] = expnpmf.vector()[:] * u0.vector()[:]
      up.vector()[:] = expnpmf.vector()[:] * u.vector()[:]
      file << (up,t) 

      # report on prev iter
      uds = assemble(u0p*ds(rMarker),mesh=mesh)
      area = assemble(Constant(1.)*ds(rMarker),mesh=mesh)
      conc = uds/area
      ts.append(t0)
      concs.append(conc)


  ts = np.asarray(ts)
  concs = np.asarray(concs)

  return (ts,concs)

def valid2():
  mode = "bc"
  mesh = Mesh("m15.xml") 
  V = FunctionSpace(mesh,"CG",1)
  pmf = Function(V) 
  exprN = Expression("             0") # neutral      
  exprA = Expression("-0.1*(x[0]+10)") # attractive 
  exprR = Expression(" 0.1*(x[0]+10)") # repulsive     

  plt.figure()
  # 
  pmf.interpolate(exprN)
  (ts,concs) = tsolve(Diff=1.0,fileName="m15.xml",outName="o15n.pvd",mode=mode,pmf=pmf.vector())
  plt.plot(ts,concs,"k-",label="neutral") 
  # 
  pmf.interpolate(exprA)
  (ts,concs) = tsolve(Diff=1.0,fileName="m15.xml",outName="o15a.pvd",mode=mode,pmf=pmf.vector())
  plt.plot(ts,concs,"b-",label="attractive") 
  # 
  pmf.interpolate(exprR)
  (ts,concs) = tsolve(Diff=1.0,fileName="m15.xml",outName="o15r.pvd",mode=mode,pmf=pmf.vector())
  plt.plot(ts,concs,"r-",label="repulsive") 

  ##
  plt.legend(loc=4)
  title = "Diffusion profile: "
  plt.title(title)
  plt.ylabel("Conc (at outer bound)") 
  plt.xlabel("time") 
  figname = "comp.png"
  plt.gcf().savefig(figname)


  
  

def valid1():
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
    if(arg=="-valid1"):
      valid1()
    if(arg=="-valid2"):
      valid2()


