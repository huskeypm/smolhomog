#
# Time dependent solver, single species 
# 
from dolfin import *
import numpy as np
import matplotlib.pylab as plt


root = "example/lattice/"
class empty:pass

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

def steadysolve(Diff=1.,fileName="m25.xml",outName="output.pvd",mode="pointsource",pmf=1.): 
  problem = probdef(Diff,fileName,mode,pmf) 

  ds = problem.ds
  pmf = problem.pmf
  V = problem.V
  #du = problem.du
  u = problem.u
  q = problem.q
  D = problem.D
  bcs= problem.bcs
  
  # borrow other setup
  form = -inner(D*expnpmf*grad(u), grad(q))*dx
  a = lhs(form) 
  L = rhs(form) 

  x = Function(V) 
  solve(a==L,x,bcs=bcs)


  # 2 D: 
  # get steady soln

  # evaluate flux at RHS
  jRHS_avg = 1 
 
  # conc at LHS, LHS
  cLHS = 1

  # del x 
  del_x = 1
  del_c_del_x = (cRHS - cLHS) / del_x 

  # 1 D: 
  # jRHS = Deff * del c / del x 
  Deff = jRHS_avg / del_c_del_x
 

  # get total cell vol
  # minimax

  # get accessible vol 
  # assemble()


  lower=1 
  print "Deff  lower upper " 

  return Deff 

def probdef(Diff=1.,fileName="m25.xml.gz",mode="pointsource",pmf=1.):

  
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

  problem = empty()
  problem.V = V
  problem.ds = ds 
  problem.du = du 
  problem.pmf = pmf 
  problem.bcs = bcs 
  problem.u = u 
  problem.q = q 
  problem.D = D 
  problem.u0 = u0 
  

  return problem


def tsolve(Diff=1.,fileName="m25.xml.gz",outName="output.pvd",mode="pointsource",pmf=1.):
  problem = probdef(Diff,fileName,mode,pmf) 

  ds = problem.ds
  pmf = problem.pmf
  V  = problem.V
  u = problem.u
  du = problem.du
  q = problem.q
  D = problem.D
  u0 = problem.u0
  bcs = problem.bcs
  
  ## weak form 
  # weak form of smol is
  # Int[ u*v - u0*v - -dt(e^-*grad(e^+ u)*grad(v))] 
  # let u' = e^+ * u
  # then 
  # Int[ e^-*(u'*v - u0'*v - -dt(grad(u')*grad(v))] 
  dt=15.0   
  T=1000 
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

def DHExpression(x0=0,x1=0):
  exact=0
  if(exact==1):
    import sys
    sys.path.append("/home/huskeypm/sources/fenics-pb/pete") 
    import poissonboltzmann as pb
    pb.params.center=np.array([x0,x1,0.])
    print "need to adjust DH to use center (power(x[0]-c0,2) "
    pb.params.molrad = 1. 
    print "WARNING: what is radius from original mesh eneration? "
    exprA = pb.DebyeHuckelExpr()

  else:
    exprA = Expression("-3*exp(-(pow(x[0]-x0,2) + pow(x[1]-x1,2))/0.1)",x0=x0,x1=x1)

  return exprA


def valid2():
  mode = "bc"
  mesh = Mesh(root+"m15.xml.gz") 
  V = FunctionSpace(mesh,"CG",1)
  #exprA = Expression("-0.1*(x[0]+10)") # attractive 
  #exprR = Expression(" 0.1*(x[0]+10)") # repulsive     

  ## make pmf maps 
  pmf = Function(V) 
  mask = np.copy(pmf.vector())

  # this scoots around a Debye-Huckel-like potential over all obstactles
  for i in np.arange(8):
    for j in np.arange(8):
      x0 = -7.+2*i
      x1 = -7.+2*j
      exprA = DHExpression(x0,x1)
      pmf.interpolate(exprA)
      mask += pmf.vector()[:]

  pmf.vector()[:] = mask


  ## create pmfs going into attractive, neutral, repulseive cases
  pmfn = Function(V)
  pmfn.vector()[:] = 0.
  pmfa = Function(V)
  pmfa.vector()[:] = pmf.vector()[:] 
  pmfr = Function(V)
  pmfr.vector()[:] = -1*pmf.vector()[:] 
  
  plt.figure()
  (ts,concs) = tsolve(Diff=1.0,fileName=root+"m15.xml.gz",outName="o15n.pvd",mode=mode,pmf=pmfn.vector())
  plt.plot(ts,concs,"k--",label="neutral") 
  # 
  (ts,concs) = tsolve(Diff=1.0,fileName=root+"m15.xml.gz",outName="o15a.pvd",mode=mode,pmf=pmfa.vector())
  plt.plot(ts,concs,"b.",label="attractive") 
  # 
  (ts,concs) = tsolve(Diff=1.0,fileName=root+"m15.xml.gz",outName="o15r.pvd",mode=mode,pmf=pmfr.vector())
  plt.plot(ts,concs,"r.",label="repulsive") 

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
  (ts,concs) = tsolve(Diff=1.0,fileName="m15.xml.gz",outName="o15.pvd",mode=mode) 
  plt.plot(ts,concs,"b-",label="Diff=1.0, m15",lw=1)
  #
  (ts,concs) = tsolve(Diff=1.0,fileName="m25.xml.gz",outName="o25.pvd",mode=mode) 
  plt.plot(ts,concs,"b-",label="Diff=1.0, m25",lw=2)
  #
  (ts,concs) = tsolve(Diff=1.0,fileName="m50.xml.gz",outName="o50.pvd",mode=mode) 
  plt.plot(ts,concs,"b-.",label="Diff=1.0, m50",lw=3)
  #
  (ts,concs) = tsolve(Diff=1.0,fileName="m75.xml.gz",outName="o75.pvd",mode=mode) 
  plt.plot(ts,concs,"b-.",label="Diff=1.0, m75",lw=4)
  #
  (ts,concs) = tsolve(Diff=1.0,fileName="m85.xml.gz",outName="o85.pvd",mode=mode) 
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
  msg+="  %s -valid1/-valid2" % (scriptName)
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
    if(arg=="-valid3"):
      steadysolve(fileName=root+"m85.xml.gz")



