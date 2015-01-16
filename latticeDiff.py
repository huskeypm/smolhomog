#
# Time dependent solver, single species 
# 
from dolfin import *
import numpy as np
import matplotlib.pylab as plt
from scipy.interpolate import griddata
class empty:pass 


class LeftBoundary(SubDomain):
  def inside(self,x,on_boundary):
    edge = (np.abs(x[0]- self.mmin[0]) < DOLFIN_EPS) 
    #print x[0], edge, on_boundary
    return on_boundary and edge

class RightBoundary(SubDomain):
  def inside(self,x,on_boundary):
    edge = (np.abs(x[0]- self.mmax[0]) < DOLFIN_EPS) 
    #print x[0], edge, on_boundary
    return on_boundary and edge

def PrintSlice(mesh,u):  
#    mesh = results.mesh
    #dims = np.max(mesh.coordinates(),axis=0) - np.min(mesh.coordinates(),axis=0)
    mmin = np.min(mesh.coordinates(),axis=0)
    mmax = np.max(mesh.coordinates(),axis=0)
    
    #u = results.u_n.split()[0]
    #u = results.u_n
    up = project(u,FunctionSpace(mesh,"CG",1))
    res = 100
    #(gx,gy,gz) = np.mgrid[0:dims[0]:(res*1j),
    #                      dims[1]/2.:dims[1]/2.:1j,
    #                      0:dims[2]:(res*1j)]
    (gx,gy) = np.mgrid[mmin[0]:mmax[0]:(res*1j),
                       mmin[0]:mmax[1]:(res*1j)]
    img0 = griddata(mesh.coordinates(),up.vector(),(gx,gy))
    return img0


# Nonlinear equation 
class MyEquation(NonlinearProblem):
    def __init__(self, a, L,bcs):
        NonlinearProblem.__init__(self)
        self.L = L
        self.a = a
        self.bcs = bcs
        #self.reset_sparsity = True
    def F(self, b, x):
        assemble(self.L, tensor=b)
        for bc in self.bcs:
          bc.apply(b,x)
    def J(self, A, x):
        assemble(self.a, tensor=A)#, reset_sparsity=self.reset_sparsity)
        for bc in self.bcs:
          bc.apply(A)
        #self.reset_sparsity = False


def tsolve(Diff=1.,fileName="m25.xml.gz",\
           outName="output.pvd",mode="pointsource",pmf=0.,  
           T=1000):
  # Create mesh and define function space
  mesh = Mesh(fileName)     
  V = FunctionSpace(mesh, "Lagrange", 1)
  
  
  # Define trial and test functions
  du    = TrialFunction(V)
  q  = TestFunction(V)
  # Define functions
  u   = Function(V)  # current solution
  u0  = Function(V)  # solution from previous converged step



  ## mark boundaries 
  subdomains = MeshFunction("size_t",mesh,1)
  boundary = LeftBoundary()
  boundary.mmin = np.min(mesh.coordinates(),axis=0)
  lMarker = 2
  boundary.mark(subdomains,lMarker)
  boundary = RightBoundary()
  boundary.mmax = np.max(mesh.coordinates(),axis=0)
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
    #u.vector()[np.where(u.vector()) < vmin] = vmin
    #u.vector()[np.where(u.vector()) > vmax] = vmax
    u0.vector()[:] = u.vector()[:]
  
    print "total conc", assemble(u*dx)

  ## dirichlet on left 
  else: 
    ic = Expression("0")
    u.interpolate(ic)
    u0.interpolate(ic)

    bc = DirichletBC(V,Constant(1.0),subdomains,lMarker)
    bcs.append(bc)
   
    # visualize BC 
    f = Function(V)
    bc.apply(f.vector())
    File("bc.pvd") << f

  # define integrator 
  ds = Measure("ds")[subdomains]

  
  ## weak form 
  # params 
  D = Constant(Diff) 
  t = 0.0
  dt=15.0   
  # weak form of smol is
  # Int[ u*v - u0*v - -dt(e^-*grad(e^+ u)*grad(v))] 
  # let u' = e^+ * u
  # then 
  # Int[ e^-*(u'*v - u0'*v - -dt(grad(u')*grad(v))] 
  expnpmf=Function(V)
  expnpmf.vector()[:] = np.exp(-pmf/0.6)
  RHS = -inner(D*expnpmf*grad(u), grad(q))*dx
  L = (expnpmf*u*q*dx - expnpmf*u0*q*dx - dt*RHS)
  #no PMF RHS = -inner(D*grad(u), grad(q))*dx
  #no PMF L = u*q*dx - u0*q*dx - dt * RHS

  
  # Compute directional derivative about u in the direction of du (Jacobian)
  a = derivative(L, u, du)
  
  
  ## define solver,params 
  problem = MyEquation(a,L,bcs)
  solver = NewtonSolver()                
  solver.parameters["linear_solver"] = "gmres"
  solver.parameters["convergence_criterion"] = "incremental"
  solver.parameters["relative_tolerance"] = 1e-6
  file = File(outName, "compressed")
  
  ## Variables for storing answers (after projection) 
  # need to declare outside of loop 
  up = Function(V)
  u0p = Function(V)
  
  concs=[]
  ts = []
  us = []
  while (t < T):
      # advance 
      t0=t
      t += dt
      u0.vector()[:] = u.vector()
      solver.solve(problem,u.vector())

      # remap to u from u' = e^+ * u (see above) 
      u0p.vector()[:] = expnpmf.vector()[:] * u0.vector()[:]
      up.vector()[:] = expnpmf.vector()[:] * u.vector()[:]
      file << (up,t) 

      
      # store
      us.append(PrintSlice(mesh,up))

      # report on prev iter
      uds = assemble(u0p*ds(rMarker,domain=mesh))#,mesh=mesh)
      area = assemble(Constant(1.)*ds(rMarker,domain=mesh))#,mesh=mesh)
      conc = uds/area
      ts.append(t0)
      concs.append(conc)


  ts = np.asarray(ts)
  concs = np.asarray(concs)

  results = empty()
  results.us = us 
  results.u_n = up
  results.ts = ts
  results.concs = concs 
  results.mesh = mesh 

  
  return (results)

def DHExpression(x0=0,x1=0):
  exact=0
  if(exact==1):
    import sys
    sys.path.append("/home/huskeypm/sources/fenics-pb/pete") 
    import poissonboltzmann as pb
    pb.params.center=np.array([x0,x1,0.]) 
    print """need to adjust DH to use center (power(x[0]-c0,2) 
    pb.params.molrad = what is radius from original mesh eneration? 
    """
    exprA = pb.DebyeHuckelExpr()

  else:
    exprA = Expression("-3*exp(-(pow(x[0]-x0,2) + pow(x[1]-x1,2))/0.1)",x0=x0,x1=x1)

  return exprA


def valid2():
  mode = "bc"
  mesh = Mesh("m15.xml.gz") 
  V = FunctionSpace(mesh,"CG",1)
  #exprA = Expression("-0.1*(x[0]+10)") # attractive 
  #exprR = Expression(" 0.1*(x[0]+10)") # repulsive     

  ## make pmf maps 
  pmf = Function(V) 
  #if MPI.rank(mpi_comm_world())==0:
  if 1: # need to find away to perform this via MPI - np might not be appropriate 
    print "WARNING: not sure if this is being applied on all processors!"
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
  (ts,concs) = tsolve(Diff=1.0,fileName="m15.xml.gz",outName="o15n.pvd",mode=mode,pmf=pmfn.vector())
  plt.plot(ts,concs,"k--",label="neutral") 
  # 
  (ts,concs) = tsolve(Diff=1.0,fileName="m15.xml.gz",outName="o15a.pvd",mode=mode,pmf=pmfa.vector())
  plt.plot(ts,concs,"b.",label="attractive") 
  # 
  (ts,concs) = tsolve(Diff=1.0,fileName="m15.xml.gz",outName="o15r.pvd",mode=mode,pmf=pmfr.vector())
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

def valid():
  mode = "bc"
  (ts,concs) = tsolve(Diff=1.0,fileName="m15.xml.gz",outName="o15.pvd",mode=mode) 


import sys

if __name__ == "__main__":
  import sys
  scriptName= sys.argv[0]
  msg="""
Purpose: 
 
Usage:
"""
  msg+="  %s -valid/-valid1/-valid2\n" % (scriptName)
  msg+="or \n"
  msg+="  %s -test filename.xml" % (scriptName)
  msg+="""
  
 
Notes:

"""
  remap = "none"

  if len(sys.argv) < 2:
      raise RuntimeError(msg)

  for i,arg in enumerate(sys.argv):
    if(arg=="-valid1"):
      valid1()
    if(arg=="-valid2"):
      valid2()
    if(arg=="-valid"):
      valid()
    if(arg=="-test"):
      fileName = sys.argv[i+1]
      tsolve(Diff=1.0,fileName=fileName,outName="out.pvd",mode="bc") 
  






