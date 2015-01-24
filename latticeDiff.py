#
# Time dependent solver, single species 
# 

# Added functionality to pull out 'x' averages 
# validated results against time-indep 


from dolfin import *
import numpy as np
import sys 
sys.path.append("/home/AD/pmke226/sources/cytoplasm/")
import miscutil

class empty:pass 


class LeftBoundary(SubDomain):
  def inside(self,x,on_boundary):
    edge = (np.abs(x[0]- self.mmin[0]) < self.eps)
#    if on_boundary:
#      print "LEFT", x[0], edge, on_boundary
    return on_boundary and edge

class RightBoundary(SubDomain):
  def inside(self,x,on_boundary):
    edge = (np.abs(x[0]- self.mmax[0]) < self.eps)
    #print x[0], edge, on_boundary
    return on_boundary and edge


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
           outName="output.pvd",mode="bc",pmf=0.,  
           nxIncs=4,
           T=1000, printSlice=False,dt=15., debug=False,eps=DOLFIN_EPS):
  # Create mesh and define function space
  if debug:
    mesh = UnitCubeMesh(16,16,16)
  else:
    mesh = Mesh(fileName)     
  V = FunctionSpace(mesh, "Lagrange", 1)
  dim = np.shape(mesh.coordinates())[1]
  
  
  # Define trial and test functions
  du    = TrialFunction(V)
  q  = TestFunction(V)
  # Define functions
  u   = Function(V)  # current solution
  u0  = Function(V)  # solution from previous converged step


  ## mark boundaries 
  subdomains  = MeshFunction("size_t",mesh,dim)
  subsurfaces = MeshFunction("size_t",mesh,dim-1)
  boundary = LeftBoundary()
  boundary.mmin = np.min(mesh.coordinates(),axis=0)
  boundary.eps = eps
  lMarker = 2

  boundary.mark(subsurfaces,lMarker)
  boundary = RightBoundary()
  boundary.mmax = np.max(mesh.coordinates(),axis=0)
  boundary.eps = eps
  rMarker = 3
  boundary.mark(subsurfaces,rMarker)

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

    bc = DirichletBC(V,Constant(1.0),subsurfaces,lMarker)
    bcs.append(bc)

    bc = DirichletBC(V,Constant(0.0),subsurfaces,rMarker)
    bcs.append(bc)


    # visualize BC 
    f = Function(V)
    for i, bc in enumerate(bcs):
      bc.apply(f.vector())
    File("bc.pvd") << f

  # define integrator 
  dx = Measure("dx")[subdomains ]
  ds = Measure("ds")[subsurfaces]

  
  ## weak form [ from (u-u0)/dt = d^2/dx^2 u ]  
  # params 
  D = Constant(Diff) 
  # weak form of smol is
  # Int[ u*v - u0*v - -dt(e^-*grad(e^+ u)*grad(v))] 
  # let u' = e^+ * u
  # then 
  # Int[ e^-*(u'*v - u0'*v - -dt(grad(u')*grad(v))] 
  expnpmf=Function(V)
  expnpmf.vector()[:] = np.exp(-pmf/0.6)
  #print "DISABLED PMF"
  RHS = -inner(D*expnpmf*grad(u), grad(q))*dx
  L = (expnpmf*u*q*dx - expnpmf*u0*q*dx - dt*RHS)
  ##no PMF 
  #RHS = -inner(D*grad(u), grad(q))*dx(domain=mesh)
  #L = u*q*dx - u0*q*dx - dt * RHS

  ## print Test time indep 
  #form = L
  #form = inner(grad(du), grad(q))*dx(domain=mesh)
  #form += Constant(0.)*q*ds
  #a = lhs(form)
  #L = rhs(form)
  #u = Function(V)
  #solve(a == L, u, bcs)
  #File("test.pvd") << u
  #xMids, volConcs = integrateYZ(mesh,subdomains ,u,nIncs=4)
  #quit()

  
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
  
  t = 0.0
  concs=[]
  concsX=[]
  ts = []
  us = []
  while (t < T):
      print "t=%f" %t
      # advance 
      t0=t
      t += dt
      u0.vector()[:] = u.vector()
      solver.solve(problem,u.vector())

      # remap to u from u' = e^+ * u (see above) 
      u0p.vector()[:] = expnpmf.vector()[:] * u0.vector()[:]
      up.vector()[:] = expnpmf.vector()[:] * u.vector()[:]
      file << (up,t) 
      #print "REMOVED exp/pmf"   
      #file << (u,t) 

      
      # store
      if printSlice:
        us.append(miscutil.Store2DMesh(mesh,up))

      xMids, volConcs = integrateYZ(mesh,subdomains ,u,nIncs=nxIncs)
      concsX.append(volConcs)

      # report on prev iter
      #uds = assemble(u0p*ds(rMarker,domain=mesh))#,mesh=mesh)
      #area = assemble(Constant(1.)*ds(rMarker,domain=mesh))#,mesh=mesh)
      #conc = uds/area
      ts.append(t0)
      #concs.append(conc)


  # package 
  ts = np.asarray(ts)
  concs = np.asarray(concs)

  results = empty()
  results.xMids = xMids
  results.xConcs  = np.asarray(concsX)
  print results.xConcs
  #print np.shape(concsX)
  results.us = us 
  results.u_n = up
  results.ts = ts
  results.concs = concs 
  results.mesh = mesh 

  
  return (results)


### MISC ROUTINES 
# need to move this within loop 
# x - the function 
def integrateYZ(mesh,subdomains ,x,nIncs=4):
    coords = mesh.coordinates()
    bounds = np.max(coords,axis=0) - np.min(coords,axis=0)
    xdim = bounds[0]

    xMins = np.arange(nIncs)*xdim/nIncs
    delx = xdim/np.float(nIncs)
    xMaxs = xMins+delx
    xMids = 0.5*(xMaxs+xMins)

    volConcs = np.zeros(nIncs)
    for i in np.arange(nIncs):
      xbounds = [xMins[i],xMaxs[i]]
      volConcs[i] = miscutil.integrated(mesh,subdomains ,x,xbounds)
      print "<Conc> between x=%f/%f: %f " %(xMins[i],xMaxs[i],volConcs[i])

    return xMids,volConcs


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
  
  import matplotlib.pylab as plt
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
  import matplotlib.pylab as plt
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
  tsolve(Diff=1.0,fileName="m15.xml.gz",outName="o15.pvd",mode=mode) 


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
  outName = "out.pvd"
  Diff=1.

  if len(sys.argv) < 2:
      raise RuntimeError(msg)

  debug = False
  fileName = "none"
  T = 1000.
  for i,arg in enumerate(sys.argv):
    if(arg=="-valid1"):
      valid1()
      quit()
    if(arg=="-valid2"):
      valid2()
      quit()
    if(arg=="-valid"):
      valid()
      quit()
    if(arg=="-file"):
      fileName = sys.argv[i+1]
    if(arg=="-outName"):
      outName = sys.argv[i+1]
    if(arg=="-D"):
      Diff= np.float(sys.argv[i+1])
      print "D=%f"%Diff
    if(arg=="-debug"):
      Diff = 0.01
      #tsolve(Diff=Diff,debug=True,T=100,dt=15,eps=0.1) 
      #tsolve(Diff=Diff,debug=True,T=50,dt=1,eps=0.1) 
      tsolve(Diff=Diff,debug=True,T=10,dt=.1,eps=0.1) 
      quit()


  tsolve(Diff=Diff,fileName=fileName,outName=outName,mode="bc",debug=debug,T=T) 
  






