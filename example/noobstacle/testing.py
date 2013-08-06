import matplotlib.pyplot as plt 

from dolfin import *
import sys, math, numpy
import numpy as np
from scipy.interpolate import griddata   

tol = 1E-14   # tolerance for coordinate comparisons
EPS=tol
boundsMin = np.zeros(3)
boundsMax = np.zeros(3)


## Define layered problem over omega in [0,1] where obstacles 
## are represented by large PMF (Omega0, Omega2)
## Note that I am using a dirichlet condition of X=0 on both ends, since 
## the periodic BCs weren't working; results are independent of the value of the dirichlet cond
## so I should get the same result had periodic BCs been used 

class params():
    def __init__(self,margin=0.1):
        self.mid = 0.5
        self.margin = margin
        self.update()
    def update(self):    
        self.minSide = self.mid-self.margin
        self.maxSide = self.mid+self.margin   

parms = params()

    
def DefineBoundary(x,btype,on_boundary):
  if(not on_boundary):
    return 0

  # need this, since for some reason 'dim' is not global
  dim = np.shape(x)[0]

  tb = fb = 0
  lr = ( np.abs(x[0]-boundsMin[0]) < EPS or np.abs(x[0]-boundsMax[0]) < EPS)
  if(dim>=2):
    tb = ( np.abs(x[1]-boundsMin[1]) < EPS or np.abs(x[1]-boundsMax[1]) < EPS)
  if(dim==3):
    fb = ( np.abs(x[2]-boundsMin[2]) < EPS or np.abs(x[2]-boundsMax[2]) < EPS)

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

#print "WARNING: I 'flipped' the obstacle location here" 
class Omega1(SubDomain):
    def inside(self, x, on_boundary):
        dim = np.shape(x)[0]
        if(dim==1):
          return True if (x[0] >= parms.minSide and x[0]<=parms.maxSide) else False
          #return True if (x[0] <= parms.minSide and x[0]>=parms.maxSide) else False
        else:
          return True if (x[0] >= parms.minSide and x[0]<=parms.maxSide \
                          and x[1]>= parms.minSide and x[1]<=parms.maxSide) else False 



# barrier height is in units E/kT
# discontinuous - use 'discontinuous' PMF [e.g. 0 for free region, non-zero elsewhere]. 
#                 If False, a smoothly varying PMF is used 
# pmfScale - scales Gaussian used for pmf when  discontinuous=False
# fileIn - used with mode='hack2' to import externally-computed potential 
def doit(dim=1,margin=.1,barrierHeight=50,discontinuous=True,\
         mode="default",pmfScale=1, pmfWidth=0.02, plot=False,outName="img.png",fileIn="",potential=""):
    ## params
    parms.margin=margin
    parms.update()
    
    print mode 
    ## Define Mesh
    nx = 100; 
    #dim = 1
    if(dim==1):
        mesh = UnitInterval(nx)
    elif(mode=="hole"):
        discontinuous=False
        mesh = Mesh("cylinder.xml") 
    elif(mode=="holeBig"):
        discontinuous=False
        mesh = Mesh("cylinderBig.xml")
    elif(mode=="hack"):
        mesh = Mesh("cylinderDense.xml") 
        discontinuous=True
    elif(mode=="hack2"): 
        mesh = Mesh(fileIn)
        discontinuous = False 
    else:
        mesh = UnitSquare(nx,nx)
    #dim = mesh.ufl_cell().geometric_dimension()

    ## Define a MeshFunction over two subdomains
    subdomains = MeshFunction('uint', mesh, dim)

    # Mark subdomains with numbers 0 and 1
    markerFree=0;
    markerObstacle=1;
    for cell_no in range(len(subdomains.array())):
        subdomains.array()[cell_no]=markerFree

    if(discontinuous): 
        #subdomain0 = Omega0() 
        #subdomain0.mark(subdomains, 0)
        subdomain1 = Omega1() # obstacle
        subdomain1.mark(subdomains, markerObstacle)
        
        
        ## Define diffusion constant (D = exp(-V)) within two domains 
        # Discontinuous basis is used for boundary condition
        # LAter, continuous basis used for PDE solution 
        V0 = FunctionSpace(mesh, 'DG', 0)
        k = Function(V0)
        
        ## Define diffusion constant to assume large values at 'edges' (Omega0 and Omega2), and 
        ## small value at middle (Omega1). We assume D(y) = exp(-V(y)), so large D corresponds to a obstacle
        Vpmf = np.array([0,barrierHeight]) # pmg 
        expnv = np.exp(-Vpmf)
        k_values = expnv
        
        # populate 
        import numpy
        help = numpy.asarray(subdomains.array(), dtype=numpy.int32)
        k.vector()[:] = numpy.choose(help, k_values)
        print "exp(-V): ", k.vector().array()
    else:

        if(potential=="none"): 
          print "Using a continuously-varying PMF" 
          # exp(-V), where V = exp(-x^2/v)
          expr = Expression("exp(-A*exp(-(pow(x[0]-x0,2) + pow(x[1]-x1,2))/W))",\
                            x0=0.5,x1=0.5,A=pmfScale,W=pmfWidth)
          k=expr
          pmf = Function(FunctionSpace(mesh,"CG",1))
          pmf.interpolate(k)
          File("test.pvd") << pmf
          gx,gy = np.mgrid[0:1:100j,0:1:100j]
          interp0 = griddata(mesh.coordinates(),pmf.vector(),(gx,gy))
          plt.figure()
          plt.pcolormesh(interp0.T)
          plt.colorbar()
          plt.title("pmf") 
        else:
          print "Using imported grid" 
          print "WARNING: wrong beta"
          beta = 1/0.6
          pmf = Function(FunctionSpace(mesh,"CG",1))
          print np.shape(pmf.vector()[:] )
          print np.shape(potential.vector()[:] )
          pmf.vector()[:] = np.exp(-beta * potential.vector()[:]) 
          k = pmf  
    
    
    #plot(subdomains, title='subdomains')
    dx = Measure("dx")[subdomains] 
    

    ## Assign BC 
    for i in np.arange( dim ):
        boundsMin[i] = np.min(mesh.coordinates()[:,i])
        boundsMax[i] = np.max(mesh.coordinates()[:,i])

    bcs=[]
    if(dim==1):
        V = FunctionSpace(mesh, 'Lagrange', 1)
        bcs.append(DirichletBC(V,Constant(0),XBoundary()))
    else:
        V = VectorFunctionSpace(mesh, 'Lagrange', 1)
        bcs.append(DirichletBC(V.sub(0),Constant(0),XBoundary()))
        bcs.append(DirichletBC(V.sub(1),Constant(0),YBoundary()))

    
    ## Define weak form of homogenized smol eqn 
    D0    = Constant(1.)    # WARNING: we ignore Dij for these calcs
    Dij = D0*Identity( dim )    
    Delta = Identity( dim )
    
    #  D [del u + 1] del v = 0
    u = TrialFunction(V)
    v = TestFunction(V)
    # leaving out Dij=const, since solving steady state 
    formi = k*inner((grad(u)+Delta), grad(v)) 
    #print "CODE IS WRONG - DO NOT USE"
    #formi = inner((grad(u)+Delta), grad(v)) 
    
    form = formi*dx(markerFree) + formi*dx(markerObstacle)
    a = lhs(form)
    L = rhs(form)
    
    ## Compute solution
    u = Function(V)
    solve(a == L, u, bcs)
    File("out.pvd") << u
    z = project(k*u,V)
    File("wrong.pvd") << z
    
    ## Show chi solution  
    
    if(dim==1):
        gx = np.mgrid[0:1:100j]
        interp0 = griddata(mesh.coordinates(),u.vector(),(gx))
        interp1 = interp0
        
    else:
        Vs = FunctionSpace(mesh,"CG",1)
        up = project(u[0],V=Vs)    
       
        gx,gy = np.mgrid[0:1:100j,0.5:0.5:1j]
        interp1 = griddata(mesh.coordinates(),up.vector(),(gx,gy))
        gx,gy = np.mgrid[0:1:100j,0:1:100j]
        interp0 = griddata(mesh.coordinates(),up.vector(),(gx,gy))
        
        
        up = project(u[1],V=Vs) 
        interpx = griddata(mesh.coordinates(),up.vector(),(gx,gy))
        plt.figure()
        plt.subplot(421)
        plt.pcolormesh(interp0.T)
        plt.title("u0") 
        plt.subplot(422)
        plt.pcolormesh(interpx.T)
        #def doplt(up):
        #    interpgx = griddata(mesh.coordinates(),up.vector(),(gx,gy))
        #    plt.pcolormesh(interpgx.T)
        plt.title("u1")
            
        _component = grad(u[0])
        up = project(_component[0],V=Vs) 
        interpgx = griddata(mesh.coordinates(),up.vector(),(gx,gy))
        plt.subplot(423)
        plt.pcolormesh(interpgx.T)
        plt.title("du0")       
        _component = grad(u[1])
        up = project(_component[1],V=Vs) 
        interpgx = griddata(mesh.coordinates(),up.vector(),(gx,gy))
        plt.subplot(424)
        plt.pcolormesh(interpgx.T)
        plt.title("du1")     
        
        
        _component = grad(u[0])
        up = project(k*(_component[0]),V=Vs) 
        interpgx = griddata(mesh.coordinates(),up.vector(),(gx,gy))
        plt.subplot(425)
        plt.pcolormesh(interpgx.T)
        plt.title("du0")    
        
              
        _component = grad(u[0])
        up = project(k*(_component[0]+1.),V=Vs) 
        interpgx = griddata(mesh.coordinates(),up.vector(),(gx,gy))
        plt.subplot(427)
        plt.pcolormesh(interpgx.T)
        plt.title("k*du0+1)")    
        print "Assemble x %f" % (assemble(up*dx(0) + up*dx(1),mesh=mesh))
        
        
      
    if(plot):     
        plt.figure()
        plt.subplot(121)  
        if(dim==2):
            plt.pcolormesh(gx[:,0],gy[0,:],interp0.T)
            plt.colorbar()
        plt.subplot(122)  
        if(dim==2):
          plt.plot(gx[:,0],interp1)
        else:
          plt.plot(gx,interp1)
        plt.gcf().savefig(outName)

 
    ## Compute effective diff. tensor     
    # define form for grad(chi) + delta 
    
    omegas = np.zeros(dim)    
    forms = []
    if(dim==1):    
        grad_Xi_component = k*(grad(u)+Constant(1.0))  
        forms.append(grad_Xi_component)
        
        _component = grad(u)
        g=project(_component,FunctionSpace(mesh,"CG",1))
        interp2 = griddata(mesh.coordinates(),g.vector(),(gx))
        plt.figure()
        plt.plot(interp2,label="grad")
        plt.plot(interp1,label="sol")
        plt.legend()
 
    
    if(dim==2):
      for i in range(dim):
            v = [0,0,0]
            v[i] = 1
            grad_Xi_component = k*(inner(grad(u[i]),Constant((v[0],v[1]))) + Constant(1))
            _component = (inner(grad(u[i]),Constant((v[0],v[1]))))
            _component = grad(u[i])
            g=project(_component[1],FunctionSpace(mesh,"CG",1))
            interp2 = griddata(mesh.coordinates(),g.vector(),(gx,gy))
            print "inerpmin %f" % (np.min(interp2[:,50]))
            if(i==1):
              plt.figure()
              plt.pcolormesh(interp2.T)
              plt.colorbar()
              File("comp.pvd") << g
              plt.figure()
              gx,gy = np.mgrid[0:1:100j,0.5:0.5:1j]
              interp1 = griddata(mesh.coordinates(),up.vector(),(gx,gy))               
              plt.plot(gx[:,0],interp1,label="sol")
              interp2 = griddata(mesh.coordinates(),g.vector(),(gx,gy)) 
              plt.plot(gx[:,0],interp2,label="grad")
              plt.legend()  
            forms.append(grad_Xi_component)
            
    # Evaluate form 
    for i,grad_Xi_component in enumerate(forms):     
        print "WARNING: double check that pmf should be applied to delta"
        form = grad_Xi_component * dx(markerFree) + grad_Xi_component*dx(markerObstacle)  
        integrand = assemble(form)
        omegas[i] = integrand
    
    ## Compare diff const. with analytical bounds (2D) 
    diff=boundsMax-boundsMin
    unitCellVol =   np.prod(diff[0:dim])
    freeVol = assemble(Constant(1.)*dx(markerFree),mesh=mesh)
    phi = freeVol/unitCellVol
    Dsanal = phi/(2-phi)
    Ds = omegas/unitCellVol
    print phi
    print "Ds anal est (2D) ", Dsanal
    print "Ds pred ", Ds
    class empty:pass
    results = empty()
    results.phi = phi
    results.Ds = Ds
    results.Dsanal= Dsanal
    results.interp1 = interp1
    results.gx=gx
    return results
    

   
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
      #arg1=sys.argv[i+1] 
      doit(dim=1,barrierHeight=199,plot=True)
      results=doit(dim=2,discontinuous=False,pmfScale=0,mode="hole",plot=True)




  #:doit(fileIn)


