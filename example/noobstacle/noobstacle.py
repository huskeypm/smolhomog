# -*- coding: utf-8 -*-
# <nbformat>2</nbformat>

# <codecell>

# i tink I need to remove all mention of k from everywhere, 
# then it might work 

from dolfin import *
import sys, math, numpy
import numpy as np
from scipy.interpolate import griddata   
import matplotlib.pylab as plt
class empty:pass

tol = 1E-14   # tolerance for coordinate comparisons
EPS=tol
boundsMin = np.zeros(3)
boundsMax = np.zeros(3)

# <codecell>

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
    return False

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

# <codecell>


## computes eff diff coefficient for a variety of cases in which the barrier is defined by PMF, not actual boundary condition
# barrier height is in units E/kT
# discontinuous - use 'discontinuous' PMF [e.g. 0 for free region, non-zero elsewhere]. 
#                 If False, a smoothly varying PMF is used 
# pmfScale - scales Gaussian used for pmf when  discontinuous=False
# fileIn - used with mode='hack2' to import externally-computed potential 
# potential MUST be unitless (e.g. multiply by F/RT before passing into this function)
def doit(dim=1,margin=.1,barrierHeight=50,discontinuous=True,\
         mode="default",pmfScale=1, pmfWidth=0.02, plot=False,outName="img.png",fileIn="",potential="none"):
    ## params
    parms.margin=margin
    parms.update()
    
    ## Define Mesh
    nx = 100; 
    #dim = 1
    if(dim==1):
        mesh = UnitIntervalMesh(nx)
    elif(mode=="hole"):
        #PKH150204 discontinuous=False
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
        mesh = UnitSquareMesh(nx,nx)
    #dim = mesh.ufl_cell().geometric_dimension()

    ## Define a MeshFunction over two subdomains
    subdomains = MeshFunction('size_t', mesh, dim)

    # Mark subdomains with numbers 0 and 1
    markerFree=0;
    markerObstacle=1;
    #for cell_no in range(len(subdomains.array())):
    #    subdomains.array()[cell_no]=markerFree
    print "WARNING: not tested" 
    subdomains.set_all(markerFree)

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
        ## Not MPI safe 
        import numpy
        help = np.asarray(subdomains.array(), dtype=numpy.int32)
        k.vector()[:] = numpy.choose(help, k_values)
        #print "exp(-V): ", k.vector().array()
        kFree = k
        kObs  = k
    
        # PKH 150205
        print "STill not quite correct" 
        #kFree = Constant(1.)  # works, 
	kFree= Constant(k_values[0]) # does not work 
        kObs = Constant(k_values[1])  
        
    else:

        if(potential=="none"): 
          print "Using a continuously-varying PMF" 
          print "WRNING: this might not be correct w MPI, so double check " 
          # exp(-V), where V = exp(-x^2/v)
          expr = Expression("exp(-A*exp(-(pow(x[0]-x0,2) + pow(x[1]-x1,2))/W))",\
                            x0=0.5,x1=0.5,A=pmfScale,W=pmfWidth)
          kObs=expr; print "Do I reall want expr here?"
          kFree=expr
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
          print "WRNING: this likely will not run w MPI"
          print "Using imported grid" 
          pmfFactor = Function(FunctionSpace(mesh,"CG",1))
          print np.shape(pmfFactor.vector()[:] )
          print np.shape(potential.vector()[:] )
          # Higher (less favorable) potentials reduce np.exp toward 0
          # MUST have unitless potential here (e.g. already divided by kT) 
          pmfFactor.vector()[:] = np.exp(-1*potential.vector()[:]) 
          ar = np.array(pmfFactor.vector()[:])
          print "Exp(-p) min %f" % (np.min(ar))
          print "Exp(-p) max %f" % (np.max(ar))
          k = pmfFactor  
    
    
    #plot(subdomains, title='subdomains')
    dx = Measure("dx",domain=mesh)[subdomains] 
    

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
    #form = k*inner((grad(u)+Delta), grad(v))
    #print "CODE IS WRONG - DO NOT USE"
    #formi = inner((grad(u)+Delta), grad(v)) 
    #form = formi*dx(markerFree) + formi*dx(markerObstacle)
    form = kFree*inner((grad(u)+Delta),grad(v))*dx(markerFree) 
    form += kObs*inner((grad(u)+Delta),grad(v))*dx(markerObstacle)
    a = lhs(form)
    L = rhs(form)
    
    ## Compute solution
    u = Function(V)
    solve(a == L, u, bcs)
    File("out.pvd") << u

    ## DEBUG 
    #Vs = FunctionSpace(mesh,"CG",1)
    #up = project(u[0],V=Vs)    
    #ar=np.asarray(up.vector())
    #print "Noobs ", np.min(ar)
    #print "Noobs ", np.max(ar)
    #z = project(k*u,V)
    #File("wrong.pvd") << z
    
    ## Show chi solution  
    if(plot): 
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
          # NOt sure what this does
          # print "Assemble x %f" % (assemble(up*dx(0) + up*dx(1)))#,mesh=mesh))
        
        
      
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
        if(plot): 
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
            if(plot): 
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
            # append form 
            forms.append(grad_Xi_component)
            
    # Evaluate form 
    for i,grad_Xi_component in enumerate(forms):     
        print "WARNING: double check that pmf should be applied to delta"
        form = grad_Xi_component * dx(markerFree) + grad_Xi_component*dx(markerObstacle)  
        integrand = assemble(form)
        print "integrand %d %f" %(i,integrand)
        omegas[i] = integrand
    
    ## Compare diff const. with analytical bounds (2D) 
    diff=boundsMax-boundsMin
    unitCellVol =   np.prod(diff[0:dim])

    freeVol = assemble(Constant(1.)*dx(markerFree))#,mesh=mesh)
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
    if(plot):
      results.interp1 = interp1
      results.gx=gx
    return results
    

## I think I am just exemplifying chi distro, eff diff for different DISCONTINUOUS potentials
def DiscontTest():
	#results=doit(dim=2,discontinuous=True,barrierHeight=-2,mode="hack",plot=True)
	#results=doit(dim=2,discontinuous=True,barrierHeight=-0.5,mode="hack",plot=True)
	results=doit(dim=2,discontinuous=True,barrierHeight=-0.3,mode="hack",plot=True)
	results=doit(dim=2,discontinuous=True,barrierHeight=0,mode="hack",plot=True)
	results=doit(dim=2,discontinuous=True,barrierHeight=15,mode="hack",plot=True)

   	results=doit(dim=2,margin=0.1,plot=True,outName="discontinuousRepulsive.png")

# <codecell>
# for smoothly-varying potentials
def SmoothTest():
  results=doit(dim=2,discontinuous=False,pmfScale=1,plot=True,outName="smoothRepulsive.png")
  results=doit(dim=2,discontinuous=False,pmfScale=-1,plot=True,outName="smoothAttractive.png")


# explicit obstacle + pmf 
def HoleTest():
	# <codecell>

	results=doit(dim=2,discontinuous=True,pmfScale=-1,mode="hole",plot=True,outName="holeAttractive.png")
	#results=doit(dim=2,discontinuous=False,pmfScale=-1,mode="hole",plot=True,outName="holeAttractive.png")
	#results=doit(dim=2,discontinuous=False,pmfScale=0,mode="hole",plot=True,outName="holeNeutral.png")
	#results=doit(dim=2,discontinuous=False,pmfScale=1,mode="hole",plot=True,outName="holeRepulsive.png")

def validationDiscontVsRepulsive():
	margins=np.linspace(0.05,0.48,20);
	Dss = []
	Dsanals=[]
	phis=[]

	for i,margin in enumerate(margins):
	    results = doit(dim=2,margin=margin,plot=False)
	    Dss.append(results.Ds[0])
	    Dsanals.append(results.Dsanal)
	    phis.append(results.phi)
	    
	Dss = np.array(Dss) 
	Dsanals = np.array(Dsanals)
	phis = np.array(phis)

	# <codecell>

	# plot and process
	plt.figure()
	plt.plot(phis,Dss,'k.',label="Predicted")
	plt.plot(phis,Dsanals,'k--',label="HS bound")
	plt.legend(loc=0)
	plt.ylabel("D")
	plt.xlabel("Accessible volume fraction, $\phi$")
	plt.title("Homogenized smoluchowski equation with Box PMF")
	plt.gcf().savefig("validationDiscontinuousRepulsive.png")
	results = empty()
	results.Dss = Dss
	results.Dsanals = Dsanals
	results.phis = phis
	return results

# <codecell>
# gives chi solutions and eff. diff for a given vol frac wrt pmf ampl. 
def pmfAmplitudes(): 
	allResults = []

	kT = 0.6 # [kcal/mol]
	barrierHeights = np.linspace(-5,5,9)
	#barrierHeights = np.array([-6,-1,0,1,2])
	cols=['b-','b--','b-.','b.','k-','r.','r-.','r--','r-']
	margin=0.2

	plt.figure()
	plt.subplot(121)
	Dss = []
	for i,barrierHeight in enumerate(barrierHeights):
	  VkT = barrierHeight/kT
	  results=doit(dim=2,margin=margin,barrierHeight=VkT,plot=True)
	  Dss.append(results.Ds[0])
	  label = "V=%3.1f" % barrierHeight
	  plt.plot(results.gx[:,0],results.interp1,cols[i],label=label)
	  allResults.append(results)
	    

	plt.legend(bbox_to_anchor=(2.5, -.2),ncol=5)
	plt.xlabel("x")
	plt.ylabel("$\chi_x$") 
	plt.title("$\chi$")

	maxHeight = 10
	results=doit(dim=2,margin=margin,barrierHeight=maxHeight,plot=True)
	Dss = np.asarray(Dss)
	plt.subplot(122)
	plt.plot(barrierHeights,Dss,'k-',label="Predicted")
	plt.scatter(0, 1.0, s=80, facecolors='none', edgecolors='k',label="Free")
	plt.scatter(maxHeight, results.phi/(2-results.phi), s=80, facecolors='k', edgecolors='k',label="HS")
	plt.legend(loc=0)
	plt.ylabel("D")
	plt.xlabel("V [kcal/mol]")
	plt.title("Effective Diff")
	plt.gcf().savefig("discontinuousPmfRange.png")
	#plt.title("Homogenized smoluchowski equation with Box PMF")

# <codecell>
# plot deff vs vol frac and pmf
def dEffsSizePMF():
	kT = 0.6 # [kcal/mol]
	#nBarriers=5
	#nMargins=8
	nBarriers=11
	nMargins=8

	barrierHeights = np.linspace(-5,5,nBarriers)
	margins = np.linspace(0,0.48,nMargins)

	Ds = np.zeros([nBarriers,nMargins])
	phis = np.zeros(nMargins)
	for i,barrierHeight in enumerate(barrierHeights):
	  for j,margin in enumerate(margins):
	    VkT = barrierHeight/kT
	    results=doit(dim=2,margin=margin,barrierHeight=VkT,plot=False)
	    Ds[i,j]=results.Ds[0]
	    phis[j] = results.phi
	    

	    

	# <codecell>

	# for the life of me I could not get pcolormesh to display correctly 
	plt.figure()
	#plt.pcolormesh(G,cmap=plt.cm.jet)
	subplot(121)
	plt.pcolormesh(Ds.T)#,cmap=plt.cm.jet)
	plt.colorbar()
	plt.clim([0,2])
	subplot(122)
	#plt.pcolormesh(np.arange(5),np.arange(8),G,cmap=plt.cm.jet)
	#plt.pcolormesh(margins,barrierHeights,Ds,cmap=plt.cm.jet)
	js = np.arange(nMargins)
	js = [0,2,4,5,6]
	ms=['r-','k-','k--','k-.','k.']
	for j,i in enumerate(js):
	  plt.plot(barrierHeights,Ds[:,i],ms[j],linewidth=2,label='$\phi$=%2.1f' % phis[i])

	    
	    
	plt.xlabel("V [kcal/mol]")
	plt.ylabel("D")
	plt.ylim([0,2])
	plt.xlim([-10,10])

	plt.legend(bbox_to_anchor=(2.5, -.2),ncol=5)
	plt.gcf().savefig("dEffsSizePMF.png")


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
      validationDiscontVsRepulsive()   
      
      # PKH - broken, not sure why I would want this anyway 
      # doit(dim=1,barrierHeight=199,plot=True)
      results=doit(dim=2,discontinuous=False,pmfScale=0,mode="hole",plot=True)

    if(arg=="-pmfAmplitudes"): 
      pmfAmplitudes()
   
    if(arg=="-dEffsSizePMF"): 
      dEffsSizePMF()

    if(arg=="-Tests"):
      DiscontTest()
      SmoothTest()
      HoleTest()
    if(arg=="-Hole"):
      HoleTest()


# <codecell>



# <markdowncell>

# Misc messing around 

# <codecell>

#mesh = UnitSquare(100,100)
##V = FunctionSpace(mesh,"CG",1)
#pmf = Function(V) 

#expr = Expression("exp(-(pow(x[0]-x0,2) + pow(x[1]-x1,2))/0.1)",x0=0.5,x1=0.5)
#pmf.interpolate(expr)
#expr = Expression("exp(-sqrt(pow(x[0]-0.5,2) + pow(x[1]-0.5,2)/2.)")

# <codecell>

#pmf(0.0,0.0)

# <codecell>


