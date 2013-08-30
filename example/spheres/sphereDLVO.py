## worked example 14.5, Israelachvilie Surface Forces book
##
## PARAMS 
##
# W(D) = 0.5 R Z exp(-kappa D) - A*R/12D
print "WARNING: NEED TO DOUBLE CHECK THAT THIS IS THE CORRECT VDW TERM"
print "WARNING: DLVO stuff is for LIKE charged surface, not opposite"

import sys
prefix = "/net/home/huskeypm/Sources/"
prefix = "/home/huskeypm/sources/"
sys.path.append(prefix+"/smolhomog/example/")
sys.path.append(prefix+"/modified-pb/example/")
sys.path.append("/net/home/huskeypm/Sources/homogenization/example/volfracs/")
sys.path.append("/net/home/huskeypm/bin/grids/")
#import createVolumeFractionMeshes as cvm
path = prefix+"/homogenization/example/volfracs/"
import matplotlib.pylab as plt
import numpy as np
import poissonboltzmann as pb
import homoglight as hl
import caseRunner as cr
debug = 0
cr.dims = 3


J_to_kT = 1/4.114e-21 # kT per J at 298 based on Iraelachvili
kT_to_J = 1/J_to_kT
nm_to_m = 1e-9
m_to_nm=1/nm_to_m
nm_to_Ang = 10.
Ang_to_m = 1e-10
m_to_Ang = 1/Ang_to_m
Ang_to_nm = 1e-1

parms = pb.parms
print "WARNING: this is a bit of a debugging hack to prevent divide by zeros/exploding potential "
parms.dtol = 0.5e-10 *m_to_nm   
#mine
#ionC = 0.1 # [M] 
#iKappa =0.304/np.sqrt(ionC); kappa = 1/iKappa

# Might be a good idea to move hmogenization/exmaple/volcracs into smolhomog
# or even just call gamer on the fly

# <codecell>

def mergeplot():
  fig = plt.figure()
  #axplot = fig.add_axes([0.07,0.25,0.90,0.70])
  #axplot.
  plt.plot(np.random.randn(100))
  
  #axicon = fig.add_axes([0.07+0.11*k,0.05,0.1,0.1])
  axicon = fig.add_axes([0.1,0.1,0.3,0.3])
  axicon.imshow(np.zeros((2,2)),interpolation='nearest')
  axicon.set_xticks([])
  axicon.set_yticks([])

# <markdowncell>

# DLVO section 

# <codecell>


#Dnm= 0.25
# A = 1e-20 # Hamaker [J]
# R=0.1e-6; # particle diam, [m]
# minDLVO [J] minimum value returned for W 
def DLVO(Dnm,psi0=24.5,z=1,A = 1e-20, R=0.1e-6,minDLVO=-3*kT_to_J):
  D = Dnm*nm_to_m # D [m]
  kappa = parms.kappa*nm_to_Ang # 1/[A] --> 1/[nm]  

  absz = np.abs(z)
  if(z==0): 
    signz=0.
  else:
    signz=np.sign(psi0/z)  # THIS IS NOT CORRECT!!!! (sign=1 like charges, sign=-1 opposite)
  # if validating against Israelach, use 9.38e-11 and psi0/107 for T=38
  Z = signz*9.22e-11 *np.tanh(absz*psi0/103)**2   # [J/m] at T = 25 C 
  #Z = signz*9.38e-11 *np.tanh(absz*psi0/107)**2   # [J/m] at T = 38 C 
  W1 = 0.5*R*Z*np.exp(-kappa*Dnm)

  # NEED TO DOUBLE CHECK THAT THIS IS THE CORRECT VDW TERM 
  #W2 = A*R/12e-9 #/(12*D)
  W2 = A*R/(12*D)

  #W1 = 0.   # for VDW only 
  #W2 = 0.   # for electro only 
  W = W1-W2
  #print W1,W2,W
    
  # add fake repulsive  
  W = np.max((W,minDLVO)) 
  return W

vW=np.vectorize(DLVO,otypes=[float])


# For validating DLVO model against Israelachvili
# Note that the model is NOT valid for oppositely charged
# surfaces, so this is more of an illustration 
def validateDLVO():
  Dnm = np.linspace(0.1,8,80) # [nm]
  A = 1e-20 # Hamaker [J]  for vesicles (Israelachvili)
  R=0.1e-6; # particle diam, [m]
  parms.ionC = 0.1 # [M]
  parms.update()

  #parms.kappa = 1/0.95
  israelachVal=vW(1.0,14.5,1,1e-20, 0.1e-6,-1.5e19)      
  israelachVal*=1e20
  assert(np.abs(israelachVal - -5.137)<0.001), "DLVO broken %f" % israelachVal
    
  runDLVO(A=A,R=R,minDLVO=-1.5e-19,name="dlvo_valid.png")
  print "Compare dlvo_valid.png with Fig 14.14 of Israelachvili"
    

# Evaluates DLVO potential at several surface potentials from Israelachvili
def runDLVO(A=1e-20,R=0.1e-6,minDLVO=-2*kT_to_J,name="dlvo.png"):    
  Dnm = np.linspace(0.1,8,80) # [nm]
  ## get sigmas from potentials so we can
  # adjust poetntial based on ionic str.
  israelsPsi0s =np.array([0,14.5,24.5,34.5])
  israelsPsi0s =np.array([14.5,24.5,34.5])
  # Graham equation
  sigmas= 0.117*np.sqrt(parms.ionC)*np.sinh(israelsPsi0s/51.4)

  fig = plt.figure()
  ax = fig.add_subplot(111)

  units = "kT"
  #units = "J"  
  cols=['r','b']
  styles=['-','--','-.','.']
  #for j,z in enumerate([-1,1]):
  for j,z in enumerate([1,-1]):
    for i,sigma in enumerate(sigmas):
      psi0=pb.Grahame(sigma,parms.ionC)
        
      if(units=="J"):   
            #Dnm=0.95
            #wa=DLVO(Dnm,psi0=34.5,z=1,A = 1e-20, R=0.1e-6,minDLVO=-50*kT_to_J)
            #print psi0,z,A,R,minDLVO
            vals =vW(Dnm,psi0,z,A,R,minDLVO)
            
      else:
            vals =vW(Dnm,psi0,z,A,R,minDLVO) * J_to_kT
      ax.plot(Dnm,vals,cols[j]+styles[i],\
              label="%4.1f [mV], z=%d" %(psi0,z))

  if(units=="J"):            
    ax.set_ylim([-2.e-19,1.5e-19])
    ax.set_ylabel("W [J]")
  else:
    ax.set_ylabel("W [kT]")
            
  ax.set_xlabel("D [nm]")
  plt.legend(loc=1)

  #ax1=twinx(ax)  
  #for j,z in enumerate([-1,1]):
  #  for i,sigma in enumerate(sigmas):
  #    psi0=pb.Grahame(sigma,ionC)
  #    ax1.plot(Dnm,vW(Dnm,psi0,z)*J_to_kT,cols[j]+styles[i],\
  #            label="%4.1f [mV], z=%d" %(psi0,z))

  plt.title("DLVO (R=%5.1f [nm], A=%5.1f [kT])" % (R*m_to_nm,A*J_to_kT))
  #ax.xaxis.set_ticks([1.,2.,3.,10.])
  ax.grid(True)
  plt.gcf().savefig(name)

# <codecell>


# <markdowncell>

# Interpolation section 

# <codecell>

# most simple way of implementing this DLVO expression might be to 
# use Expression to compute distance from sphere/plane, then pass 
# the result to the vectorized DLVO function


# psi0 [mV]  - potential at boundary 
def ApplyDLVO(case="unitsphere",mesh="none",psi0=24.5): # Israelachvili ,psi0 = 24.5,  z = 1.,R=1e-7,A=1e-20):
# R - particle radius [m]
# A - Hamakaer [J] 
  ## Decide on coordinates, distance expression
  # assuming [A] here 
  if(case=="unitline"):
    mesh = UnitInterval(100)
    mesh.coordinates()[:] = mesh.coordinates()[:]*8*nm_to_Ang
    dexpr = Expression("x[0]-px",px=-0.1*nm_to_Ang)

  if(case=="unitplane"):
    mesh = UnitSquare(100,100)
    cs = mesh.coordinates()[:]
    cs[:,0]*= 20.*nm_to_Ang
    cs[:,1]*= 8*nm_to_Ang

    mesh.coordinates()[:] = cs
    dexpr = Expression("x[1]-py",py=-0.1*nm_to_Ang)
    print "WARNING: coded DLVO expression is NOT correct for planes!!"

  if(case=="unitsphere"):
    mesh = UnitCube(20,20,20)
    cs = mesh.coordinates()[:]
    cs *= 2*8.*nm_to_Ang
    cs -= np.array([8,8,8])*nm_to_Ang

    mesh.coordinates()[:] = cs
    sphereCenter=np.array([0,0,0])*nm_to_Ang
    dexpr = Expression("sqrt(pow(x[0]-cx,2)+pow(x[1]-cy,2)+pow(x[2]-cz,2))",\
                        cx=sphereCenter[0],cy=sphereCenter[1],cz=sphereCenter[2])

  # for actual geometry         
  if(case=="sphere"):
    sphereCenter=np.array([0.,0.,0.])
    #sphereCenter=np.array([-21.7,-21.7,-21.7])
    # Get distance from sphere center, then substrate sphere radius to get
    # distance from surface
    print parms.molRad
    dexpr = Expression("sqrt(pow(x[0]-cx,2)+pow(x[1]-cy,2)+pow(x[2]-cz,2))-R",\
                        cx=sphereCenter[0],cy=sphereCenter[1],cz=sphereCenter[2],R= parms.molRad)   
    boxDiam = np.max(mesh.coordinates())-np.min(mesh.coordinates())
    left = np.max(mesh.coordinates())
    #print "Dexpr ", dexpr(left,0,0)
    #print "Dexpr ", dexpr(left,left,left)
    #print "Dexpr ", dexpr(parms.molRad,0,0)
    #print "BoxRad %f R %f " % (boxDiam/2., parms.molRad)
    #quit()

  # for actual geometry         
  if(case=="layer"):
    yt = 1 # y at top [A] 
    yb = 0 # y at bottom         
    dexpr = Expression("yt - x[1]",yt=yt)
    dexpr2 = Expression("x[1]-yb",yb)        
        
  # Get distances from mesh
  #print np.min(mesh.coordinates())
  #print np.max(mesh.coordinates())
  V = FunctionSpace(mesh,"CG",1)
  ds = interpolate(dexpr,V).vector().array()        # in [A]
  #print np.max(ds)
  ds /= nm_to_Ang                                   # in [nm] 
  #print np.max(ds)
  
  # we do an additional step for layered, sine there are two boundaries       
  if(case=="layer"): 
    ds2 = interpolate(dexpr2,V).vector().array()        # in [A]  
    ds2 /= nm_to_Ang                                    # in [nm] 
    ds = np.min((ds,ds2),axis=0)
    
  if(0):  
    d = Function(V)
    d.vector()[:] = ds*nm_to_Ang
    File("distances.pvd") << d
        
  #
  ds[ np.where(ds < parms.dtol) ] = parms.dtol
  plt.plot(ds,np.zeros(np.shape(ds)[0]),"k.")
  plt.ylim([-3,3])
  plt.xlabel("[nm]")
  # evaluate DLVO expression
  #ws = vW(ds.vector()[:],psi0,z)
  #print "z=",z
  #print "psi0",psi0  
  #z *= -1
  #psi0 *= -1
 
  R_m = parms.molRad*Ang_to_m
  ws = vW(ds,psi0,parms.zLig,parms.A,R_m) * J_to_kT
  pmf = Function(V)
  pmf.vector()[:] = ws
  plt.plot(ds,ws,"r.")
  plt.gcf().savefig("distances.png") 
  #print np.sort(ds)[0:10]
  #quit()


  ds*= nm_to_Ang
  return (ds,pmf)

from dolfin import *


# molRad - radius of enzeyme 
# z - enzyme charge
# q - substrate charge 
# A - Hamaker const [J] 
# R - molRad [A] 

def CalcPMF(mesh,meshType="dolfin",pmfType="DebyeHuckel",case="sphere"):               
   pb.parms.sphericalDomainBoundary=False # i think this is for the outer domain??
   pb.parms.update()
   V=FunctionSpace(mesh,"CG",1) 
    
   #
   if(pmfType=="DebyeHuckel"):  
     (Vdumm,psi) = pb.SolvePoissonBoltzmann(mesh,meshType=meshType)
     pmf = Function(V)
   
     # convert from [mV] to [kT] 
     pmf.vector()[:]= psi.vector()[:]*pb.parms.F_o_RT*pb.parms.kT
     pmf.vector()[:]= psi.vector()[:]*zLig
     pmfar = np.asarray(pmf.vector()[:])
     print "Pmf min/max %f/%f [kT]" %(np.min(pmfar),np.max(pmfar))
     (ds,dummypmf)= ApplyDLVO(case=case,mesh=mesh)#,R=R*Ang_to_m) 
   
   elif(pmfType=="DLVO"): 
     #f = pb.DLVOExpr(dim=3)
     #f.interpolate(pmf.vector()[:])
     #pmf = interpolate(f,V)       
     expr = pb.DebyeHuckelExpr(dim=3)
     # get potential at 'left' edge of sphere
     psi0 = expr(parms.molRad,0.,0.)
     (ds,pmf)= ApplyDLVO(case=case,mesh=mesh,psi0=psi0)#,psi0=psi0,z=zLig,R=R*Ang_to_m,A=A)


   ## for debugging PMFs
   if(0):
     print "compare", np.min(ds), np.max(ds)
     # this is a bit of a hack to get the distances from the ApplyDLVO function 
     # plot
     plt.figure()
     #plt.plot(ds.vector().array(),pmf.vector().array(),'k.')
     plt.plot(ds*Ang_to_nm,pmf.vector().array(),'k.')
     plt.title("%s interaction energy \n(R=%4.2f [nm], H=%4.2f [kT],psi=%4.2f [mV],z=%d)" % (pmfType,\
       parms.molRad*Ang_to_nm,parms.A*J_to_kT,psi0,parms.zLig))
     plt.xlabel("Dist from surface[nm]") 
     #plt.xlim([0,8])
     plt.ylim([-7,1])     
     plt.grid(True)
     plt.ylabel("Energy [kT]") 
     plt.gcf().savefig("pmfWRTDist.png") 
     File("testpmf.pvd") << pmf
     quit()

   return pmf
            
# A Hamaker [J] 
#def call(meshName,zLig,molRad,zProt,ionC=0.15,A=1e-20,molGamer=0,debug=0,pmfType="DLVO"):
def call(meshName,molGamer=0,debug=0,pmfType="DLVO"):

  # mesh 
  mesh = Mesh(meshName)

  # recaling all meshes to be much larger, otherwise PMF is wayy too attractive
  #mesh.coordinates()[:] *= 10.


  if(debug): 
    mesh = UnitCube(20,20,20)
    cs = mesh.coordinates()[:]
    cs *= 2*8.*nm_to_Ang
    cs -= np.array([8,8,8])*nm_to_Ang
    mesh.coordinates()[:] = cs


  #  potential of mean force [kT] 
  #psi = CalcPMF(mesh,molRad,z,q,meshType="gamer")
  #print pb.parms.zProt 
  #parms.molRad = molRad          
  #parms.zLig = zLig
  #parms.zProt = zProt
  #parms.ionC = ionC                
  #parms.A = A
  pmf = CalcPMF(mesh,pmfType=pmfType) 
     
  
  # double check size just in case
  #mn = np.min(mesh.coordinates(),axis=0)
  #mx = np.max(mesh.coordinates(),axis=0)
  #boxVol = np.prod(mx-mn)
  #sphereVol = 4/3.*np.pi*(molRad**3)
  #print meshName
  #print "Vol frac %f vs %f " % (sphereVol/boxVol,volFrac)


  #molDomUnit.smolMode = smolMode
  #solve_homogeneous_unit(molDomUnit,type="field",debug=debug,smolMode=smolMode)
  #homog.SolveHomogSystem(molPrefix=meshPrefix,molGamer=molGamer)  
  #quit()
  #results = homog.SolveHomogSystem(molPrefix=meshPrefix,molGamer=molGamer,\
  #  smolMode=True,smolPsi=pmf,smolq = zLig)
    
  # Note: We use q=1 here since we are expecting potentials in units
  #       of [kT], e.g. already multiplied by zLig 
  results = hl.runHomog(fileXML=meshName,psi=pmf,q=1.,smolMode=True)
  Dx = results.d_eff[0]
  volFrac = results.volFrac
    
  #(V,x) = pb.SolvePoissonBoltzmann(mesh)
  # do electro homog
  # store value
  #results[i,j] = np.random.rand() +  i*j

  # moleculardomain.q =q 
  # moleculardomain.psi = psi
  
  #Dx = hwe.doit(meshFile=meshName, meshType="gamer")
  return Dx,volFrac
  
  
def validation():
  meshPrefix = "example/volfracs/volFrac_0.10" 
  meshPrefix = path+"volFrac_0.10" 
  meshName = meshPrefix+".xml.gz"
  
  zProt = -3
  molRad = 12.5 # [A] 
  parms.A = 0.2 * kT_to_J
  parms.ionC = 0.1 # [M] 

  parms.zLig=1
  pmfType="DebyeHuckel"
  pmfType="DLVO"          
  valuep,phi =  call(meshName,molGamer=0,pmfType=pmfType)             
  parms.zLig=0.
  #value0,phi = results.d_eff[0] =  call(meshPrefix,zLig,molRad,zProt,molGamer=0,pmfType=pmfType)             
  value0,phi =call(meshName,molGamer=0,pmfType=pmfType)             
  parms.zLig=-1
  valuen,phi =  call(meshName,molGamer=0,pmfType=pmfType)             
  #value130404=0.535335
  #value130520=0.55234253        

  print "sphere, neutral ", value0
  print "sphere, positive ", valuep
  print "sphere, negative ", valuen

  assert(np.abs(value0-0.8748) < 0.001), "RESULT CHANGED. DO NOT COMMIT"
  assert(np.abs(valuep-0.9592) < 0.001), "RESULT CHANGED. DO NOT COMMIT"
  assert(np.abs(valuen-0.7989) < 0.001), "RESULT CHANGED. DO NOT COMMIT"


# Semi-validated
# Compared dlvo_valid.png with dlvotest.png for mesh and UnitCube 
#  
def test():
  molRad = 12.5  # [A] assuming thsi is correct for all meshes, and it is only the box size that varies 


  #meshName= "tmp.xml"
  #makeMesh=True
  #if(makeMesh):
  #  molRad = 12.5 # [A]
  #  boxDiam = 65. # [A] 
  #  cvm.Make3DMesh(rSphere=molRad,dBox=boxDiam,name=meshName)
    

  # semi-validated 
  zLig=-1
  #zProt =-4 # -25 mV
  #zProt =-5.4  # -25 mV
  zProt = -3
  ionC = 0.1
  debug = 0

  ## test 
  pb.parms.zProt = zProt

  pmfType = "DLVO" 
  
  ## check that repolsive electro + attract VDW 'beats' purely neutral case 
  if(0): 
    meshName= path+"volFrac_0.10.xml.gz" 
    # neutral 
    parms.zLig = 0.
    parms.A = 0.
    valueN,phi = call(meshName,zLig,molRad,zProt,ionC=ionC,A=A,molGamer=0,debug=debug,pmfType=pmfType)             
    
    # -1 with A   
    A = 3 * kT_to_J # [J] Hamakaer
    zLig=-1  
    valueM,phi = call(meshName,zLig,molRad,zProt,ionC=ionC,A=A,molGamer=0,debug=debug,pmfType=pmfType)             

  if(0): 
    meshName= path+"volFrac_0.50.xml.gz" 
    # neutral 
    zLig = 0.
    A = 0.
    valueN,phi = call(meshName,zLig,molRad,zProt,ionC=ionC,A=A,molGamer=0,debug=debug,pmfType=pmfType)             
    
    # -1 with A   
    A = 0.5 * kT_to_J # [J] Hamakaer
    zLig=-1  
    valueM,phi = call(meshName,zLig,molRad,zProt,ionC=ionC,A=A,molGamer=0,debug=debug,pmfType=pmfType)             

  if(1): 
    meshName= path+"volFrac_0.34.xml.gz" 

    valueMs=[]
    valueNs=[]
    As = np.arange(0.01,1.1,0.1) *kT_to_J
    print As
    for i,A in enumerate(As): 
      zLig=-1  
      valueM,phi = call(meshName,zLig,molRad,zProt,ionC=ionC,A=A,molGamer=0,debug=debug,pmfType=pmfType)             
      valueMs.append(valueM)

      zLig = 0.
      valueN,phi = call(meshName,zLig,molRad,zProt,ionC=ionC,A=A,molGamer=0,debug=debug,pmfType=pmfType)             
      valueNs.append(valueN)


    plt.figure()
    plt.plot(As,np.asarray(valueMs),label="-1")
    plt.plot(As,np.asarray(valueNs),label="-1")
    plt.gcf().savefig("range.png") 
    
    
    

  print valueM,valueN,phi
  quit()
  return valueN

## params 
def runner():            
  meshes = np.array([0.05,0.1,0.2,0.27,0.34,0.5])
  # kappa hard coded into PB solver
  #qs = np.array([-2,-1,0,1,2])
  qs = np.array([-1,0,1])
  #qs = np.array([-1,0])           
  parms.molRad = 12.5  ; # [A] assumed when making original meshes  
  parms.zProt = -3     
  molRad = parms.molRad
  zProt = parms.zProt
  pmfType = "DLVO" 
  
  
  #if(debug==1):
  #  qs=[0]
  #  volFracs=[0.1]
  #  #qs=[2]
  #  volFracs=[0.27]
#    meshes = volFracs
#
  results = np.zeros([ np.shape(meshes)[0], np.shape(qs)[0] ])
  phis = np.zeros(np.shape(meshes)[0])
  for i, mesh in enumerate(meshes):
    for j, q in enumerate(qs):
      #meshName = "meshes/volFrac_%4.2f.xml" % meshes[i]
      #meshName = "meshes/volFrac_%4.2f_mesh.xml.gz" % meshes[i]
      #meshPrefix= "example/volfracs/volFrac_%4.2f" % meshes[i]  
      meshPrefix= path+"/volFrac_%4.2f" % mesh
      meshName=meshPrefix+".xml.gz"
      #molGamer = 1
      molGamer = 0
      parms.zLig = q
      results[i,j],phis[i] = call(meshName,molGamer=molGamer,pmfType=pmfType)
      #print results[i,j],phis[i]
      #quit()


  return results,phis
  
  
def doit(asdf):
  parms.A = 1e-20 # [J] Hamaker   
  results,volFracs=runner()
    
  plt.figure()
  col = ["r--","r-","k-","b-","b--"]
  for j, q in enumerate(qs):
    label = "zLig= %d " % q
    plt.plot(volFracs,results[:,j], col[j],label=label)
  
  title="Protein with DLVO interactions (zProtein=%d)"\
     % zProt
  plt.title(title)
  plt.xlabel("$\phi$")
  plt.ylabel("D")
  plt.legend(loc=0)
  
  plt.gcf().savefig("volfrac.png")
  

# shows difference between Israelach, and more protein-sized obstacles 
def test1(A_kT=5.): 

  validateDLVO()

  # for protein 
  A = A_kT * kT_to_J # 5 [kT] --> J
  R = 12.5e-10 # [nm]
  runDLVO(A=A,R=R,name="smallsphere.png")
  plt.ylim([-2,1])

  A = 0.001 * kT_to_J # 5 [kT] --> J
  runDLVO(A=A,R=R,name="smallsphere_novdw.png")
  plt.ylim([-2,1])




# Generates figure comparing different ligand charges, protein sizes and 
# PMF representations (purely electro vs. DLVO) 
def final(): 
  parms.A = 0.2 * kT_to_J # Guess at Hamakar constant for small ligand with protein 
  resultsDLVO,volFracs=runner()
  parms.A = 0.00001 * kT_to_J # No VDW
  resultsElectroOnly,volFracs=runner()         
    
  #volFracs = 1-np.array([0.1,0.2,0.27,0.34,0.5])
  qs=np.array([-1,0,1])
  plt.figure()
  plt.subplot(211)
  #col = ["r--","r-","k-","b-","b--"]
  col = ["r-","k-","b-"]  
  for j, q in enumerate(qs):
    label = "zLig= %d " % q
    plt.plot(volFracs,resultsElectroOnly[:,j], col[j],label=label)
    
  #col = ["r--","r-","k-","b-","b--"]
  col = ["r-.","k-.","b-."]  
  for j, q in enumerate(qs):
    label = "zLig= %d +VDW" % q
    plt.plot(volFracs,resultsDLVO[:,j], col[j],label=label)    
    
  phi = volFracs
  HS = phi/(2-phi)
  plt.plot(volFracs,HS, "k.",label="HS (cylinder)")

  title="Protein with DLVO interactions (zProtein=%d)"\
     % parms.zProt
    
  plt.title(title)
  plt.xlabel("$\phi$")
  plt.ylabel("D")
  #plt.legend(bbox_to_anchor = (1.5, 0.7),ncol=1)
  plt.legend(bbox_to_anchor = (1.0,-0.2),ncol=3)
  plt.gcf().savefig("final.png") 
    

# Shows the role of chemica specificity in changing diffusion rates
def chemSpecificity():
  # get data for non-specific
  parms.A =0.2*kT_to_J
  results,phis = runner()
  
  # at molFrac 0.5 (12.5*2 = 25 + margin)
  # 27*27*27 cubic angstroms into decimeters cubed 
  # 1.9e-23
  B = 1.9683e-23*6.02e23 # seems high, but that about right for 55 waters
  
  # use a more modest value
  B=1e-6 # [M] 
  
  # store Ds for z=-1, z=1
  Dn=results[:,0]
  Dp=results[:,2]
  
  plt.figure()
  plt.plot(phis,Dn,"r-",label="z=-1")
  exps=np.array([0,-6,-7])
  KDs = 10.**(exps )
  cols = ['b-','b--','b-.']
  for i,KD in enumerate(KDs):
    # from Cheng, PLOS Comp Bio
    Dp_KD =  Dp / (1 + B/KD)
    plt.plot(phis,Dp_KD,cols[i],label="z= 1, $K_D$=10$^{%d}$"%exps[i])
      
  plt.legend(loc=0)    
  plt.ylim([0,2.5])    
  plt.title("Chemical specificity") 
  plt.xlabel("$\phi$") 
  plt.ylabel("D") 
  plt.gcf().savefig("chemspec.png") 





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
  Compute diffusion constant for a variety of geometries 
 
Usage:
"""
  msg+="  %s -debug/-validation" % (scriptName)
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
    if(arg=="-debug"):
      debug=1
    if(arg=="-validation"):
      validateDLVO()
      validation()
      quit()
    if(arg=="-run"): 
      doit(fileIn)
      quit()  
    if(arg=="-test1"): 
      test1()
      quit()  
    if(arg=="-test"): 
      test()
      quit()
    if(arg=="-final"): 
      chemSpecificity()
      final()
      quit()




