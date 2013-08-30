## worked example 14.5, Israelachvilie Surface Forces book
##
## PARAMS 
##
# W(D) = 0.5 R Z exp(-kappa D) - A*R/12D
print "WARNING: NEED TO DOUBLE CHECK THAT THIS IS THE CORRECT VDW TERM"
print "WARNING: DLVO stuff is for LIKE charged surface, not opposite"
print "So what i have is a hack"

import sys
prefix = "/net/home/huskeypm/Sources/"
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
#mine
ionC = 0.1 # [M] 
iKappa =0.304/np.sqrt(ionC); kappa = 1/iKappa

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
  # if validating against Israelach, use 9.38e-11 and psi0/107 for T=38
  absz = np.abs(z)
  if(z==0): 
    signz=0.
  else:
    signz=np.sign(psi0/z)  # THIS IS NOT CORRECT!!!! (sign=1 like charges, sign=-1 opposite)
  Z = signz*9.22e-11 *np.tanh(absz*psi0/103)**2   # [J/m] at T = 25 C 
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
    
  runDLVO(A=A,R=R,minDLVO=-1.5e-19,name="dlvo_valid.png")
    
def runDLVO(A=1e-20,R=0.1e-6,minDLVO=-2*kT_to_J,name="dlvo.png"):    
  Dnm = np.linspace(0.1,8,80) # [nm]
  ## get sigmas from potentials so we can
  # adjust poetntial based on ionic str.
  israelsPsi0s =np.array([0,14.5,24.5,34.5])
  # Graham equation
  sigmas= 0.117*np.sqrt(ionC)*np.sinh(israelsPsi0s/51.4)

  fig = plt.figure()
  ax = fig.add_subplot(111)

  units = "kT"
  cols=['r','b']
  styles=['-','--','-.','.']
  for j,z in enumerate([-1,1]):
    for i,sigma in enumerate(sigmas):
      psi0=pb.Grahame(sigma,ionC)
        
      if(units=="J"):   
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


# psi0 [mV]
# R - particle radius [m]
# A - Hamakaer [J] 
def ApplyDLVO(case="unitsphere",mesh="none",psi0 = 24.5,  z = 1.,R=1e-7,A=1e-20):
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
    R_Ang = R*m_to_Ang
    dexpr = Expression("sqrt(pow(x[0]-cx,2)+pow(x[1]-cy,2)+pow(x[2]-cz,2))-R",\
                        cx=sphereCenter[0],cy=sphereCenter[1],cz=sphereCenter[2],R= R_Ang)          
    #boxDiam = np.max(mesh.coordinates())-np.min(mesh.coordinates())
    #left = np.max(mesh.coordinates())
    #print "Dexpr ", dexpr(left,0,0)
    #print "Dexpr ", dexpr(left,left,left)
    #print "Dexpr ", dexpr(R_Ang,0,0)
    #print "BoxRad %f R %f " % (boxDiam/2., R_Ang)
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
    
  if(1):  
    d = Function(V)
    d.vector()[:] = ds*nm_to_Ang
    File("distances.pvd") << d
        
  #
  print "WARNING: this is a bit of a debugging hack to prevent divide by zeros"
  tol = 1e-19
  ds[ np.where(ds < tol) ] = tol
  plt.plot(ds,np.zeros(np.shape(ds)[0]),"k.")
  plt.xlabel("[nm]")
  # evaluate DLVO expression
  #ws = vW(ds.vector()[:],psi0,z)
  print "z=",z
  print "psi0",psi0  
  #z *= -1
  #psi0 *= -1
 
  ws = vW(ds,psi0,z,A,R) * J_to_kT
  pmf = Function(V)
  pmf.vector()[:] = ws
  plt.plot(ds,ws,"r.")
  plt.gcf().savefig("distances.png") 


  ds*= nm_to_Ang
  return (ds,pmf)

from dolfin import *


# molRad - radius of enzeyme 
# z - enzyme charge
# q - substrate charge 
# A - Hamaker const [J] 
# R - molRad [A] 
def CalcPMF(mesh,R,zLig,zProt,ionC=0.15, meshType="dolfin",pmfType="DebyeHuckel",case="sphere",A=1e-20):
   pb.parms.zLig = zLig
   pb.parms.zProt = zProt
   pb.parms.molRad = R               
   pb.parms.ionC = ionC                
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
     (ds,dummypmf)= ApplyDLVO(case=case,mesh=mesh,R=R*Ang_to_m) 
   
   elif(pmfType=="DLVO"): 
     #f = pb.DLVOExpr(dim=3)
     #f.interpolate(pmf.vector()[:])
     #pmf = interpolate(f,V)       
     expr = pb.DebyeHuckelExpr(dim=3)
     # get potential at 'left' edge of sphere
     psi0 = expr(R,0.,0.)
     print psi0
     (ds,pmf)= ApplyDLVO(case=case,mesh=mesh,psi0=psi0,z=zLig,R=R*Ang_to_m,A=A)


   # compare
   if(1):
     print "compare", np.min(ds), np.max(ds)
     # this is a bit of a hack to get the distances from the ApplyDLVO function 
     # plot
     plt.figure()
     #plt.plot(ds.vector().array(),pmf.vector().array(),'k.')
     plt.plot(ds*Ang_to_nm,pmf.vector().array(),'k.')
     plt.title("%s interaction energy \n(R=%4.2f [nm], H=%4.2f [kT],psi=%4.2f [mV],z=%d)" % (pmfType,R*Ang_to_nm,A*J_to_kT,psi0,zLig))
     plt.xlabel("Dist from surface[nm]") 
     #plt.xlim([0,8])
     plt.ylim([-7,1])     
     plt.grid(True)
     plt.ylabel("Energy [kT]") 
     plt.gcf().savefig("pmfWRTDist.png") 
     File("testpmf.pvd") << pmf

   return pmf
            
# A Hamaker [J] 
def call(meshName,zLig,molRad,zProt,ionC=0.15,A=1e-20,molGamer=0,debug=0,pmfType="DLVO"):
#  hwe.params.z = -1
#  hwe.pb.molRad = 12.5
#  hwe.pb.sphericalDomainBoundary=False
#  hwe.params.molRad = 12.5
#  hwe.params.sphericalDomainBoundary=False
#  hwe.params.q=1.


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
  print pb.parms.zProt 
  pmf = CalcPMF(mesh,molRad,zLig,zProt,A=A,ionC=ionC,pmfType=pmfType) 
     
  
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
  molRad = 12.5 # A !! need to determine on the fly!!

  zLig=1
  pmfType="DebyeHuckel"
  pmfType="DLVO"          
  valuep,phi =  call(meshPrefix,zLig,molRad,zProt,molGamer=0,pmfType=pmfType)             
  zLig=0.
  value0,phi = results.d_eff[0] =  call(meshPrefix,zLig,molRad,zProt,molGamer=0,pmfType=pmfType)             
  zLig=-1
  valuen,phi =  call(meshPrefix,zLig,molRad,zProt,molGamer=0,pmfType=pmfType)             
  #value130404=0.535335
  value130520=0.55234253        

  print "sphere, neutral ", value0
  print "sphere, positive ", valuep
  print "sphere, negative ", valuen

  assert(np.abs(valuep-value130520) < 0.001), "RESULT CHANGED. DO NOT COMMIT"


# Semi-validated
# Compared dlvo_valid.png with dlvotest.png for mesh and UnitCube 
#  
def test():
  molRad = 12.5  # assuming thsi is correct for all meshes, and it is only the box size that varies 


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
    zLig = 0.
    A = 0.
    valueN,phi = call(meshName,zLig,molRad,zProt,ionC=ionC,A=A,molGamer=0,debug=debug,pmfType=pmfType)             
    
    # -1 with A   
    A = 3 * kT_to_J # [J] Hamakaer
    zLig=-1  
    valueM,phi = call(meshName,zLig,molRad,zProt,ionC=ionC,A=A,molGamer=0,debug=debug,pmfType=pmfType)             

  if(1): 
    meshName= path+"volFrac_0.50.xml.gz" 
    # neutral 
    zLig = 0.
    A = 0.
    valueN,phi = call(meshName,zLig,molRad,zProt,ionC=ionC,A=A,molGamer=0,debug=debug,pmfType=pmfType)             
    
    # -1 with A   
    A = 0.5 * kT_to_J # [J] Hamakaer
    zLig=-1  
    valueM,phi = call(meshName,zLig,molRad,zProt,ionC=ionC,A=A,molGamer=0,debug=debug,pmfType=pmfType)             
    
    

  print valueM,valueN,phi
  return valueN

## params 
def runner(A=1e-20):
  meshes = np.array([0.1,0.2,0.27,0.34,0.5])
  #meshes = np.array([0.1,0.2,0.27,0.34,0.5])
  # kappa hard coded into PB solver
  #qs = np.array([-2,-1,0,1,2])
  qs = np.array([-1,0,1])
  #qs = np.array([-1,0])           
  molRad = 12.5  ; print "WARNING: put this in common place" 
  zProt = -3
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
      print "WARNING: molRads are NOT correct" 
      #molGamer = 1
      molGamer = 0
      results[i,j],phis[i] = call(meshName,q,molRad,zProt,A=A,molGamer=molGamer,pmfType=pmfType)


  return results,phis
  
  
def doit(asdf):
  A = 1e-20 # [J] Hamaker   
    
    
  results,volFracs=runner(A=A)
    
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
def TestComparisons(A_kT=5.): 

  validateDLVO()

  # for protein 
  A = A_kT * kT_to_J # 5 [kT] --> J
  R = 12.5e-10 # [nm]
  runDLVO(A=A,R=R,name="smallsphere.png")
  plt.ylim([-2,1])

  A = 0.001 * kT_to_J # 5 [kT] --> J
  runDLVO(A=A,R=R,name="smallsphere_novdw.png")
  plt.ylim([-2,1])




def final(): 
  A = 1e-20 # [J] Hamaker   
  resultsDLVO,volFracs=runner(A=A)
  resultsElectroOnly,volFracs=runner(A=0.)
    
  print "WARNING: vol fracs are wrong, need to check on comuations"
  #volFracs = 1-np.array([0.1,0.2,0.27,0.34,0.5])
  qs=np.array([-1,0,1])
  plt.figure()
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
     % -1000 
    
  plt.title(title)
  plt.xlabel("$\phi$")
  plt.ylabel("D")
  plt.legend(bbox_to_anchor = (1.5, 0.7),ncol=1)
  plt.gcf().savefig("final.png") 
    




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
      validation()
      quit()
    if(arg=="-run"): 
      doit(fileIn)
      quit()  
    if(arg=="-test1"): 
      TestComparisons()
      quit()  
    if(arg=="-test"): 
      test()
      quit()
    if(arg=="-final"): 
      final()
      quit()




