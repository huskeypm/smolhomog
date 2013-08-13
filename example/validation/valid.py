
import sys
sys.path.append("/net/home/huskeypm/Sources/modified-pb/example")
sys.path.append("/home/huskeypm/sources/modified-pb/example")
sys.path.append("/home/huskeypm/sources/smolhomog/example/noobstacle/")
sys.path.append("/net/home/huskeypm/Sources/smolhomog/example/noobstacle/")
sys.path.append("/home/huskeypm/sources/smolhomog/example/noobstacle/")
sys.path.append("/net/home/huskeypm/Sources/homogenization/")
sys.path.append("/net/home/huskeypm/Sources/smolfin/")           
import matplotlib.gridspec as gridspec
import noobstacle as test
import homoglight as hl

import poissonboltzmann as pb

#
# Still in development 
#
debug = 0
from dolfin import *
import matplotlib.pyplot as plt
import sys
import numpy as np
import poissonboltzmann as pb
parms = pb.parms

parms.res =3       
m_to_A = 1e10
 
parms.z = -1.  # lig charge    [Chloride]
    
#F=96485.3365 # Faradays constant [C/mol]
#R = 8.3143   # Gas const [J/mol K] 
#T = 298.     # Temp [K] 
#RT_o_F=(R*T)/F*1e3  # [mV]   
#print RT_o_F
#Fz_o_RT=z/RT_o_F     # [1/mV] 

import buildMesh   

from scipy.interpolate import griddata

# for interpolating solutions onto 2d mesh 
def interp2d(mesh,x,mode="line"):
      maxr = parms.domRad    
      minr = -parms.domRad    
      if(mode=="line"):   
        res = 10000
        (gx,gy) = np.mgrid[0:0:1j,minr:maxr:res*1j]
        interp = griddata(mesh.coordinates(),x.vector(),(gx,gy))
        interp[np.isnan(interp)]=0
        interp = np.reshape(interp,res)
        gy = np.reshape(gy,res)

        return gy,interp

      if(mode=="plane"):# quarter plane
        res = 1000
        (gx,gy) = np.mgrid[0:maxr:res*1j,0:maxr:res*1j]
        interp = griddata(mesh.coordinates(),x.vector(),(gx,gy))
        interp[np.isnan(interp)]=0
        return gx,gy,interp

def fig7Murad():
  print "WARNING: should do for bilayer, not cylinder [graham eqn, etc]"
  parms.res=0.75 # slow
  parms.res=0.50 # slow
  #parms.res=3 # fast  
  # separate spheres by 2nm [equiv to H=1nm in Fig 7]
    
  domRads=[20.,30.] # A 
  parms.molRad=10. # A   
  abssigma = 0.01 #  \cite{Anonymous:pvWPv9jI} for lipid bilayer 
  nCbs=15
  cbs =10.**(-1*np.linspace(0,6,nCbs))
  nDomRads = np.shape(domRads)[0] 

  #cbs=[1,1]
  #domRads=[1,1]   
  #nDomRads=2
  #nCbs = 2
    
  Dms = np.zeros([nCbs,nDomRads])  
  Dps = np.zeros([nCbs,nDomRads])  
    
  #molRads=[molRads[8]]; 
  #sigmas=[sigmas[0] ]
  for i,cb in enumerate(cbs):
    for j,domRad in enumerate(domRads):
        class empty:pass

        #results = empty()
        #results.Ds=[1,1]
        #Dms[i,j] = results.Ds[0]
        #Dps[i,j] = results.Ds[0]
        #continue

        parms.domRad=domRad # A 
        parms.ionC = cb
        
        parms.sigma = -1.*abssigma
        results = runCase()
        Dms[i,j] = results.Ds[0]
        print results.Ds[0]    

        parms.sigma = abssigma
        results = runCase()
        Dps[i,j] = results.Ds[0]
        print results.Ds[0]    
        
  ps=['b','b--']
  ms=['r','r--']
    
  fig=plt.figure()
  #ax = plt.subplot(211)  
  # make plot with different sizes of subplots 
  gs = gridspec.GridSpec(2, 1, height_ratios=[5, 1]) 
  ax = plt.subplot(gs[0])
  #plt.subplots_adjust(hspace= 0)    

  for j in range(nDomRads):
    ax.plot(cbs,Dms[:,j],ms[j],label="$\sigma=%5.3f$ $[C/m^2]$, H=%3.1f [A]" %\
         (-1*abssigma,domRads[j]))
    ax.plot(cbs,Dps[:,j],ps[j],label="$\sigma=%5.3f$ $[C/m^2]$, H=%3.1f [A]" %\
         ( 1*abssigma,domRads[j]))

  #plt.plot(cbs,Dms,'r',label="$\sigma=%5.3f$ $[C/m^2]$" % (-1*abssigma))
  #plt.plot(cbs,Dps,'b',label="$\sigma=%5.3f$ $[C/m^2]$" % (abssigma))
  ax.set_yscale('log')
  ax.set_xscale('log')
  ax.set_ylabel('D')
  ax.set_xlabel('[cb] [M]')
  ax.legend(loc=0,ncol=2,bbox_to_anchor=(1.10, -0.2))        
  plt.gcf().savefig("fig7muradvalid.png")        
  
        



# for validating PB equation against Fig 3 of Bourbatacha
def fig3():
  parms.res=3
  plt.figure()
  parms.mode='linear'
  # Mesh is made according to Bourbatache
  # R = 10^-8 [M] --> 100 nm --> 1000A
  # separation is 
  # 2e-8 m --> 20 nm, which implies R + 10 nm per box half-width --> 1100A
  fileIn= "/net/home/huskeypm/Sources/smolhomog/example/validation/cylinder.xml"
  fileIn= "/home/huskeypm/sources/smolhomog/example/validation/cylinder.xml"
  parms.molRad=1e-8*m_to_A  # A 
  parms.domRad=parms.molRad+1e-8*m_to_A # A 
    
  fileIn = buildMesh.makeGmshMesh(parms.domRad,parms.molRad,parms.res)            
  mesh = Mesh(fileIn)
  
  # DH   
  #parms.ionC = 0.1 # M 
  #parms.z = -2
  #boundaryPotential="DH"
  #(V,x)= SolvePoissonBoltzmann(mesh,boundaryPotential=boundaryPotential)
  #gy,interp=interp2d(mesh,x)
  #plt.plot(gy,interp,"b",label="DH [I] %3.1f M" % parms.ionC)
  
  # VALIDATED Potential looks validated against Fig 3 of Bourbatche 2012 paper
  sigma = -0.01 # C/m^2
  cb = 50 * 1e-3 # (no lower, since otherwise Grahame invalid) mol/m^3 --> M
  
  parms.ionC = cb
  boundaryPotential = pb.Grahame(sigma,cb)
  (V,x)= pb.SolvePoissonBoltzmann(mesh,boundaryPotential=boundaryPotential)
  gy,interp=interp2d(mesh,x)
  plt.plot(gy,interp,"r",label="$\psi_0$=%4.1f mV [I] %4.2f M" % (boundaryPotential,parms.ionC))

  cb = 200 * 1e-3 # mol/m^3 --> M
  parms.ionC = cb
  boundaryPotential = pb.Grahame(sigma,cb)
  (V,x)= pb.SolvePoissonBoltzmann(mesh,boundaryPotential=boundaryPotential)
  gy,interp=interp2d(mesh,x)
  plt.plot(gy,interp,"g",label="$\psi_0$=%4.1f mV [I] %4.2f M" % (boundaryPotential,parms.ionC))
   
  plt.legend(loc=4)
  plt.xlim([parms.molRad,parms.domRad])
  #plt.ylim([0,-45])
  plt.xlabel("[A]")     
  plt.ylabel("$\psi$ [mV]") 
  #plt.gca().invert_yaxis()
  plt.gcf().savefig("fig3valid.png")


# for validated PB eqn on mesh 
def fig7():   
    #res = 10 # min res for cb=86, sigma = 0.0185
    parms.molRad=30; # 0.3 e-8 ==> A
    parms.domRad=50; # 0.5e-8 m ==> 
    fileIn = buildMesh.makeGmshMesh(parms.domRad,parms.molRad,parms.res) 
 
    cb=86e-3; sigma = -0.0185; # for comparing against Fig 7 ESP 
    boundaryPotential = pb.Grahame(sigma,cb)
    mesh = Mesh(fileIn)
    
    (V,potential)= pb.SolvePoissonBoltzmann(mesh,boundaryPotential=boundaryPotential)
    (gy,interp) = interp2d(mesh,potential)
    plt.figure()
    plt.plot(gy,interp)
    plt.xlim([parms.molRad,parms.domRad])
    plt.figure()
    plt.subplot(111,aspect='equal')
    (gx,gy,interp) = interp2d(mesh,potential,mode="plane")
    plt.pcolormesh(gx,gy,interp)
    plt.xlabel("A")
    clb=plt.colorbar()
    clb.set_label("[mV]")
    plt.gcf().savefig("fig7valid.png")

    

  #Fig. 9. Relative homogenized diffusion coefficient 
  #report Dhom/D for ep=0 -- 1.
  # for sig = 0.07, 0.05,  0.001
  #might also compare w Fig 7 of Murad, but should probably suffice to show opposite
  #tends noted for D- vs D+ [actually, they aren't opposites, so look closer]  

def runCase(engine="testing"):     
    if(debug):
      fileIn = "cylinderTmp.xml"
    else:
      fileIn = buildMesh.makeGmshMesh(parms.domRad,parms.molRad,parms.res)          
    mesh = Mesh(fileIn)
    
    # get electro potential (OVERRIDING) 
    boundaryPotential = pb.Grahame(parms.sigma,parms.ionC)
    print "boundaryPotential ", boundaryPotential
    # potential is mV 
    (V,potential)= pb.SolvePoissonBoltzmann(mesh,boundaryPotential=boundaryPotential)
    
    #(gy,interp) = interp2d(mesh,potential)
    #plt.figure()
    #plt.plot(gy,interp)
    #plt.figure()
    #(gx,gy,interp) = interp2d(mesh,potential,mode="plane")
    #plt.pcolormesh(gx,gy,interp)
    #plt.figure()
    
    #potential = Function(FunctionSpace(mesh,"CG",1))
    #potential.vector()[:] = 0.  
    
    # multiply potential [mV] by F/RT [1/mV] s.t. potential is unitless 
    unitlessPotential = Function(V)
    unitlessPotential.vector()[:] = parms.Fz_o_RT*potential.vector()[:] 
    
    parms.update()
    if(engine=="testing"):
      results = test.doit(mode="hack2",discontinuous=False,dim=2,fileIn=fileIn,potential=unitlessPotential)
    elif(engine=="homog.py"):
      scaledPotential = Function(V)  
      scaledPotential.vector()[:] = parms.Fz_o_RT*potential.vector()[:]  
      scaledPotential.vector()[:] = parms.kT * scaledPotential.vector()[:]  # (z=+,psi=+) --> V
      results = hl.runHomog(fileXML=fileIn,psi=scaledPotential,q=1,smolMode=True)
      results.Ds =results.d_eff 
    
    return results 

# for validating against fig 9 of bourbatache 
def fig9ops():
  # all cases 
  if(debug):
    parms.res = 10 # quick 
    nSigma = 2 
    nMolRads = 3; # 8  
  
  else: 
    parms.res=0.5   
    nSigma = 7 ## DO NOT CHANGE ME SINCE LABELS ARE HARD CODED  
    nMolRads = 8; # 8  
  
  parms.domRad=0.5e-8 * m_to_A # A 
  sigmas = -1*np.array([0,0.01,0.05,0.1, 0.15]) # [C/m^2]
  #sigmas = np.array([-0.07, -0.05,  -0.001 ]) 
  sigmas = np.linspace(-0.15,0.15,nSigma)       
  molRads = np.linspace(parms.domRad*0.01,parms.domRad*0.9,nMolRads)
  cb = 500 * 1e-3 # mol/m^3 --> M  
  parms.ionC = cb

  outs = np.zeros([nSigma,nMolRads])
  vFracs = np.zeros(nMolRads)  
    
  #molRads=[molRads[8]]; 
  #sigmas=[sigmas[0] ]
  for i,molRad in enumerate(molRads): 
        print molRad
        for j,sigma in enumerate(sigmas):
          print sigma
          parms.sigma = sigmas[j]   
          if(1):  
            # make Mesh 
            parms.molRad=molRad  # A 
            results = runCase()
            vFracs[i] = results.phi
            outs[j,i] = results.Ds[0]
            print results.Ds[0]
            
  return vFracs,outs,sigmas
        
def fig9():
    vFracs,outs,sigmas = fig9ops()
    #vFracs,outs,sigmas = test2()
    ns = np.shape(outs)[0]

    gs = gridspec.GridSpec(2, 1, height_ratios=[2, 1]) 
    ax = plt.subplot(gs[0])
    labs = ['b-','b-.','b--','k-','r--','r-.','r-']
    
    for j in range(ns):
      ax.plot(vFracs, outs[j,:],labs[j],label="$\sigma$=%5.3f $[C/m^2]$" % sigmas[j] )
    
    hs = vFracs/(2-vFracs)      
    ax.plot(vFracs, hs,'k.',label="HS")      
    ax.legend(loc=0,ncol=2,bbox_to_anchor=(0.9, -0.13))        
    ax.set_xlabel("$\phi$") 
    ax.set_ylabel("D")
    ax.set_xlim([0.0,1])
    ax.grid(True)  
    #plt.title("My version of fig 9, bourb")
    plt.gcf().savefig("fig9valid.png")     

def fig8():
  parms.res=1
  parms.domRad=0.5e-8 * m_to_A # A 
  cb = 86 * 1e-3 # mol/m^3 --> M  
  #sigmas = np.array([-0.07, -0.05,  -0.001 ]) 
  sigmas =np.linspace(0.1,-0.1,9)
  parms.molRad = 0.4e-8*m_to_A
  parms.ionC = cb

  nSigma = np.shape(sigmas)[0]
  Ds = np.zeros(nSigma)  
  for i, sigma in enumerate(sigmas):
    parms.sigma = sigma    
    results = runCase() 
    Ds[i] = results.Ds[0]

  print Ds
  idx=np.where(sigmas>=0)
  Dps=Ds[idx]
  sps=sigmas[idx]

  idx=np.where(sigmas<=0)
  Dms=Ds[idx]
  sms =np.abs(sigmas[idx])

  plt.figure()  
  plt.plot(sms,Dms,'b',label="$\sigma< 0$")
  plt.plot(sps,Dps,'r',label="$\sigma> 0$")
  plt.xlabel("$|\sigma|$ $[C/m^2]$")
  plt.ylabel("D")
  plt.legend(loc=0)
  plt.grid(True)
  plt.gcf().savefig("fig8valid.png") 
#test1()
#test2()


# basically a unit test for a single case 
def validation(runSeveral=True):
  parms.res=0.5
  parms.domRad=0.5e-8 * m_to_A # A 
  cb = 86 * 1e-3 # mol/m^3 --> M  
  #sigmas = np.array([-0.07, -0.05,  -0.001 ]) 
  parms.molRad = 0.4e-8*m_to_A
  parms.ionC = cb

  sigmas = np.array([ 0.07,0.,-0.07]) # C/m^2
  Ds130813 = np.array([1.3284,0.3221,0.1310])
  for i,sigma in enumerate(sigmas):
    parms.sigma = sigma # C/m^2

    ## using 'runCase'
    results = runCase() 
    Drun1 = results.Ds[0]
    print "runCase:sigma %f z %f D %f " %(sigma,parms.z,Drun1)

    assert(np.abs(results.Ds[0] - Ds130813[i])<0.01), "Validation case FAILED! Do not commit!!"
    print "Success!"

    ## using original Homog code 
    # smol.ElectrostaticsPMF expects to multiply psi*q/kT  
    # so we first multiply potential [mV] by F/RT [1/mV] s.t. potential is unitless 
    # then multiply by -kT and assume that 'q=1' (since already included in Fz/RT term)
    results2 = runCase(engine="homog.py")
    Drun2 = results2.Ds[0]
    print "homog: sigma %f z %f D %f " %(sigma,parms.z,Drun2)

    assert(np.abs(Drun1-Drun2)<0.001)


def test1():
  parms.res=0.2
  parms.domRad=21.7 # A            
  parms.molRad = 12.0 # A           

  #parms.res=2.0 
  #parms.domRad=217 # A            
  #parms.molRad = 120 # A           

  cb = 0.15 # [M]                       
  #sigmas = np.array([-0.07, -0.05,  -0.001 ]) 
  sigmas =np.linspace(0.1,-0.1,3)
  parms.ionC = cb

  nSigma = np.shape(sigmas)[0]
  Ds = np.zeros(nSigma)
  for i, sigma in enumerate(sigmas):
    parms.sigma = sigma
    results = runCase(engine="homog.py")
    Ds[i] = results.Ds[0]

  plt.figure()
  plt.plot(sigmas,Ds,'b',label="z=%3.1f"%parms.z)           
  plt.xlabel("$\sigma$ $[C/m^2]$")
  plt.ylabel("D")
  plt.legend(loc=0)
  plt.grid(True)
  plt.gcf().savefig("test1.png")



  

  





      

#test1()
#!/usr/bin/env python
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
      validation()

    if(arg=="-fig3"):
      fig3()
    if(arg=="-fig7"):
      fig7()
    if(arg=="-fig7Murad"):
      fig7Murad()
    if(arg=="-fig8"):
      fig8()
    if(arg=="-fig9"):
      fig9()
    if(arg=="-debug"): 
      debug=1
    if(arg=="-test1"):
      test1()


