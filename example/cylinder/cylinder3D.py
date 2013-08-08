import sys 
import matplotlib.pylab as plt
#sys.path.append("/net/home/huskeypm/Sources/homogenization/")
import homoglight as hl
import numpy as np

from dolfin import *



#print "WARNING: I 'flipped' the obstacle location here" 
class Omega1(SubDomain):
    def inside(self, x, on_boundary):
      return True

# Create mesh and define function space
def testCase(filename="none",mode="cylinder",qSubstrate=1.,dim=3,qenz=-3):

  # Create mesh and define function space
  mesh = Mesh(filename)
  V = FunctionSpace(mesh, "Lagrange", 1)

  # mostly done just to make sure all the cell elements are marked 
  if(mode=="cylinder"):  
    domains = MeshFunction('uint', mesh, dim)   
    markerFree=1;
    for cell_no in range(len(domains.array())):
      domains.array()[cell_no]=markerFree
    # save 
    fileXML = "cylinder.xml.gz"  
    File(fileXML) << mesh  
  else:
    fileXML = filename
  #File("subdomain.xml.gz") << subdomains

  # potential # KIND OF DEBYEHUCKEL LIKE< BUT NOT CORRCT 
  #qenz = -1. 
  #  
  #qenz = qSubstrate
  #qSubstrate = -1  
    
  if(mode=="cylinder"):  
    expr = Expression("qenz*prefac*exp(-1*sqrt(x[0]*x[0]+x[1]*x[1])/w)",prefac=3,qenz=qenz,w=10)
  elif(mode=="sphere"):
    expr = Expression("qenz*prefac*exp(-1*sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])/w)",prefac=3,qenz=qenz,w=10)  
  psi=Function(V)  
  psi.interpolate(expr)
    
  import numpy as np
  ar =   np.asarray(psi.vector()[:])
  print "Fake psi min/max %f/%f [mV]" % (np.min(ar),np.max(ar))
  psi.vector()[:] *= (1/25.) *0.593 # 1/RT [mV] * kT 
  ar =   np.asarray(psi.vector()[:])
  print "Fake pmf min/max %f/%f [mV]" % (np.min(ar),np.max(ar))
  File("psi.pvd") << psi

    
  ## run homogenization   
  results = hl.runHomog(fileXML=fileXML,psi=psi,q=qSubstrate,smolMode=True)
  

  return results


   
def doit():
  qs = np.array([-1.,0.,1.])
  nqs = 3      
        
  fileIn="cylinder2D.xml"
  Dcyl2ds=np.zeros(nqs)
  for i,q in enumerate(qs):      
    resultsCylinder = testCase(fileIn,mode="cylinder",qSubstrate=q,dim=2)
    Dcyl2ds[i] = resultsCylinder.d_eff[0]
    
  fileIn="cylinder3D.xml"
  Dcyls=np.zeros(nqs)
  for i,q in enumerate(qs):      
    resultsCylinder = testCase(fileIn,mode="cylinder",qSubstrate=q)
    Dcyls[i] = resultsCylinder.d_eff[0]

  Dsphs=np.zeros(nqs)
  fileIn="sphere.xml.gz"
  for i,q in enumerate(qs):      
    resultsSphere = testCase(fileIn,mode="sphere",qSubstrate=q)
    Dsphs[i] = resultsSphere.d_eff[0]
    
  plt.plot(qs,Dcyl2ds,'k.',label="cylinder/2D")
  plt.plot(qs,Dcyls,'b-.',label="cylinder/3D")  
  plt.plot(qs,Dsphs,'r-.',label="sphere/3D")  
  plt.title("Sphere vs cylinder") 
  plt.ylabel("D")
  plt.xlabel("q") 
  plt.xlim([-1.1,1.1])  
  plt.legend(loc=0)

  plt.gcf().savefig("sphere_vs_cylinder.png")  

def test():
  fileIn="/net/home/huskeypm/Sources/homogenization/example/volfracs/volFrac_0.10_mesh.xml.gz"
  fileIn="/home/huskeypm/sources/homogenization/example/volfracs/volFrac_0.10_mesh.xml.gz"
  q=1
  #q=0
  q=-1
  resultsSphere = testCase(fileIn,mode="sphere",qSubstrate=q)
    

#sphere
if __name__ == "__main__":
  msg="Purpose: To compare sphere vs cylinder w fake pmf "
  msg=msg+"Usage: "
  msg=msg+".py <file.xml/-validate>"
  msg=msg+"Notes:"
  remap = "none"



  mode = "none"
  import sys
  if len(sys.argv) < 2:
      raise RuntimeError(msg)

  fileIn= sys.argv[1]
  if(len(sys.argv)==3):
    print "arg"

  dolf =0
  for i,arg in enumerate(sys.argv):
    if(arg=="-gamer"):
      dolf=0

    if(arg=="-run"):
      doit()
 
    if(arg=="-test"):
     test()

