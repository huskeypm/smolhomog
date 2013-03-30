import sys
import numpy as np
import matplotlib.pylab as plt
sys.path.append("/home/huskeypm/sources/smolfin/")

from dolfin import *
meshName = "meshes/volFrac_0.10.xml"

# quick test
#mesh = Mesh(meshName)
#import poissonboltzmann as pb
#pb.params.molRad = 12.5
#pb.params.domRad = 200
#(V,x) = pb.SolvePoissonBoltzmann(mesh,meshType="gamer")

import homog_w_electro as hwe
hwe.pb.molRad = 12.5
hwe.pb.domRad = 200
hwe.q=1.
Dx = hwe.doit(meshFile="meshes/volFrac_0.10.xml",meshType="gamer")



volFracs = np.array([0.1,0.2,0.27,0.34,0.5])
meshes = np.array([0.1,0.2,0.27,0.34,0.5])
# kappa hard coded into PB solver
qs = np.array([-1,0,1,2])


#qs=[0]
#volFracs=[0.1]

results = np.zeros([ np.shape(volFracs)[0], np.shape(qs)[0] ])

for i, volFrac in enumerate(volFracs):
  for j, q in enumerate(qs):
    meshName = "meshes/volFrac_%4.2f.xml" % meshes[j]
    mesh = Mesh(meshName)

    # double check size just in case
    mn = np.min(mesh.coordinates(),axis=0)
    mx = np.max(mesh.coordinates(),axis=0)
    boxVol = np.prod(mx-mn)
    sphereVol = 4/3.*np.pi*(hwe.pb.molRad**3)
    print "Vol frac %f vs %f " % (sphereVol/boxVol,volFrac)

    #(V,x) = pb.SolvePoissonBoltzmann(mesh)
    # do electro homog
    # store value
    #results[i,j] = np.random.rand() +  i*j
    hwe.q=1.
    Dx = hwe.doit(meshFile=meshName, meshType="gamer")
    results[i,j] = Dx


plt.figure()
col = ["r-","k-","b-","b--"]
for j, q in enumerate(qs):
  label = "qSubstrate = %d " % q
  plt.plot(volFracs,results[:,j], col[j],label=label)

plt.title("Diffusion hindrance with respect to volume fraction")
plt.xlabel("Volume fraction")
plt.ylabel("Dx/Dbulk")
plt.legend()

plt.gcf().savefig("volfrac.png")









