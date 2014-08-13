"""
For making simple sphere-centric meshes with prescribed volume
fractions
"""

import sys
#sys.path.append("/home/huskeypm/sources/smolfin/")
from meshmanipulate import *

#
# Revisions
#       10.08.10 inception
#

import numpy as np
def Make3DMesh(volFrac=0.27,rSphere=12.5,dBox=65.,name="test"):
  from gamer import SurfaceMesh, GemMesh
  sm = SurfaceMesh("./meshes/cube_high.off")
  center(sm)

  rBox = dBox/2.
  scaleGeom(sm,rBox)

  inner = SurfaceMesh(4)
  scaleGeom(inner,rSphere)

  inner.as_hole = True
  gem_mesh = GemMesh([sm, inner])
  gem_mesh.write_dolfin(name)

import sys
def MakeGmsh(volFrac=0.27,rSphere=12.5,dBox=65.,name="test"):
  """
  Creates gmsh.geo files for 2D meshes according to vol fraction 
  (only valid for densely packed cylinders) 
  foreach geo, convert to twoD mesh
               dolfin-convert volFrac_0.27.msh  volFrac_0.27.xml
  """
  center = np.zeros(3)
  rBox = dBox/2. 
  cont="Point(1) = {%f,%f,%f};\n" % (center[0],center[1],center[2])
  cont+="Point(2) = {%f, 0, 0};\n" % (rSphere)
  cont+="Point(3) = {%f, 0, 0};\n" % (-rSphere)
  cont+="Point(4) = {0, %f, 0};\n" % (rSphere)
  cont+="Point(5) = {0, %f, 0};\n" % (-rSphere)
  cont+="Point(7) = {%f, %f, 0};\n" % (rBox,rBox)
  cont+="Point(8) = {%f, %f, 0};\n" % (rBox,-rBox)
  cont+="Point(9) = {%f, %f, 0};\n" % (-rBox,rBox)
  cont+="Point(10) = {%f, %f, 0};\n" % (-rBox,-rBox)
  cont+="""
  Circle(5) = {4, 1, 2};
  Circle(6) = {2, 1, 5};
  Circle(7) = {5, 1, 3};
  Circle(8) = {3, 1, 4};
  Line(9) = {7, 8};
  Line(10) = {8, 10};
  Line(11) = {10, 9};
  Line(12) = {9, 7};
  Line Loop(13) = {12, 9, 10, 11};
  Line Loop(14) = {5, 6, 7, 8};
  Plane Surface(15) = {13, 14};
  """
  
  textFile = open(name,"w")             
  textFile.write(cont)
  textFile.close()        



def doit(fileIn):
  rSphere = 12.5

  volFracs = [0.10,.20,.27,.34,.50]

  for i,volFrac in enumerate(volFracs):
    vSphere = (4/3.*np.pi * rSphere**3)
    vBox = vSphere/volFrac
    dBox = vBox**(1/3.)
    #name = "volFrac_%4.2f.geo" % volFrac
    # test with 2D MakeGmsh(volFrac=volFrac,rSphere=rSphere,dBox=dBox)
    name = "./meshes/volFrac_%4.2f.xml" % volFrac
    Make3DMesh(volFrac=volFrac,rSphere=rSphere,dBox=dBox,name=name)
    



if __name__ == "__main__":
  msg="""
Purpose: 
  Creates volumetric meshes for different vol fractions  

Usage:
  script.py run   

Notes:

"""
  remap = "none"


  import sys
  if len(sys.argv) < 2:
      raise RuntimeError(msg)

  fileIn= sys.argv[1]
  if(len(sys.argv)==3):
    print "arg"

  for i,arg in enumerate(sys.argv):
    if(arg=="-arg1"):
      arg1=sys.argv[i+1] 




  doit(fileIn)


