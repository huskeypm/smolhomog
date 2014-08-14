
boxRad = 1000 # A 
sphereRad = 800 # A 
res = 100 # division size 

import numpy as np 

# creates regular 8x8 lattice with 64 points therein 
def makeIncrLattice(inc):
  bigInc= 2
  nstart = np.array([-7,-7])
#  nref = nstart + [0,0]
#  n0 = 5
#  circles =""
  locs = []
  for i in np.arange(8):
    for j in np.arange(8):
      nref = nstart + [i*bigInc,j*bigInc]
      locs.append(nref)

  rads = np.ones( len(locs) ) * inc
  locs = np.asarray(locs,dtype=float) #np.reshape( (len(locs), 2 ) )  
  makeLattice(locs,rads,boxMin=[-8,-8], boxMax=[8,8])


def addCircle(n0,nref,inc,ll):
  n0c = nref
  n1=n0+1
  n1c = n0c+np.array([-1*inc,0]) 
  n2=n0+2
  n2c = n0c+np.array([0,-1*inc]) 
  n3=n0+3
  n3c = n0c+np.array([1*inc,0]) 
  n4=n0+4
  n4c = n0c+np.array([0,1*inc]) 
  
  line="Point(%d) = {%f,%f,0,1.};\n" % (n0,n0c[0],n0c[1])
  line+="Point(%d) = {%f,%f,0,1.};\n" % (n1,n1c[0],n1c[1])
  line+="Point(%d) = {%f,%f,0,1.};\n" % (n2,n2c[0],n2c[1])
  line+="Point(%d) = {%f,%f,0,1.};\n" % (n3,n3c[0],n3c[1])
  line+="Point(%d) = {%f,%f,0,1.};\n" % (n4,n4c[0],n4c[1])
  line+="Circle(%d) = {%d,%d,%d};\n" % (n0,n1,n0,n2)            
  line+="Circle(%d) = {%d,%d,%d};\n" % (n1,n2,n0,n3)            
  line+="Circle(%d) = {%d,%d,%d};\n" % (n2,n3,n0,n4)            
  line+="Circle(%d) = {%d,%d,%d};\n" % (n3,n4,n0,n1)            
  ll+=1
  line+="Line Loop(%d) = {%d,%d,%d,%d};\n" % (ll,n0,n1,n2,n3)

  nn = n0+5
  return(line,nn,ll)


# make lattice given list of location and radii 
def makeLattice(locs,rads,fileName="out.geo",boxMin=[-8,-8], boxMax=[8,8],writeMesh=True): 
  ## print header 
  head="Point(1) = {%f, %f, 0, 1.0};\n" % (boxMin[0],boxMax[1])
  head+="Point(2) = {%f, %f, 0, 1.0};\n" % (boxMin[0],boxMin[1])
  head+="Point(3) = {%f, %f, 0, 1.0};\n" % (boxMax[0],boxMin[1])
  head+="Point(4) = {%f, %f, 0, 1.0};\n" % (boxMax[0],boxMax[1])
  head+="""
Line(1) = {2, 1};
Line(2) = {1, 4};
Line(3) = {4, 3};
Line(4) = {3, 2};
Line Loop(1) = {2, 3, 4, 1};
  """ 
  n0 = 5
  
  
  ## Add all circles 
  ll = 1 
  circles =""
  nCircles = np.shape(locs)[0]
  for i in np.arange(nCircles):
      nref = locs[i,]
      rad = rads[i]
      (line,n0,ll) = addCircle(n0,nref,rad,ll)
      circles+=line
  
  ## add plane tag
  plane="Plane Surface(14) = {"
  for i in np.arange(ll):
    if(i==0):
      pre=""
    else:
      pre=","
    plane+=pre+"%d" % (i+1)
  plane+="};"
  #print plane 

  fileLine = head + circles + plane

  with open(fileName, "w") as text_file:
      text_file.write(fileLine)
  print "Wrote " , fileName

  if writeMesh:
    makeDolfinMeshFromGeo(fileName)



## make box via gmsh and write 
def makeGmshFile(boxRad,sphereRad,res,caseName="cylinderTmp.geo"): 
  fileLine="""
b=%f;
c=%f;
l=%f;
Point(1) = {-b, -b, 0, l};
Point(2) = {b, -b, 0, l};
Point(3) = {b, b, 0, l};
Point(4) = {-b, b, 0, l};
Point(5) = {0.0, 0.0, 0, l};
Point(6) = {-c, 0.0, 0, l};
Point(7) = {c, 0.0, 0, l};
Point(8) = {0.0, -c, 0, l};
Point(9) = {0.0, c, 0, l};

Line(1) = {4, 1};
Line(3) = {1, 2};
Line(4) = {2, 3};
Line(5) = {3, 4};
Circle(6) = {6, 5, 8};
Circle(7) = {8, 5, 7};
Circle(8) = {7, 5, 9};
Circle(9) = {9, 5, 6};
Line Loop(10) = {5, 1, 3, 4};
Line Loop(11) = {8, 9, 6, 7};
Plane Surface(12) = {10, 11};
""" % (boxRad,sphereRad,res)
#print fileLine
  with open(caseName, "w") as text_file:
      text_file.write(fileLine)



def makeCylMesh(boxRad,sphereRad,res): 
  ## build mesh 
  caseName = "cylinderTmp.geo"

  print res
  makeGmshFile(boxRad,sphereRad,res,caseName=caseName)
  makeDolfinMeshFromGeo(caseName)

def makeDolfinMeshFromGeo(caseName):
  msh = caseName.replace(".geo",".msh") 
  xml = caseName.replace(".geo",".xml") 
  from subprocess import call,check_output
  cmd1="gmsh -2 %s"%caseName
  cmd2 ="dolfin-convert %s %s" % (msh,xml) 
  # only way this seems to work 
  call([cmd1],shell=True)
  call([cmd2],shell=True)
  print "Wrote " , xml         
  return xml


import sys
#
# Revisions
#       10.08.10 inception
#

def doit(fileIn):
  xmlName = makeCylMesh(boxRad,sphereRad,res)
  1

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
      arg1=sys.argv[i+1] 
      doit(fileIn)
      quit() 
  





  raise RuntimeError("Arguments not understood")




