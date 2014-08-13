
boxRad = 1000 # A 
sphereRad = 800 # A 
res = 100 # division size 

import numpy as np 


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
  
  print "Point(%d) = {%f,%f,0,1.};" % (n0,n0c[0],n0c[1])
  print "Point(%d) = {%f,%f,0,1.};" % (n1,n1c[0],n1c[1])
  print "Point(%d) = {%f,%f,0,1.};" % (n2,n2c[0],n2c[1])
  print "Point(%d) = {%f,%f,0,1.};" % (n3,n3c[0],n3c[1])
  print "Point(%d) = {%f,%f,0,1.};" % (n4,n4c[0],n4c[1])
  print "Circle(%d) = {%d,%d,%d};" % (n0,n1,n0,n2)            
  print "Circle(%d) = {%d,%d,%d};" % (n1,n2,n0,n3)            
  print "Circle(%d) = {%d,%d,%d};" % (n2,n3,n0,n4)            
  print "Circle(%d) = {%d,%d,%d};" % (n3,n4,n0,n1)            
  ll+=1
  print "Line Loop(%d) = {%d,%d,%d,%d};" % (ll,n0,n1,n2,n3)

  nn = n0+5
  return(nn,ll)

def makeInc(inc): 
  ## print header 
  head="""
Point(1) = {-8, 8, 0, 1.0};
Point(2) = {-8, -8, 0, 1.0};
Point(3) = {8, -8, 0, 1.0};
Point(4) = {8, 8, 0, 1.0};
Line(1) = {2, 1};
Line(2) = {1, 4};
Line(3) = {4, 3};
Line(4) = {3, 2};
Line Loop(1) = {2, 3, 4, 1};
  """
  ll=1
  print head
  
  
  ## Add all circles 
  bigInc= 2
  nstart = np.array([-7,-7])
  nref = nstart + [0,0] 
  n0 = 5
  for i in np.arange(8): 
    for j in np.arange(8): 
      nref = nstart + [i*bigInc,j*bigInc] 
      (n0,ll) = addCircle(n0,nref,inc,ll)
      #print nref
  
  ## add plane tag
  plane="Plane Surface(14) = {"
  for i in np.arange(ll):
    if(i==0):
      pre=""
    else:
      pre=","
    plane+=pre+"%d" % (i+1)
  plane+="};"
  print plane 



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



def makeGmshMesh(boxRad,sphereRad,res): 
  ## build mesh 
  caseName = "cylinderTmp.geo"
  msh = caseName.replace(".geo",".msh") 
  xml = caseName.replace(".geo",".xml") 

  print res
  makeGmshFile(boxRad,sphereRad,res,caseName=caseName)

  from subprocess import call,check_output
  cmd1="gmsh -2 %s"%caseName
  cmd2 ="dolfin-convert %s %s" % (msh,xml) 
  # only way this seems to work 
  call([cmd1],shell=True)
  call([cmd2],shell=True)
  return xml


import sys
#
# Revisions
#       10.08.10 inception
#

def doit(fileIn):
  xmlName = makeGmshMesh(boxRad,sphereRad,res)
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
  





  raise RuntimeError("Arguments not understood")




