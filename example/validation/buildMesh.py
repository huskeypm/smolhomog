
boxRad = 1000 # A 
sphereRad = 800 # A 
res = 100 # division size 

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




