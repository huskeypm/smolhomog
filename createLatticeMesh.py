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

def doit(inc): 
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
  For creating meshes with many obstacles (2d)  
 
Usage:
"""
  msg+="  %s <inc> " % (scriptName)
  msg+="""
  
 
Notes:

"""
  remap = "none"

  if len(sys.argv) < 2:
      raise RuntimeError(msg)

  inc = sys.argv[1]
  if(len(sys.argv)==3):
    print "arg"

  doit(np.float(inc))


