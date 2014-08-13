import numpy as np 


import sys
root = "/home/AD/pmke226/sources/"
sys.path.append(root+"/smolhomog/example/validation/")
import buildMesh as bm

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
  msg+="  %s <increment> " % (scriptName)
  msg+="""
  
  where inc = 0.1-0.9
 
Notes:

"""
  remap = "none"

  if len(sys.argv) < 2:
      raise RuntimeError(msg)

  inc = sys.argv[1]
  if(len(sys.argv)==3):
    print "arg"

  bm.makeInc(np.float(inc))


