import numpy as np
# Cretates 'currogated' mesh for inestigating slower diffusion 
  
nm_to_Ang = 10.
def doit(outName):
  nRidges = 20
  periodRidge = 0.4 * nm_to_Ang
  heightRidge = 0.75 * nm_to_Ang
  layerHeight = 2. * nm_to_Ang
  res = layerHeight/30.
  
  ## define  'ridge' positions on base layer 
  ridgeXCoords = np.arange(2*nRidges)/(2/periodRidge)
  ridgeYCoords = np.ones(np.shape(ridgeXCoords)[0]) * heightRidge
  
  ## define base layer 
  # add 'half width' on each side 
  baseXCoords = ridgeXCoords
  baseXCoords = np.concatenate((baseXCoords, [baseXCoords[-1]+periodRidge/4.]))
  baseXCoords = np.concatenate(([baseXCoords[0]-periodRidge/4.],baseXCoords ))
  baseYCoords = np.zeros(np.shape(baseXCoords)[0]) 
  
  ## define top ridges 
  topRidgeXCoords = np.copy(ridgeXCoords)
  topRidgeYCoords = ridgeYCoords + layerHeight- 2 * heightRidge
  
  ## define top layer 
  topXCoords = np.copy(baseXCoords)
  topYCoords = baseYCoords + layerHeight
  
  ## combine into single list of points 
  #np.concatenate((np.arange(5).T,np.arange(5).T),axis=0)
  #np.vstack((np.arange(5),np.arange(5))).T
  allXCoords = np.concatenate((baseXCoords,ridgeXCoords,topRidgeXCoords,topXCoords))
  allYCoords = np.concatenate((baseYCoords,ridgeYCoords,topRidgeYCoords,topYCoords))
  coords = np.vstack((allXCoords,allYCoords)).T
  np.shape(coords)
  
  
  
  
  ## add points 
  nPoints = np.shape(coords)[0]
  lines=""
  for i in np.arange(nPoints):
      line = "Point(%d) = {%f, %f, 0, %f};\n" % (i+1,coords[i,0],coords[i,1],res)
      lines += line
  
  ## add lines (empirical so far)
  lastBasePoint = 2*nRidges + 2
  gap = lastBasePoint -1 
  basePoints  = 2+np.arange(2*nRidges)
  baseRidgePoints  = basePoints + gap
  topRidgePoints = basePoints + 2*gap-1
  topPoints  = topRidgePoints + gap
  connect = np.concatenate((basePoints,topRidgePoints))
  # vertical 
  ctr=0
  for i,c in enumerate(connect):
      ctr+=1
      line = "Line(%d) = {%d,%d};\n" % (ctr,c, c+gap)
      lines += line
      
  # horizontal ridges
  baseTmp = np.arange(nRidges)*2+baseRidgePoints[0]
  topTmp = np.arange(nRidges)*2+topRidgePoints[0]
  for i  in range(nRidges):
      ctr+=1
      line = "Line(%d) = {%d,%d};\n" % (ctr,baseTmp[i], baseTmp[i]+1)
      lines += line
      ctr+=1
      line = "Line(%d) = {%d,%d};\n" % (ctr,topTmp[i], topTmp[i]+1)
      lines += line
  
  # horizontal top/bottom
  baseLeft = 1
  topLeft = topPoints[0] -1
  baseRight = 2*nRidges+baseLeft + 1
  topRight = 2*nRidges+topLeft+1
  baseTmp = np.arange(nRidges+1)*2+baseLeft
  topTmp = np.arange(nRidges+1)*2+topLeft
  for i  in range(nRidges+1):
      ctr+=1
      line = "Line(%d) = {%d,%d};\n" % (ctr,baseTmp[i], baseTmp[i]+1)
      lines += line
      ctr+=1
      line = "Line(%d) = {%d,%d};\n" % (ctr,topTmp[i], topTmp[i]+1)
      lines += line
  
  # corners   
  ctr+=1
  line = "Line(%d) = {%d,%d};\n" % (ctr,baseLeft,topLeft)
  lines += line
  ctr+=1
  line = "Line(%d) = {%d,%d};\n" % (ctr,baseRight,topRight)
  lines += line
  
      
  f = open(outName,'w+')
  f.write(lines)
  f.close()    
  
  print "WARNING: for now the defiition of the outer surface/must b edone manually in gmsh, since i think i need to reorder the lines"

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

 # for i,arg in enumerate(sys.argv):
    #if(arg=="-validation"):
    #  arg1=sys.argv[i+1] 
    #  doit(fileIn)
 # 
  #outName = "/home/huskeypm/test.geo"
  outName = "mesh.geo"
  doit(outName)
  quit()





  raise RuntimeError("Arguments not understood")




