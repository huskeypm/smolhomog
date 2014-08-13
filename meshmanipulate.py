
import numpy as np
class empty: pass

def center(mesh):
  min = empty()
  max = empty()
  min.x=9999
  min.y=9999
  min.z=9999
  max.x=-min.x
  max.y=-min.y
  max.z=-min.z
#  print mesh.vertex(1).x
  
  for i in range(mesh.num_vertices):
      min.x = np.min([min.x, mesh.vertex(i).x])
      max.x = np.max([max.x, mesh.vertex(i).x])
      min.y = np.min([min.y, mesh.vertex(i).y])
      max.y = np.max([max.y, mesh.vertex(i).y])
      min.z = np.min([min.z, mesh.vertex(i).z])
      max.z = np.max([max.z, mesh.vertex(i).z])

  xlat =empty()
  xlat.x=-(max.x+min.x)/2
  xlat.y=-(max.y+min.y)/2
  xlat.z=-(max.z+min.z)/2
  #print min.y
  #print max.y
  #print xlat.z
  #print xlat.y
  #print xlat.z

  translate(mesh,xlat)
  print (min.x,min.y,min.z)
  print (max.x,max.y,max.z)

def scaleGeom(mesh,ratio):
  rat = np.zeros(3)
  if(np.size(ratio)==1):
    rat[:] = ratio
  else:
    rat[:] = ratio[:]

  for i in range(mesh.num_vertices):
      mesh.vertex(i).x = mesh.vertex(i).x*rat[0]
      mesh.vertex(i).y = mesh.vertex(i).y*rat[1]
      mesh.vertex(i).z = mesh.vertex(i).z*rat[2]
#      mesh.vertex(i).x = mesh.vertex(i).x*atom.radius*ratio+atom.x;
#      mesh.vertex(i).y = mesh.vertex(i).y*atom.radius*ratio+atom.y;
#      mesh.vertex(i).z = mesh.vertex(i).z*atom.radius*ratio+atom.z;

def translate(mesh, transvec):
  for i in range(mesh.num_vertices):
      mesh.vertex(i).x = mesh.vertex(i).x+transvec.x
      mesh.vertex(i).y = mesh.vertex(i).y+transvec.y
      mesh.vertex(i).z = mesh.vertex(i).z+transvec.z
#      print "%f %f %f " % (mesh.vertex(i).x,mesh.vertex(i).y,mesh.vertex(i).z)


