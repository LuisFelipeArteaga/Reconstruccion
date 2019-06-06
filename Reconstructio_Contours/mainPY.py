import pygmsh
import numpy as np
import meshio
from scipy.spatial import Delaunay

def write_ply(ply_file,vertex,faces):
  file = open(ply_file,"w")
  file.write('''ply
format ascii 1.0
element vertex %d
property float x
property float y
property float z
element face %d
property list uchar int vertex_index
end_header
'''%(len(vertex),len(faces)))

  for row in vertex:
    file.write('%f %f %f \n' %(row[0],row[1],row[2]))
  for row in faces:
    file.write('3 %i %i %i \n' %(row[0],row[1],row[2]))
  file.close()


if __name__ == "__main__":
    
    point_cloud = np.loadtxt('Files/point_cloud.txt')#.tolist()
    geom = pygmsh.built_in.Geometry()

    a = [[ -100, -100, -100],[-100, -100,  100], [-100,  100, -100],[100,  100,  100],[  100, -100, -100],[ 100, -100,  100],[ 100, -100,  100],[ 100,  100,  100]]
    tri = Delaunay(point_cloud)
    print tri.simplices
    write_ply('ZZZ.ply',point_cloud,tri.simplices)


# Draw a cross.
    poly = geom.add_polygon(point_cloud.tolist())
    

    axis = [0, 0, 1]

    geom.extrude(
        poly,
        translation_axis=axis,
        rotation_axis=axis,
        point_on_axis=[0, 0, 0],
        angle=0
        )

    points, cells, point_data, cell_data, field_data = pygmsh.generate_mesh(geom)
    meshio.write_points_cells('point_cloud.vtu', points, cells)
