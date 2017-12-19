#include <pcl/point_types.h>
#include <pcl/io/pcd_io.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/features/normal_3d.h>
#include <pcl/surface/gp3.h>
#include <pcl/io/vtk_io.h>
#include <pcl/io/ply_io.h>
#include <pcl/io/obj_io.h>
//#include <pcl/io/stl_io.h>
#include <pcl/PolygonMesh.h>
#include <pcl/TextureMesh.h>
#include <pcl/visualization/cloud_viewer.h>

#include <stdlib.h>  
#include <stdio.h>  

int
main (int argc, char** argv)
{
  std::string ply_file = argv[1];

  // Carga la nube de puntos en formato punto PLY
  pcl::PointCloud<pcl::PointXYZ>::Ptr cloud (new pcl::PointCloud<pcl::PointXYZ>);
  pcl::PCLPointCloud2 cloud_blob;
  pcl::io::loadPLYFile (ply_file, cloud_blob);
  pcl::fromPCLPointCloud2 (cloud_blob, *cloud);
 
  // Estimacion de Normales
  pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> n;
  pcl::PointCloud<pcl::Normal>::Ptr normals (new pcl::PointCloud<pcl::Normal>);
  pcl::search::KdTree<pcl::PointXYZ>::Ptr tree (new pcl::search::KdTree<pcl::PointXYZ>);
  tree->setInputCloud (cloud);
  n.setInputCloud (cloud);
  n.setSearchMethod (tree);
  n.setKSearch (20);
  n.compute (*normals);

  //Como las coordenadas y las normales deben estar en el mismo PointCloud, creamos una nube de puntos de tipo PointNormal.
  pcl::PointCloud<pcl::PointNormal>::Ptr cloud_with_normals (new pcl::PointCloud<pcl::PointNormal>);
  pcl::concatenateFields (*cloud, *normals, *cloud_with_normals);
  
  //crear arbol binario de busqueda
  pcl::search::KdTree<pcl::PointNormal>::Ptr tree2 (new pcl::search::KdTree<pcl::PointNormal>);
  tree2->setInputCloud (cloud_with_normals);

  // Inicializar objetos
  pcl::GreedyProjectionTriangulation<pcl::PointNormal> gp3;
  pcl::PolygonMesh triangles;

  
  gp3.setSearchRadius(atof(argv[2])); //  Longitud de arista máxima para cada triángulo. 
  //Parametros que controlan el tamaño de la vecindad
  gp3.setMu(atof(argv[3])); // Especifica la distancia máxima aceptable para un punto a considerar vecino

  
  gp3.setMaximumNearestNeighbors(atoi(argv[4])); // Define en cuántos vecinos se buscan.

  /************Valores tipicos de Parametros************/
  // Hace frente a los casos en que hay bordes o esquinas afilados y donde dos 
  // lados de una superficie se ejecutan muy cerca uno del otro. 
 
  //son el mínimo y ángulos máximos en cada triángulo.
 
  gp3.setMaximumSurfaceAngle(M_PI/4); // 45 degrees
  gp3.setMinimumAngle(M_PI/18); // 10 degrees
  gp3.setMaximumAngle(2*M_PI/3); // 120 degrees
  gp3.setNormalConsistency(false);


  // Resultado
  gp3.setInputCloud (cloud_with_normals);
  gp3.setSearchMethod (tree2);
  gp3.reconstruct (triangles);

  // Informacion adicional de los vertices
  std::vector<int> parts = gp3.getPartIDs();
  std::vector<int> states = gp3.getPointStates();

  // Almacenar superficie
  pcl::io::saveOBJFile("surface/surface3.obj", triangles);
  return (0);
}