/***************Suavizado y remuestreo*****************/

//Moving Least Squares (MLS) el metodo de reconstruccion de superficie se 
//usa para suavizar y volver a muestrear los datos ruidosos.
#include <string>
#include <iostream>
#include <pcl/point_types.h>
#include <pcl/io/pcd_io.h>
#include <pcl/io/ply_io.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/surface/mls.h>

#include <ctime>
int
main (int argc, char *argv[])
{
  
  // Carga la nube de puntos la cual se desea suavizar 
  pcl::PointCloud<pcl::PointXYZ>::Ptr cloud (new pcl::PointCloud<pcl::PointXYZ> ());
  pcl::io::loadPLYFile ("point_cloud/noise_point_cloud.ply", *cloud);

  // Crea  un  KD-Tree(es una estructura de datos de particionado del espacio que organiza 
  //los puntos en un Espacio euclideo de k dimensiones.)
  pcl::search::KdTree<pcl::PointXYZ>::Ptr tree (new pcl::search::KdTree<pcl::PointXYZ>);

  // La salida tiene el tipo PointNormal para almacenar las normales calculadas por MLS
  pcl::PointCloud<pcl::PointNormal> mls_points;

  // Init object (second point type is for the normals, even if unused)
  pcl::MovingLeastSquares<pcl::PointXYZ, pcl::PointNormal> mls;
  
  // Parametros por defecto
  mls.setComputeNormals (true);
  mls.setInputCloud (cloud);
  mls.setPolynomialFit (true);
  mls.setSearchMethod (tree);

  // Parametro
 
  mls.setSearchRadius(atof(argv[1])); // 0.3

  // Reconstruccion nube de puntos
  mls.process (mls_points);

  // Se almecena la nube de puntos suavizada
  pcl::io::savePLYFile ("point_cloud/outliers_point_cloud.ply", mls_points);
}
