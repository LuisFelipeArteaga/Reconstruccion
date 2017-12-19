/****************Remover outliers************************/

// Elimina puntos aislados; es decir, aquellos puntos que tienen pocos vecinos dentro de un radio especificado.
// El usuario debe especificar el radio que define la region local, asi como el umbral de 
// aislamiento es decir, el numero de puntos vecinos requeridos para que el punto se considere aislado.
#include <string>
#include <vtkSmartPointer.h>
#include <vtkRadiusOutlierRemoval.h>
 
#include <vtkPLYReader.h>
#include <vtkPLYWriter.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkPointSource.h>
 
#include <vtkSphereSource.h>
#include <vtkGlyph3D.h>

#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
 
#include <vtksys/SystemTools.hxx>


 
int main (int argc, char *argv[])
{

// Carga la nube de puntos en formato PLY 
  vtkSmartPointer<vtkPolyData> polyData;
  std::string extension = vtksys::SystemTools::GetFilenameExtension("point_cloud/outliers_point_cloud.ply");
  vtkSmartPointer<vtkPLYReader> reader =
    vtkSmartPointer<vtkPLYReader>::New();
  reader->SetFileName ("point_cloud/outliers_point_cloud.ply");
  reader->Update();
  polyData = reader->GetOutput();


  double bounds[6];
  polyData->GetBounds(bounds);
  double range[3];
  for (int i = 0; i < 3; ++i)
  {
    range[i] = bounds[2*i + 1] - bounds[2*i];
  }
  vtkSmartPointer<vtkRadiusOutlierRemoval> removal =
    vtkSmartPointer<vtkRadiusOutlierRemoval>::New();
  removal->SetInputData(polyData);
  removal->SetRadius(range[0]/atof(argv[1])); // Radio 
  removal->SetNumberOfNeighbors(atof(argv[2])); // Numero de vecinos
  removal->GenerateOutliersOn();
  removal->Update();
  std::cout << "   Numero de puntos removidos: " << removal->GetNumberOfPointsRemoved() << std::endl;
 
 //Guardar la nube de puntos sin outliers
  vtkSmartPointer<vtkPLYWriter> plyWriter1 =
  vtkSmartPointer<vtkPLYWriter>::New();
  plyWriter1->SetFileName("point_cloud/point_cloud.ply");
  plyWriter1->SetInputConnection(removal->GetOutputPort());
  plyWriter1->Write();
  plyWriter1->Update();

 
  return EXIT_SUCCESS;
}
 



