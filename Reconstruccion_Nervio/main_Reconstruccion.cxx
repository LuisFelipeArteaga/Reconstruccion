#include <vtkSmartPointer.h>

#include <vtkBYUReader.h>
#include <vtkPLYReader.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkPolyDataReader.h>
#include <vtkOBJReader.h>
#include <vtkSTLReader.h>
#include <vtkPointSource.h>

#include <vtkPCANormalEstimation.h>
#include <vtkSignedDistance.h>
#include <vtkExtractSurface.h>
#include <vtkExtractSurface.h>
#include <vtkPoissonReconstruction.h>
#include <vtkPowerCrustSurfaceReconstruction.h>

#include <vtkPointData.h>

#include <vtkCoordinate.h>
#include <vtkTextMapper.h>
#include <vtkTextProperty.h>
#include <vtkActor2D.h>

#include <vtkPolyLine.h>
#include <vtkActor.h>
#include <vtkPolyDataMapper2D.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkProperty2D.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkCamera.h>
#include <vtkSTLWriter.h>

#include <vtkNamedColors.h>
#include <vtkTimerLog.h>
#include <vtksys/SystemTools.hxx>

#include <sstream>
#include <vector>
#include <string>
#include <ctime>

#include <sys/types.h>


// Funciones para visualizar las superficies
namespace
{
void MakeViewportGrid(
  std::vector<vtkSmartPointer<vtkRenderer>> &renderers,
  unsigned int renderersize, // Tamaño de ventana
  unsigned int xGridDimensions,  
  unsigned int yGridDimensions);
void ViewportBorder(vtkSmartPointer<vtkRenderer> &renderer,
                    double *color,
                    bool last = false);
}

int main(int argc, char *argv[])
{
  // Se genera un objeto de datos para almacenar la nube de puntos 
  // vtkSmartPointer es una plantilla de clase que proporciona conversión automática para objetos 
  vtkSmartPointer<vtkPolyData> polyData;

  if(argc < 2)
  {
  	/*************Se realiza el registro de las imagenes  de profundidad*****************/
  	std::string number_images = "10"; // Divide el numero total de imagenes en este number_images. int
  	std::string downsample = "10"; // Submuestrea la nube de puntos. int
  	std::string boundary = "5000"; // Se utliza para selecionar el objeto de interes en las imagenes de profundidad. int
  	std::stringstream command1;
  	command1 <<"python image_registration_python/registered_points.py "<<" "<<number_images<<" "<<downsample<<" "<<boundary;
  	

  	std::cout<<"Realizando registro de imagenes..."<<std::endl;
    // Se corre el script en python que realiza el registro
    system (&command1.str()[0]);


    // Los siguientes pasos deben realizarse cunando la nube de puntos tiene un numero considerable de puntos.
    // (Mayor a 10.000 puntos)

    /*******************Suaviza la nube de puntos*********************/

    // Como el registro de las imagenes no es preciso, es necesario suvaizar la nube de puntos.
    std::string SearchRadius = " 0.03"; // Radio de busqueda  
    std::stringstream command2;
    command2 << "./smoothing/build/resampling"<<" "<<SearchRadius;

    std::cout<<"Suavizando Nube de Puntos..."<<std::endl;
    // Se corre el scrip en C++ que realiza el suavizado. Con la libreria PCL(Point Clound Library)
    system (&command2.str()[0]);
    

    /*********************Removiendo Outliers*******************/
    // Remover partes de objetos que no se quieren en la nube de puntos 
    std::string RemovingRadius = " 30.0"; // Radio 
    std::string NumberNeighbors = " 25.0"; // Numero de vecinos
    std::stringstream command3;
    command3 << "./removing_outliers/build/RadiusOutlierRemoval"<<" "<<RemovingRadius<<" "<<NumberNeighbors;

    std::cout<<"Removiendo outliers de nube de Puntos..."<<std::endl;
    // Se corre el script en C++ que remove los outliers. Con la librreria VTK.
    system (&command3.str()[0]);


    // Carga la nube de puntos despues de suavisar la nube y remover outliers	
    std::string extension = vtksys::SystemTools::GetFilenameExtension("point_cloud/point_cloud.ply");
    vtkSmartPointer<vtkPLYReader> reader =
      vtkSmartPointer<vtkPLYReader>::New();
    reader->SetFileName ("point_cloud/point_cloud.ply");
    reader->Update();
    polyData = reader->GetOutput();
    std::cout << "   Numero de puntos de la nube de puntos: " << polyData->GetNumberOfPoints() << std::endl;
  }
  else
  {
  	/*************Carga nube de puntos***************/

  	// Carga la nube de puntos procesada; es decir, La nube de puntos no es necesario realizar 
  	// suavizado y remover outliers. 
  	// El numero de la nube de puntos no es mayor a 10.000 puntos

  	std::cout<<"Cargando nube de Puntos..."<<std::endl;
	std::string extension = vtksys::SystemTools::GetFilenameExtension(argv[1]);
	vtkSmartPointer<vtkPLYReader> reader =
	  vtkSmartPointer<vtkPLYReader>::New();
	reader->SetFileName (argv[1]);
	reader->Update();
	polyData = reader->GetOutput();
	std::cout << "   Numero de Puntos: " << polyData->GetNumberOfPoints() << std::endl;
  }
  

  std::vector<std::string> strVec;
  strVec.push_back("Método 1");
  strVec.push_back("Método 2");
  strVec.push_back("Método 3");


  ////////////////////////////////// Metodo 1 /////////////////////////////////////////
  // ExtractSurface: Extrae la isosuperficie de cruce por cero a partir de una función de distancia con signo.
  // Primero, calcula las normales a partir de los puntos, luego crea un campo de distancia con signo, 
  // seguido de la extracción de superficie del conjunto de nivel cero del campo de distancia.

  unsigned  t01 , t11;
  t01 = clock();
  
  std::cout<<"Reconstruyendo superficie con Método 1..."<<std::endl;

  double bounds[6];
  polyData->GetBounds(bounds);
  double range[3];
  for (int i = 0; i < 3; ++i)
  {
    range[i] = bounds[2*i + 1] - bounds[2*i];
  }

  int sampleSize = polyData->GetNumberOfPoints()*.00005;
  if (sampleSize < 10)
  {
    sampleSize = 10;
  }
  //vtkSignedDistance es un filtro que calcula las distancias firmadas sobre un volumen desde una nube de puntos de entrada.
  //La nube de puntos de entrada debe tener puntos normales definidos.
  vtkSmartPointer<vtkSignedDistance> distance =
    vtkSmartPointer<vtkSignedDistance>::New();
  if (polyData->GetPointData()->GetNormals())
  {
    distance->SetInputData (polyData);
  }
  else
  {
  	// vtkPCANormalEstimation: Genera puntos normales utilizando PCA, se estima un plano tangente local alrededor de cada punto de muestra  considerando 
  	//un pequeño vecindario de puntos alrededor , y ajustando un plano al vecindario a través de PCA.
    vtkSmartPointer<vtkPCANormalEstimation> normals =
      vtkSmartPointer<vtkPCANormalEstimation>::New();
    normals->SetInputData (polyData);
    normals->SetSampleSize(sampleSize);
    normals->SetNormalOrientationToGraphTraversal();
    normals->FlipNormalsOn();
    distance->SetInputConnection (normals->GetOutputPort());
  }
  int dimension = 256;
  double radius;
  radius = std::max(std::max(range[0], range[1]), range[2])
    / static_cast<double>(dimension) * 4; //

  distance->SetRadius(radius);
  distance->SetDimensions(dimension, dimension, dimension);
  distance->SetBounds(
    bounds[0] - range[0] * .10,
    bounds[1] + range[0] * .10,
    bounds[2] - range[1] * .10,
    bounds[3] + range[1] * .10,
    bounds[4] - range[2] * .10,
    bounds[5] + range[2] * .10);

  vtkSmartPointer<vtkExtractSurface> surface1 =
    vtkSmartPointer<vtkExtractSurface>::New();
  surface1->SetInputConnection (distance->GetOutputPort());
  surface1->SetRadius(radius * 10);
  surface1->Update();

  // Almacenar superficie
  vtkSmartPointer<vtkSTLWriter> stlWriter1 =
    vtkSmartPointer<vtkSTLWriter>::New();
  stlWriter1->SetFileName("surface/surface1.stl");
  stlWriter1->SetInputConnection(surface1->GetOutputPort());
  stlWriter1->Write();

  t11 = clock();
  float time1 = (float)(t11-t01)/CLOCKS_PER_SEC;

  ///////////////////////////// Metodo 2 ///////////////////////////////////
  unsigned  t02 , t12;
  t02 = clock();

  // El algoritmo de reconstrucción de superficies PowerCrust reconstruye superficies a 
  // partir de datos de puntos no organizados.

  std::cout<<"Reconstruyendo superficie con Metodo 2..."<<std::endl;
  vtkSmartPointer<vtkPowerCrustSurfaceReconstruction> surface2 =
    vtkSmartPointer<vtkPowerCrustSurfaceReconstruction>::New();
  surface2->SetInputData (polyData);


  //Almacenar superficie
  vtkSmartPointer<vtkSTLWriter> stlWriter2 =
    vtkSmartPointer<vtkSTLWriter>::New();
  stlWriter2->SetFileName("surface/surface2.stl");
  stlWriter2->SetInputConnection(surface2->GetOutputPort());
  stlWriter2->Write();

  t12 = clock();
  float time2 = (float)(t12-t02)/CLOCKS_PER_SEC;


 /////////////////////////// Metodo 3 ///////////////////////////////////	
  unsigned  t03 , t13;
  t03 = clock();
  // Recontruye la  superficie a partir de nube de puntos.
  // El metodo 3 fue implemnatado en la libreria PCL.
  std::cout<<"Reconstruyendo superficie con Método 3..."<<std::endl;

  std::string triangular_radius = "40.0"; // Establezca la distancia máxima entre los puntos conectados.
  std::string Mu = "7.0";// Especifica la distancia máxima aceptable para un punto a considerar vecino.
  std::string MaximunNeighbors = "200";  // Define en cuántos vecinos se buscan.
  std::stringstream command4;
  
  if(argc < 2)
  {
  	command4 << "./surface_triangulo/build/greedy_projection point_cloud/point_cloud.ply" <<" "<<triangular_radius<<" "<<Mu<<" "<<MaximunNeighbors;
  	system(&command4.str()[0]);
  }
  else
  { 
    std::string name2 = argv [1];
  	command4<<"./surface_triangulo/build/greedy_projection"<<" "<<name2<<" "<<triangular_radius<<" "<<Mu<<" "<<MaximunNeighbors;
  	system(&command4.str()[0]);
  }

  
  // Carga la superficie reconstruidad por el metodo 3
  vtkSmartPointer<vtkPolyData> surface3;
  std::string extension2 = vtksys::SystemTools::GetFilenameExtension(std::string("surface/surface3.obj"));
  vtkSmartPointer<vtkOBJReader> surfaceObj =
    vtkSmartPointer<vtkOBJReader>::New();
  surfaceObj->SetFileName ("surface/surface3.obj");
  surfaceObj->Update();
  surface3 = surfaceObj->GetOutput();
  
  t13 = clock();
  float time3 = (float)(t13-t03)/1000;



  //////////////////////// Graficas /////////////////////////////////////

  // Almacena las superficies de cada metodo en un vector
  std::vector<vtkSmartPointer<vtkPolyDataAlgorithm> > surfaceObjects;
  surfaceObjects.push_back(surface1);
  surfaceObjects.push_back(surface2);
  surfaceObjects.push_back(surface2); 


  std::vector<vtkSmartPointer<vtkRenderer>> renderers;

  // vtkNamedColors: Una clase con colores y sus nombres.
  vtkSmartPointer<vtkNamedColors> colors =
    vtkSmartPointer<vtkNamedColors>::New();

  // vtkRenderWindow es un objeto abstracto para especificar el comportamiento de una ventana de representación.
  vtkSmartPointer<vtkRenderWindow> renderWindow =
    vtkSmartPointer<vtkRenderWindow>::New();


  // vtkCamera es una cámara virtual para renderización 3D
  vtkSmartPointer<vtkCamera> camera =
    vtkSmartPointer<vtkCamera>::New();


   // For que recorre las tres superficies a visualizar
  for (size_t i = 0; i < surfaceObjects.size(); ++i)
  {
    //vtkPolyDataMapper es una clase que asigna datos poligonales PolyData a gráficos primitivos
    vtkSmartPointer<vtkPolyDataMapper> surfaceMapper =
      vtkSmartPointer<vtkPolyDataMapper>::New(); 
    // Se extrae la superficie del vector	
    if(i == 2 )
    {
      surfaceMapper->SetInputData(surface3);
    } 
    else
    {
      surfaceMapper->SetInputConnection(surfaceObjects[i]->GetOutputPort());	
    }
    // vtkProperty es un objeto que representa la iluminación y otras propiedades de superficie de un objeto geométrico.
    vtkSmartPointer<vtkProperty> back =
      vtkSmartPointer<vtkProperty>::New();
    back->SetDiffuseColor(colors->GetColor3d("alizarin_crimson").GetData());
    back->SetSpecular(.6);
    back->SetSpecularPower(50.0);

    //vtkActor se usa para representar una entidad en una escena de renderizado.
    vtkSmartPointer<vtkActor> surfaceActor =
      vtkSmartPointer<vtkActor>::New();
    surfaceActor->SetMapper(surfaceMapper);
    surfaceActor->GetProperty()->SetDiffuseColor(colors->GetColor3d("alizarin_crimson").GetData());
    surfaceActor->GetProperty()->SetSpecular(.6);
    surfaceActor->GetProperty()->SetSpecularPower(50.0);
    surfaceActor->SetBackfaceProperty(back);
    
    // vtkRenderer proporciona una especificación abstracta para renderizadores
    vtkSmartPointer<vtkRenderer> renderer =
      vtkSmartPointer<vtkRenderer>::New();
    renderer->AddActor(surfaceActor);
    renderer->SetBackground(colors->GetColor3d("cobalt").GetData());
    renderer->SetActiveCamera(camera);
    renderer->GetActiveCamera()->SetPosition (-1, 0, 0);
    renderer->GetActiveCamera()->SetFocalPoint (0, 1, 0);
    renderer->GetActiveCamera()->SetViewUp (0, 0, 1);
    renderer->GetActiveCamera()->Dolly(5);
    renderer->ResetCamera();
    renderer->ResetCameraClippingRange();
    renderers.push_back(renderer);
    renderWindow->AddRenderer(renderer);

    //vtkTextProperty es un objeto que representa propiedades de texto.
    vtkSmartPointer<vtkTextProperty> textProperty =
      vtkSmartPointer<vtkTextProperty>::New();
    textProperty->SetFontSize(15);
    textProperty->SetJustificationToCentered();


    std::stringstream ss;
    ss << strVec[i] << std::endl;
    if(i == 0){
      ss << "Number of Polys: " << surfaceObjects[i]->GetOutput()->GetNumberOfPolys() << std::endl;
      ss << "Time: " << time1 << "s" << std::endl;
    }
    if(i == 1){
      ss << "Number of Polys: " << surfaceObjects[i]->GetOutput()->GetNumberOfPolys() << std::endl;
      ss << "Time: " << time2 << "s" << std::endl;
    }
    if(i == 2){
      ss << "Number of Polys: " << surface3->GetNumberOfPolys() << std::endl;
      ss << "Time: " << time3 << "s" << std::endl;
    }
  
    

    //vtkTextMapper proporciona soporte de anotación de texto 2D para VTK.
    vtkSmartPointer<vtkTextMapper> textMapper =
      vtkSmartPointer<vtkTextMapper>::New();
    textMapper->SetInput(ss.str().c_str());
    textMapper->SetTextProperty(textProperty);

    //vtkActor2D es similar a vtkActor,
    vtkSmartPointer<vtkActor2D> textActor =
      vtkSmartPointer<vtkActor2D>::New();
    textActor->SetMapper(textMapper);
    textActor->SetPosition(200, 0);
    textActor->GetProperty()->SetLineWidth(4.0); // Line Width 

    renderer->AddViewProp(textActor);
  }

  unsigned int rendererSize = 450;
  unsigned int xGridDimensions = 3;
  unsigned int yGridDimensions = 1;
  renderWindow->SetSize(
    rendererSize * xGridDimensions,
    rendererSize * yGridDimensions);

  MakeViewportGrid(renderers,
                   rendererSize, 
                   xGridDimensions,
                   yGridDimensions);
  for (size_t i = 0; i < renderers.size(); ++i)
  {
    ViewportBorder(renderers[i],
                  colors->GetColor3d("Blackclose").GetData(),
                  i == renderers.size() - 1);
  }

  //vtkRenderWindowInteractor proporciona un mecanismo de interacción independiente de la plataforma para eventos de mouse/key/tiempo
  vtkSmartPointer<vtkRenderWindowInteractor> iren =
    vtkSmartPointer<vtkRenderWindowInteractor>::New();
  iren->SetRenderWindow(renderWindow);

  iren->Initialize();
  iren->Start();

  return EXIT_SUCCESS;

}

namespace {
void MakeViewportGrid(
  std::vector<vtkSmartPointer<vtkRenderer>> &renderers,
  unsigned int rendererSize,
  unsigned int xGridDimensions,
  unsigned int yGridDimensions)
{
  //
  //Configuración de viewports para los renderizadores
  int index = 0;
  for (int row = 0; row < static_cast<int>(yGridDimensions); row++)
  {
    for (int col = 0; col < static_cast<int>(xGridDimensions); col++)
    {
      int index = row * xGridDimensions + col;
      double viewport[4] = {
        static_cast<double>(col) * rendererSize /
        (xGridDimensions * rendererSize),
        static_cast<double>(yGridDimensions - (row + 1)) * rendererSize /
        (yGridDimensions * rendererSize),
        static_cast<double>(col + 1) * rendererSize /
        (xGridDimensions * rendererSize),
        static_cast<double>(yGridDimensions - row) * rendererSize /
        (yGridDimensions * rendererSize)};
      renderers[index]->SetViewport(viewport);
    }
  }
}
// dibujar los bordes de la ventana gráfica de un renderizador
void ViewportBorder(vtkSmartPointer<vtkRenderer> &renderer,
                    double *color,
                    bool last)
{
  // los puntos comienzan en la esquina superior derecha y avanzan en sentido antihorario
  vtkSmartPointer<vtkPoints> points =
    vtkSmartPointer<vtkPoints>::New();
  points->SetNumberOfPoints(4);
  points->InsertPoint(0, 1, 1, 0);
  points->InsertPoint(1, 0, 1, 0);
  points->InsertPoint(2, 0, 0, 0);
  points->InsertPoint(3, 1, 0, 0);

  // vtkCellArray es un objeto de soporte que representa explícitamente la conectividad celular.
  vtkSmartPointer<vtkCellArray> cells =
    vtkSmartPointer<vtkCellArray>::New();
  cells->Initialize(); 

  vtkSmartPointer<vtkPolyLine> lines =
    vtkSmartPointer<vtkPolyLine>::New();
  if (last)
  {
    lines->GetPointIds()->SetNumberOfIds(5);
  }
  else
  {
  lines->GetPointIds()->SetNumberOfIds(4);
  }
  for(unsigned int i = 0; i < 4; ++i)
  {
    lines->GetPointIds()->SetId(i,i);
  }
  if (last)
  {
    lines->GetPointIds()->SetId(4, 0);
  }
  cells->InsertNextCell(lines);

 
  vtkSmartPointer<vtkPolyData> poly =
    vtkSmartPointer<vtkPolyData>::New();
  poly->Initialize(); 
  poly->SetPoints(points); 
  poly->SetLines(cells); 

  // usa coordenadas de ventana gráfica normalizadas ya que son independientes del tamaño de la ventana
  vtkSmartPointer<vtkCoordinate> coordinate =
    vtkSmartPointer<vtkCoordinate>::New();
  coordinate->SetCoordinateSystemToNormalizedViewport(); 

  vtkSmartPointer<vtkPolyDataMapper2D> mapper =
    vtkSmartPointer<vtkPolyDataMapper2D>::New();
  mapper->SetInputData(poly); 
  mapper->SetTransformCoordinate(coordinate); 

  vtkSmartPointer<vtkActor2D> actor =
    vtkSmartPointer<vtkActor2D>::New();
  actor->SetMapper(mapper); 
  actor->GetProperty()->SetColor(color);
  actor->GetProperty()->SetLineWidth(4.0); 

  renderer->AddViewProp(actor);
}
}