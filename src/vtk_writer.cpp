#include "vtk_writer.hpp"

void af_to_vtk(const af::array field, const Mesh& mesh, std::string outputname){
//void af_to_vtk(const State& a){

  const double d[3]={mesh.dx,mesh.dy,mesh.dz};
  af::dim4 dims(field.dims(0),field.dims(1),field.dims(2),field.dims(3));

  double* host_a = field.host<double>();

  std::cout<<"vtk_writer: Number of points:"<< dims[0]*dims[1]*dims[2]*dims[3]<<std::endl;

  //------------------------------------------------------------------------------
  //VKT grid
  //------------------------------------------------------------------------------
  vtkSmartPointer<vtkRectilinearGrid> grid =
    vtkSmartPointer<vtkRectilinearGrid>::New();

  grid->SetDimensions(dims[0], dims[1], dims[2]);

  vtkDataArray* coords[3];

  // compute & populate coordinate vectors
  for(int i=0; i < dims[3]; ++i)//i^= x,y,z
  {
    coords[i] = vtkDataArray::CreateDataArray(VTK_DOUBLE);
    coords[i]->SetNumberOfTuples(dims[i]);

    for(int j=0; j < dims[i]; ++j)
    {
      double val = (double)j*d[i];
      coords[ i ]->SetTuple( j, &val );
    } // END for all points along this dimension

  } // END for all dimensions

  grid->SetXCoordinates( coords[0] );
  grid->SetYCoordinates( coords[1] );
  grid->SetZCoordinates( coords[2] );
  coords[0]->Delete();
  coords[1]->Delete();
  coords[2]->Delete();

  //VKT value grid
  vtkSmartPointer<vtkRectilinearGrid> value_grid =
    vtkSmartPointer<vtkRectilinearGrid>::New();
  value_grid->SetDimensions(1,1,1);
  vtkDataArray* value_coords[3];

  for(int i=0; i < 3; ++i)//i^= x,y,z
  {
    value_coords[i] = vtkDataArray::CreateDataArray(VTK_DOUBLE);
    value_coords[i]->SetNumberOfTuples(1);
  } // END for all dimensions

  // compute & populate XYZ field
  vtkIdType npoints = grid->GetNumberOfPoints();
  vtkDoubleArray* vkt_xyz = vtkDoubleArray::New();
  vkt_xyz->SetName("m");
  //vkt_xyz->SetName( outputname.c_str() );
  vkt_xyz->SetNumberOfComponents(3);
  vkt_xyz->SetNumberOfTuples( npoints );

  for(vtkIdType pntIdx=0; pntIdx < npoints; ++pntIdx )
  {
    double af_vals[3];
    for(int i=0; i < 3; ++i)
    {
      af_vals[i]=host_a[pntIdx+i*npoints];
      value_coords[i]->SetTuple(0, &af_vals[i] );
    }
      value_grid->SetXCoordinates( value_coords[0] );
      value_grid->SetYCoordinates( value_coords[1] );
      value_grid->SetZCoordinates( value_coords[2] );

      vkt_xyz->SetTuple(pntIdx, value_grid->GetPoint(0) );
  } // END for all points
  grid->GetPointData()->AddArray( vkt_xyz );

//   vtkNew<vtkPointDataToCellData> pd2cd;
//   pd2cd->PassPointDataOn();
//   pd2cd->SetInputDataObject(grid);
//   pd2cd->Update();
//   grid->ShallowCopy(pd2cd->GetOutputDataObject(0));
  //Write into Test.vtk
  //
  vtkRectilinearGridWriter* writer = vtkRectilinearGridWriter::New();
  
  writer->SetFileName((outputname.append(".vtk")).c_str());
  //writer->SetFileName( oss.str().c_str() );
  writer->SetInputData( grid );
  writer->Write();
  writer->Delete();
  vkt_xyz->Delete();

  value_coords[0]->Delete();
  value_coords[1]->Delete();
  value_coords[2]->Delete();

  delete[] host_a;
}


//3D Image data
void af_to_vti(const af::array field, const Mesh& mesh, std::string outputname){
 
   double* host_a = field.host<double>();
 
   std::cout<<"vtk_writer: af_to_vti: Number of points:"<< field.dims(0)* field.dims(1)* field.dims(2) * field.dims(3)<<std::endl;
 
   vtkSmartPointer<vtkImageData> imageData =
     vtkSmartPointer<vtkImageData>::New();
   imageData->SetDimensions(field.dims(0), field.dims(1), field.dims(2));
    imageData->SetSpacing(mesh.dx,mesh.dy,mesh.dz);
 #if VTK_MAJOR_VERSION <= 5
   imageData->SetNumberOfScalarComponents(field.dims(3));
   imageData->SetScalarTypeToDouble();
 #else
   imageData->AllocateScalars(VTK_DOUBLE, field.dims(3));
 #endif
   int* dims = imageData->GetDimensions();
 
   for (int z = 0; z < dims[2]; z++)
     {
     for (int y = 0; y < dims[1]; y++)
       {
       for (int x = 0; x < dims[0]; x++)
         {
         for (int im=0; im < field.dims(3); im++)
           {
           double* pixel = static_cast<double*>(imageData->GetScalarPointer(x,y,z));
           pixel[im] = host_a[x+dims[0]*(y+dims[1]*(z+ dims[2] * im))];
           //int idx = x+dims[0]*(y+dims[1]*(z+ dims[2] * im));
           //std::cout<<"idx= "<< idx<< "\t" <<host_a[idx]<<std::endl;                                               
           }
         }
       }
     }
 
   vtkSmartPointer<vtkXMLImageDataWriter> writer =
     vtkSmartPointer<vtkXMLImageDataWriter>::New();
   writer->SetFileName((outputname.append(".vti")).c_str());
 #if VTK_MAJOR_VERSION <= 5
   writer->SetInputConnection(imageData->GetProducerPort());
 #else
   writer->SetInputData(imageData);
 #endif
   writer->Write();
}
