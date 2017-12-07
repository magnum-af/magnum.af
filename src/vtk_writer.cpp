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


//3D Image data Cell Centered
//Optimization should avoid generation of two vtkImageData objects
void vti_writer(const af::array field, const Mesh& mesh, std::string outputname){
  
    //print("FLAT",af::flat(field));
    double* host_a = field.host<double>();
    //for(int i=0; i < field.dims(0)* field.dims(1)* field.dims(2) * field.dims(3); i++){
    //    std::cout<<"host: "<< i << " = " << host_a[i]<<std::endl;
    //}
    //af::array Test(field.dims(0), field.dims(1), field.dims(2) , field.dims(3),host_a);
    //af::print("Test=",Test);
  
    std::cout<<"vtk_writer: af_to_vti: Number of af::array elements:"<< field.dims(0)* field.dims(1)* field.dims(2) * field.dims(3)<<std::endl;
  
    vtkSmartPointer<vtkImageData> imageDataPointCentered = vtkSmartPointer<vtkImageData>::New();
    imageDataPointCentered->SetDimensions(field.dims(0), field.dims(1), field.dims(2));
    imageDataPointCentered->SetSpacing(mesh.dx,mesh.dy,mesh.dz);
    imageDataPointCentered->SetOrigin(0,0,0);
    #if VTK_MAJOR_VERSION <= 5
       imageDataPointCentered->SetNumberOfScalarComponents(field.dims(3));
       imageDataPointCentered->SetScalarTypeToDouble();
    #else
       imageDataPointCentered->AllocateScalars(VTK_DOUBLE, field.dims(3));
    #endif
    int* dims = imageDataPointCentered->GetDimensions();
  
    for (int z = 0; z < dims[2]; z++)
      {
      for (int y = 0; y < dims[1]; y++)
        {
        for (int x = 0; x < dims[0]; x++)
          {
          for (int im=0; im < field.dims(3); im++)
            {
            double* pixel = static_cast<double*>(imageDataPointCentered->GetScalarPointer(x,y,z));
            pixel[im] = host_a[x+dims[0]*(y+dims[1]*(z+ dims[2] * im))];
            //int idx = x+dims[0]*(y+dims[1]*(z+ dims[2] * im));
            //std::cout<<"idx= "<< idx<< "\t" <<host_a[idx]<<std::endl;                                               
            }
          }
        }
      }

    vtkSmartPointer<vtkImageData> imageDataCellCentered = vtkSmartPointer<vtkImageData>::New();
    imageDataCellCentered->SetDimensions(field.dims(0)+1, field.dims(1)+1, field.dims(2)+1);
    imageDataCellCentered->SetOrigin(0,0,0);
    imageDataCellCentered->SetSpacing(mesh.dx,mesh.dy,mesh.dz);
    imageDataCellCentered->GetCellData()->SetScalars (imageDataPointCentered->GetPointData()->GetScalars());
  
    vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
    writer->SetFileName((outputname.append(".vti")).c_str());
    #if VTK_MAJOR_VERSION <= 5
        writer->SetInputConnection(imageDataCellCentered->GetProducerPort());
    #else
        writer->SetInputData(imageDataCellCentered);
    #endif
    writer->Write();
}

//https://www.vtk.org/gitweb?p=VTK.git;a=blob;f=Examples/DataManipulation/Cxx/Arrays.cxx

    //USEAGE:
    //array Aout = array();
    //Mesh meshout = Mesh(NaN,NaN,NaN,NaN,NaN,NaN);
    //vti_to_af(Aout,meshout,"/home/pth/git/pth-mag/Data/Testing/minit.vti");
    //print("Aout",Aout);
    //std::cout << "  " << meshout.n0 << "  " << meshout.n1 << "  " << meshout.n2 << "  " << meshout.dx << "  " << meshout.dy << "  " << meshout.dz <<  std::endl;
void vti_reader(af::array& field, Mesh& mesh, std::string filepath){
    int dim4th = 3;//This is the number of the components of the 3D Field (until now only 3)
    //Candidate is imageData->GetScalarSize()<<std::endl;
    vtkSmartPointer<vtkXMLImageDataReader> reader = vtkSmartPointer<vtkXMLImageDataReader>::New();
    reader->SetFileName(filepath.c_str());
    reader->Update(); 
    vtkSmartPointer<vtkImageData> imageData = vtkSmartPointer<vtkImageData>::New();
    imageData=reader->GetOutput();

    vtkSmartPointer<vtkDoubleArray> temp = vtkSmartPointer<vtkDoubleArray>::New();
    imageData->GetCellData()->GetScalars()->GetData(0,imageData->GetNumberOfCells()-1,0,dim4th-1,temp);

    int* dims_reduced = imageData->GetDimensions();
    double* spacing = imageData->GetSpacing();
    for(int i=0; i < dim4th; i++){
        dims_reduced[i]--;
    }
    mesh=Mesh(dims_reduced[0],dims_reduced[1],dims_reduced[2],spacing[0],spacing[1],spacing[2]);
    double* A_host = NULL;
    A_host = new double[dim4th*imageData->GetNumberOfCells()];

    //Two Versions to sort CellData to double* A_host 
    //Perform equally both on CPU and OpenCL

    //VERSION SORT WITH ARRAYFIRE
    for(int i=0; i < dim4th* imageData->GetNumberOfCells(); i++){
        A_host[i]=temp->GetValue(i);
    }
    af::array A(dim4th*imageData->GetNumberOfCells(),1,1,1,A_host);
    delete [] A_host;
    A=af::moddims(A,af::dim4(dim4th,dims_reduced[0],dims_reduced[1],dims_reduced[2]));
    A=af::reorder(A,1,2,3,0);
    field=A;

    ////VERSION SORT WITH C++ LOOP
    //int A_host_idx=0;
    //for(int i=0; i < dim4th; i++){
    //    for (int j = i; j < dim4th * imageData->GetNumberOfCells(); j = j + dim4th){
    //        A_host[A_host_idx]=temp->GetValue(j);
    //        A_host_idx++;
    //    }
    //}
    //af::array A(dim4th*imageData->GetNumberOfCells(),1,1,1,A_host);
    //delete [] A_host;
    //A=af::moddims(A,af::dim4(dims_reduced[0],dims_reduced[1],dims_reduced[2],dim4th));
    //field=A;
}


//3D Image data
void atomistic_vti_writer(const af::array field, const Mesh& mesh, std::string outputname){
 
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
 
//TODO    std::cout<<"---->VTK: GetScalarSize()"<<imageData->GetScalarSize()<<std::endl;
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

void atomistic_vti_reader(af::array& field, Mesh& mesh, std::string filepath){
    int dim4th = 3;//This is the number of the components of the 3D Field (until now only 3)
    //Candidate is imageData->GetScalarSize()<<std::endl;
    vtkSmartPointer<vtkXMLImageDataReader> reader = vtkSmartPointer<vtkXMLImageDataReader>::New();
    reader->SetFileName(filepath.c_str());
    reader->Update(); 
    vtkSmartPointer<vtkImageData> imageData = vtkSmartPointer<vtkImageData>::New();
    imageData=reader->GetOutput();
    //std::cout<<reader->GetOutput();

    //vtkSmartPointer<vtkDoubleArray> temp = vtkSmartPointer<vtkDoubleArray>::New();
    //imageData->GetCellData()->GetScalars()->GetData(0,imageData->GetNumberOfCells()-1,0,dim4th-1,temp);

    int* dims = imageData->GetDimensions();
    double* spacing = imageData->GetSpacing();
    mesh=Mesh(dims[0],dims[1],dims[2],spacing[0],spacing[1],spacing[2]);
    double* A_host = NULL;
    A_host = new double[dim4th*imageData->GetNumberOfPoints()];
    for (int z = 0; z < dims[2]; z++)
      {
      for (int y = 0; y < dims[1]; y++)
        {
        for (int x = 0; x < dims[0]; x++)
          {
          for (int im=0; im < dim4th; im++)
            {
            double* pixel = static_cast<double*>(imageData->GetScalarPointer(x,y,z));
            A_host[x+dims[0]*(y+dims[1]*(z+ dims[2] * im))] = pixel[im] ;
            int idx = x+dims[0]*(y+dims[1]*(z+ dims[2] * im));
            std::cout<<"idx= "<< idx<< "\t" <<A_host[idx]<<std::endl;                                               
            }
          }
        }
      }

    //for(int i=0; i < dim4th* imageData->GetNumberOfPoints(); i++){
    //    std::cout<<"A_host[i] = "<<A_host[i]<<std::endl;
    //}
    //std::cout<<imageData->GetNumberOfPoints()<<std::endl;
    //Two Versions to sort CellData to double* A_host 
    //Perform equally both on CPU and OpenCL

    ////VERSION SORT WITH ARRAYFIRE
    //for(int i=0; i < dim4th* imageData->GetNumberOfCells(); i++){
    //    A_host[i]=temp->GetValue(i);
    //}
    af::array A(dim4th*imageData->GetNumberOfPoints(),1,1,1,A_host);
    delete [] A_host;
    A=af::moddims(A,af::dim4(dims[0],dims[1],dims[2],dim4th));
    //A=af::reorder(A,1,2,3,0);
    field=A;
}
