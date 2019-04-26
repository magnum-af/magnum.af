#include "vtk_IO.hpp"


//3D vtkImageData vtkCellData writer
//Optimization should avoid generation of two vtkImageData objects
void vti_writer_micro(const af::array field, const Mesh& mesh, std::string outputname){
  
    double* host_a = field.host<double>();
  
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
    //std::cout<<"vti_writer_micro: Writing vtkCellData with "<< field.dims(0)* field.dims(1)* field.dims(2) 
    //    << " Cells in file "<<outputname<<std::endl;
    #if VTK_MAJOR_VERSION <= 5
         writer->SetInputConnection(imageDataCellCentered->GetProducerPort());
    #else
         writer->SetInputData(imageDataCellCentered);
    #endif
    writer->Write();
}


//3D vtkImageData vtkPointData writer
void vti_writer_atom(const af::array field, const Mesh& mesh, std::string outputname){
  
    double* host_a = field.host<double>();
  
    vtkSmartPointer<vtkImageData> imageData = vtkSmartPointer<vtkImageData>::New();
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
            }
          }
        }
      }
    vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
    writer->SetFileName((outputname.append(".vti")).c_str());
    //std::cout<<"vti_writer_atom: Writing vtkPointData with "<< field.dims(0)* field.dims(1)* field.dims(2)  
    //     << " Points in file "<<outputname<<std::endl;
    #if VTK_MAJOR_VERSION <= 5
        writer->SetInputConnection(imageData->GetProducerPort());
    #else
        writer->SetInputData(imageData);
    #endif
    writer->Write();
}


void vti_reader(af::array& field, Mesh& mesh, std::string filepath){
    int dim4th = 3;//This is the number of the components of the 3D Field (until now only 3)
    vtkSmartPointer<vtkXMLImageDataReader> reader = vtkSmartPointer<vtkXMLImageDataReader>::New();
    reader->SetFileName(filepath.c_str());
    reader->Update(); 
    vtkSmartPointer<vtkImageData> imageData = vtkSmartPointer<vtkImageData>::New();
    imageData=reader->GetOutput();

    //Check whether input is vtkCellData or vtkPointData
    bool celldata=false;// If true, vtkCellData, if false, vtkPointData
    double* test_pixel = static_cast<double*>(imageData->GetScalarPointer(0,0,0));
    if(test_pixel==NULL){
        std::cout<<"vti_reader: Reading vtkCellData from "<<filepath<<std::endl;
        celldata=true;
    }
    else{
        std::cout<<"vti_reader: Reading vtkPointData from "<<filepath<<std::endl;
    }

    int* dims = imageData->GetDimensions();
    double* spacing = imageData->GetSpacing();

    if(celldata){
        for(int i=0; i < dim4th; i++){
            dims[i]--;
        }
        //IF Celldata:
        vtkSmartPointer<vtkDoubleArray> temp = vtkSmartPointer<vtkDoubleArray>::New();
        imageData->GetCellData()->GetScalars()->GetData(0,imageData->GetNumberOfCells()-1,0,dim4th-1,temp);
        double* A_host = NULL;
        A_host = new double[dim4th*imageData->GetNumberOfCells()];
        //VERSION SORT WITH ARRAYFIRE
        for(int i=0; i < dim4th* imageData->GetNumberOfCells(); i++){
            A_host[i]=temp->GetValue(i);
        }
        af::array A(dim4th*imageData->GetNumberOfCells(),1,1,1,A_host);
        delete [] A_host;
        A=af::moddims(A,af::dim4(dim4th,dims[0],dims[1],dims[2]));
        A=af::reorder(A,1,2,3,0);
        field=A;
        //Two Versions to sort CellData to double* A_host 
        //Perform equally both on CPU and OpenCL

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
    else{
        double* A_host = NULL;
        A_host = new double[dim4th*imageData->GetNumberOfPoints()];
        for (int z = 0; z < dims[2]; z++)
          {
          for (int y = 0; y < dims[1]; y++)
            {
            for (int x = 0; x < dims[0]; x++)
              {
              double* pixel = static_cast<double*>(imageData->GetScalarPointer(x,y,z));
              for (int im=0; im < dim4th; im++)
                {
                A_host[x+dims[0]*(y+dims[1]*(z+ dims[2] * im))] = pixel[im] ;
                }
              }
            }
          }

        af::array A(dim4th*imageData->GetNumberOfPoints(),1,1,1,A_host);
        delete [] A_host;
        A=af::moddims(A,af::dim4(dims[0],dims[1],dims[2],dim4th));
        field=A;
    }
    mesh=Mesh(dims[0],dims[1],dims[2],spacing[0],spacing[1],spacing[2]);
}


////vtkRectilinearGrid writer (currently not needed as we only compute on regular grids
//void vtr_writer(const af::array field, const Mesh& mesh, const std::vector<double> z_spacing, std::string outputname){
//  
//    const double d[3]={mesh.dx,mesh.dy,mesh.dz};
//    const long long int dims[4] = {field.dims(0),field.dims(1),field.dims(2),field.dims(3)};
//  
//    double* host_a = field.host<double>();
//  
//    std::cout<<"vtk_writer: Number of points:"<< dims[0]*dims[1]*dims[2]<<std::endl;
//  
//    //------------------------------------------------------------------------------
//    //VKT grid
//    //------------------------------------------------------------------------------
//    vtkSmartPointer<vtkRectilinearGrid> grid = vtkSmartPointer<vtkRectilinearGrid>::New();
//  
//    grid->SetDimensions(dims[0], dims[1], dims[2]);
//  
//    vtkDataArray* coords[3];
//  
//    // compute & populate coordinate vectors
//    for(int i=0; i < dims[3]; ++i){
//        coords[i] = vtkDataArray::CreateDataArray(VTK_DOUBLE);
//        coords[i]->SetNumberOfTuples(dims[i]);
//    }
//
//    //x
//    for(int j=0; j < dims[0]; ++j){
//        double val = (double)j * d[0];
//        coords[0]->SetTuple(j, &val);
//    }
//    //y
//    std::cout << "test" << std::endl;
//    for(int j=0; j < dims[1]; ++j){
//        double val = (double)j * d[1];
//        std::cout << "test" << j << std::endl;
//        coords[1]->SetTuple(j, &val);
//        std::cout << "test" << j << std::endl;
//    }
//    
//    //z
//    double add_val = 0;
//    for(int j=0; j < dims[2]; ++j)
//    {
//        if ( j == 0){
//            add_val = 0;
//        }
//        else{
//          add_val += z_spacing.at(j-1);
//          std::cout << "val=" << add_val << std::endl;
//        }
//        coords[2]->SetTuple(j, &add_val);
//    }
//    grid->SetXCoordinates( coords[0] );
//    grid->SetYCoordinates( coords[1] );
//    grid->SetZCoordinates( coords[2] );
//
//    // compute & populate XYZ field
//    vtkIdType ncells = grid->GetNumberOfCells();
//
//    vtkDoubleArray* data = vtkDoubleArray::New();
//    data->SetName("m");
//    data->SetNumberOfComponents(3);
//    data->SetNumberOfTuples( ncells );
//    //TODO//data->SetArray();
//  
//    for(vtkIdType cellId=0; cellId < 3; ++cellId){
//    //for(vtkIdType cellId=0; cellId < ncells; ++cellId){
//        for(int i=0; i < 3; ++i){
//            //data->SetValue(cellId + i * ncells, host_a[cellId+i*ncells]);
//            //data->SetValue(3 * cellId + i, 2.);
//            data->SetValue(3 * cellId + i, host_a[cellId+i*ncells]);
//            //data->SetValue(cellId + i * ncells, host_a[cellId+i*ncells]);
//        }
//        //data->SetValue(3 * cellId+0, 0.1);
//        //data->SetValue(3 * cellId+1, 0.2);
//        //data->SetValue(3 * cellId+2, 0.3);
//    }
//  
//    grid->GetCellData()->AddArray(data);
//    coords[0]->Delete();
//    coords[1]->Delete();
//    coords[2]->Delete();
//  
//    vtkRectilinearGridWriter* writer = vtkRectilinearGridWriter::New();
//    writer->SetFileName((outputname.append(".vtr")).c_str());
//    std::cout<<"vtr_writer: Writing vtkRectilinearGrid Data with "<< field.dims(0)* field.dims(1)* field.dims(2) << " Points in file "<<outputname<<std::endl;
//    writer->SetInputData( grid );
//    writer->Write();
//  
//    writer->Delete();
//    data->Delete();
//    std::cout << "finish vtr" << std::endl;
//    delete[] host_a;
//}

//vtkRectilinearGrid writer (currently not needed as we only compute on regular grids
void vtr_writer(const af::array field, const Mesh& mesh, const std::vector<double> z_spacing, std::string outputname, const bool verbose){
  
    const double d[3]={mesh.dx,mesh.dy,mesh.dz};
    const long long int dims[4] = {field.dims(0),field.dims(1),field.dims(2),field.dims(3)};
    //std::cout << "dims = " << dims[3] << std::endl;
  
    double* host_a = field.host<double>();
  
    if(verbose) std::cout<<"vtk_writer: Number of points:"<< dims[0]*dims[1]*dims[2]<<std::endl;
  
    //------------------------------------------------------------------------------
    //VKT grid
    //------------------------------------------------------------------------------
    vtkSmartPointer<vtkRectilinearGrid> grid = vtkSmartPointer<vtkRectilinearGrid>::New();
  
    grid->SetDimensions(dims[0]+1, dims[1]+1, dims[2]+1);
  
    vtkDataArray* coords[3];
  
    // compute & populate coordinate vectors
    for(int i=0; i < 3; ++i)//i^= x,y,z
    {
        coords[i] = vtkDataArray::CreateDataArray(VTK_DOUBLE);
        coords[i]->SetNumberOfTuples(dims[i]+1);
    
        double val = 0;
        for(int j=0; j < dims[i]+1; ++j)
        {
            if (i != 2){
                // x and y
                val = (double)j*d[i];
            }
            else{
                // z dimension
                if ( j == 0){
                    val = 0;
                }
                else{
                  val += z_spacing.at(j-1);
                }

                //sum_over_z_spacings += z_spacing.at(j);
                //val = sum_over_z_spacings;

            }
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
    vtkSmartPointer<vtkRectilinearGrid> value_grid = vtkSmartPointer<vtkRectilinearGrid>::New();
    value_grid->SetDimensions(1,1,1);
    vtkDataArray* value_coords[3];
  
    for(int i=0; i < 3; ++i)//i^= x,y,z
    {
        value_coords[i] = vtkDataArray::CreateDataArray(VTK_DOUBLE);
        value_coords[i]->SetNumberOfTuples(1);
    } // END for all dimensions
  
    // compute & populate XYZ field
    vtkIdType npoints = grid->GetNumberOfCells();
    //KvtkIdType npoints = grid->GetNumberOfPoints();
    vtkDoubleArray* vtk_xyz = vtkDoubleArray::New();
    vtk_xyz->SetName("m");
    //vtk_xyz->SetName( outputname.c_str() );
    vtk_xyz->SetNumberOfComponents(dims[3]);
    vtk_xyz->SetNumberOfTuples( npoints );
  
    for(vtkIdType pntIdx=0; pntIdx < npoints; ++pntIdx )
    {
        double af_vals[3];//TODO 3->dims[3]
        for(int i=0; i < dims[3]; ++i)
        {
            af_vals[i]=host_a[pntIdx+i*npoints];
            value_coords[i]->SetTuple(0, &af_vals[i] );
        }
            value_grid->SetXCoordinates( value_coords[0] );
            value_grid->SetYCoordinates( value_coords[1] );
            value_grid->SetZCoordinates( value_coords[2] );
  
            vtk_xyz->SetTuple(pntIdx, value_grid->GetPoint(0) );
    } // END for all points
    grid->GetCellData()->AddArray( vtk_xyz );
    std::cout << "writer: GetNumberOfCells" << grid->GetNumberOfCells() << std::endl;
    std::cout << "writer: GetNumberOfPoints" << grid->GetNumberOfPoints() << std::endl;
    //grid->GetPointData()->AddArray( vtk_xyz );
  
    //vtkNew<vtkPointDataToCellData> pd2cd;
    //pd2cd->PassPointDataOn();
    //pd2cd->SetInputDataObject(grid);
    //pd2cd->Update();
    //grid->ShallowCopy(pd2cd->GetOutputDataObject(0));
    //
    vtkXMLRectilinearGridWriter* writer = vtkXMLRectilinearGridWriter::New();
    writer->SetFileName((outputname.append(".vtr")).c_str());
    if (verbose) std::cout<<"vtr_writer: Writing vtkRectilinearGrid Data with "<< field.dims(0) * field.dims(1) * field.dims(2) << " Points in file " << outputname << std::endl;
    writer->SetInputData( grid );
    writer->Write();
  
    writer->Delete();
    vtk_xyz->Delete();
    value_coords[0]->Delete();
    value_coords[1]->Delete();
    value_coords[2]->Delete();
    delete[] host_a;
}

// vtkRectilinearGrid Reader
// https://lorensen.github.io/VTKExamples/site/Cxx/IO/ReadRectilinearGrid/
// https://public.kitware.com/pipermail/paraview/2012-July/025678.html
void vtr_reader(af::array& field, Mesh& mesh, std::vector<double> z_spacing, std::string filepath){
    vtkSmartPointer<vtkXMLRectilinearGridReader> reader = vtkSmartPointer<vtkXMLRectilinearGridReader>::New();
    reader->SetFileName(filepath.c_str());
    reader->Update();
    vtkSmartPointer<vtkRectilinearGrid> output_data = vtkSmartPointer<vtkRectilinearGrid>::New();
    output_data=reader->GetOutput();

    int* dims = output_data->GetDimensions();//equival to int[3] array
    int datadim = output_data->GetDataDimension();
    //double* spacing = output_data->GetSpacing();
    //
    std::cout << dims[0] << std::endl;
    std::cout << dims[1] << std::endl;
    std::cout << dims[2] << std::endl;
    std::cout << datadim << std::endl;
    //from cell to point dims
    //for(int i=0; i < 3; i++){
    //    dims[i]--;
    //}

    //std::cout << dims[0] << std::endl;
    //std::cout << dims[1] << std::endl;
    //std::cout << dims[2] << std::endl;
    std::cout << datadim << std::endl;
    std::cout << output_data->GetNumberOfCells() << std::endl;
    std::cout << output_data->GetNumberOfPoints() << std::endl;
    std::cout << "GetPoint" << *output_data->GetPoint(0) << std::endl;
    std::cout << "GetPoint" << *output_data->GetPoint(1) << std::endl;
    //for (int i = 0; i < output_data->GetNumberOfPoints(); i++){
    //    std::cout << i << " GetPoint: " << *output_data->GetPoint(i) << std::endl;
    //}
    
    double* xcoords = (double*) output_data->GetXCoordinates()->GetVoidPointer(0);
    double* ycoords = (double*) output_data->GetYCoordinates()->GetVoidPointer(0);
    double* zcoords = (double*) output_data->GetZCoordinates()->GetVoidPointer(0);

        //std::cout << i << " GetPoint: " << xcoordinates->GetPointData() << std::endl;
        //std::cout << i << " GetVoidPointer: " << *(double*)xcoordinates->GetVoidPointer(i) << std::endl;
    for (int i = 0; i < output_data->GetDimensions()[0]; i++){
        std::cout << i << " GetVoidPointer x: " << xcoords[i] << std::endl;
    }
    for (int i = 0; i < output_data->GetDimensions()[1]; i++){
        std::cout << i << " GetVoidPointer y: " << ycoords[i] << std::endl;
    }
    for (int i = 0; i < output_data->GetDimensions()[2]; i++){
        std::cout << i << " GetVoidPointer z: " << zcoords[i] << std::endl;
    }

    //From: https://vtk.org/gitweb?p=VTK.git;a=blob;f=Filters/Extraction/Testing/Cxx/TestExtractRectilinearGrid.cxx
    vtkDoubleArray* xyz_data = vtkArrayDownCast<vtkDoubleArray>(output_data->GetCellData()->GetArray(0));///("xyz")
    double* xyz = static_cast<double*>( xyz_data->GetVoidPointer(0));
    std::cout << "xyz" << xyz[0] << std::endl;
    std::cout << "xyz" << xyz[1] << std::endl;
    std::cout << "xyz" << xyz[2] << std::endl;
    vtkIdType npoints = output_data->GetNumberOfCells();
    std::cout << "npoints "<< npoints << std::endl;
    for( vtkIdType pntIdx=0; pntIdx < npoints; ++pntIdx ){
        double* pnt = output_data->GetPoint( pntIdx );
        //std::cout << pntIdx << " pnt=" << *pnt << std::endl;
        std::cout << pntIdx << " xyz[pntIdx+0]=" << xyz[pntIdx+0] << std::endl;
        std::cout << pntIdx << " xyz[pntIdx+1]=" << xyz[pntIdx+1] << std::endl;
        std::cout << pntIdx << " xyz[pntIdx+2]=" << xyz[pntIdx+2] << std::endl;
        //std::cout << pntIdx << " xyz=" << xyz[pntIdx + 0] << ", " << xyz[pntIdx + 1] << ", " << xyz[pntIdx + 2] << ", " << std::endl;
    }

    std::cout << "test" << std::endl;

    vtkDataArray* vtk_xyz = output_data->GetCellData()->GetArray(0);//GetCellData() yields a vtkFieldData (like) object (conculded from docu)
    std::cout << "GetNumberOfArrays = " <<  output_data->GetCellData()->GetNumberOfArrays() << std::endl;
    //vtkAbstractArray* vtk_xyz = output_data->GetCellData()->GetAbstractArray(0);//GetCellData() yields a vtkFieldData like object (conculded from docu)

    //vtkDoubleArray* a = vtkDoubleArray::New();
    //output_data->GetCellData()->AddArray(a);

    //vtkDataArray* test = output_data->GetCellData()->GetScalars();
    //std::cout << *test->GetRange(0) << std::endl;
    std::cout << *vtk_xyz->GetRange(0) << std::endl;
    std::cout << *vtk_xyz->GetRange(1) << std::endl;
    std::cout << *vtk_xyz->GetRange(2) << std::endl;
    std::cout << vtk_xyz->GetActualMemorySize() << std::endl;
    std::cout << *vtk_xyz->GetRange() << std::endl;
    //vtkDataArray* vtk_xyz;// = vtkDataArray::NewInstance();
    //std::cout << "test" << std::endl;
    //vtk_xyz->NewInstance();
    //std::cout << "test" << std::endl;
    //vtk_xyz = output_data->GetCellData()->GetArray(0);
    //std::cout << "test" << std::endl;

    double* A_host = NULL;
    A_host = new double[datadim * output_data->GetNumberOfPoints()];
    for(int i=0; i < datadim * output_data->GetNumberOfCells(); i++){
        A_host[i] = xyz[i];
    }
    af::array A(datadim * output_data->GetNumberOfCells(), 1, 1, 1, A_host);
    delete [] A_host;
    A=af::moddims(A,af::dim4(datadim, dims[0]-1, dims[1]-1, dims[2]-1));
    A=af::reorder(A,1,2,3,0);
    af::print("A", A);


//    for (int x = 0; x < dims[0]-1; x++){
//        for (int y = 0; y < dims[1]-1; y++){
//            for (int z = 0; z < dims[2]-1; z++){
//                for (int im=0; im < 3; im++){//TODO 3-> dims4
//                    //std::cout << x << ", " << y << ", " << z << ", " << im << ", xyz=" << xyz[x+dims[0]*(y+dims[1]*(z+ dims[2] * im))] << std::endl; //<< ", " << xyz[pntIdx + 1] << ", " << xyz[pntIdx + 2] << ", " << std::endl;
//                    //double* test_xyz = static_cast<double*>( xyz_data->GetVoidPointer(x, y, z));
//                    //std::cout << x << ", " << y << ", " << z << ", " << im << ", xyz=" << test_xyz[im] << std::endl; //<< ", " << xyz[pntIdx + 1] << ", " << xyz[pntIdx + 2] << ", " << std::endl;
//            //    std::cout << i << j << k << output_data->GetCell(i, j, k) << std::endl;
//            //    std::cout << i << j << k << output_data->GetCell(i, j, k) << std::endl;
//                }
//            }
//         }
//    }
    //for (int i = 0; i < datadim * output_data->GetNumberOfPoints(); i++){
    //    std::cout << output_data->GetValue()
    //}
    //array_host = new double[datadim * output_data->GetNumberOfCells()];
    //delete [] array_host;

    //std::cout << dims[3] << std::endl;
    //std::cout << dims[4] << std::endl;
    //std::cout << dims[5] << std::endl;
    //std::cout << *spacing << std::endl;
}

////https://www.vtk.org/gitweb?p=VTK.git;a=blob;f=Examples/DataManipulation/Cxx/Arrays.cxx
////USEAGE:
////array Aout = array();
////Mesh meshout = Mesh(NaN,NaN,NaN,NaN,NaN,NaN);
////vti_to_af(Aout,meshout,"/home/pth/git/magnum.af/Data/Testing/minit.vti");
////print("Aout",Aout);
////std::cout << "  " << meshout.n0 << "  " << meshout.n1 << "  " << meshout.n2 << "  " << meshout.dx << "  " << meshout.dy << "  " << meshout.dz <<  std::endl;
//void vti_reader_micro(af::array& field, Mesh& mesh, std::string filepath){
//    int dim4th = 3;//This is the number of the components of the 3D Field (until now only 3)
//    //Candidate is imageData->GetScalarSize()<<std::endl;
//    vtkSmartPointer<vtkXMLImageDataReader> reader = vtkSmartPointer<vtkXMLImageDataReader>::New();
//    reader->SetFileName(filepath.c_str());
//    reader->Update(); 
//    vtkSmartPointer<vtkImageData> imageData = vtkSmartPointer<vtkImageData>::New();
//    imageData=reader->GetOutput();
//
//    vtkSmartPointer<vtkDoubleArray> temp = vtkSmartPointer<vtkDoubleArray>::New();
//    imageData->GetCellData()->GetScalars()->GetData(0,imageData->GetNumberOfCells()-1,0,dim4th-1,temp);
//
//    //Check weather cell or point data
//    double* test_pixel = static_cast<double*>(imageData->GetScalarPointer(0,0,0));
//    std::cout<<"testpixel="<<test_pixel<<std::endl;
//    assert(test_pixel==NULL);
//    std::cout<<"Assuming input is vtkCellData as testpixel is NULL"<<std::endl;
//    if(test_pixel==0) std::cout<<"testpixel is 0"<<std::endl;
//    //std::cout<<"testpixel="<<test_pixel[0]<<"  "<<test_pixel[1]<<"  "<<test_pixel[2]<<"  "<<std::endl;
//
//    int* dims_reduced = imageData->GetDimensions();
//    double* spacing = imageData->GetSpacing();
//    for(int i=0; i < dim4th; i++){
//        dims_reduced[i]--;
//    }
//    mesh=Mesh(dims_reduced[0],dims_reduced[1],dims_reduced[2],spacing[0],spacing[1],spacing[2]);
//    double* A_host = NULL;
//    A_host = new double[dim4th*imageData->GetNumberOfCells()];
//
//    //Two Versions to sort CellData to double* A_host 
//    //Perform equally both on CPU and OpenCL
//
//    //VERSION SORT WITH ARRAYFIRE
//    for(int i=0; i < dim4th* imageData->GetNumberOfCells(); i++){
//        A_host[i]=temp->GetValue(i);
//    }
//    af::array A(dim4th*imageData->GetNumberOfCells(),1,1,1,A_host);
//    delete [] A_host;
//    A=af::moddims(A,af::dim4(dim4th,dims_reduced[0],dims_reduced[1],dims_reduced[2]));
//    A=af::reorder(A,1,2,3,0);
//    field=A;
//
//    ////VERSION SORT WITH C++ LOOP
//    //int A_host_idx=0;
//    //for(int i=0; i < dim4th; i++){
//    //    for (int j = i; j < dim4th * imageData->GetNumberOfCells(); j = j + dim4th){
//    //        A_host[A_host_idx]=temp->GetValue(j);
//    //        A_host_idx++;
//    //    }
//    //}
//    //af::array A(dim4th*imageData->GetNumberOfCells(),1,1,1,A_host);
//    //delete [] A_host;
//    //A=af::moddims(A,af::dim4(dims_reduced[0],dims_reduced[1],dims_reduced[2],dim4th));
//    //field=A;
//}
//
//void vti_reader_atom(af::array& field, Mesh& mesh, std::string filepath){
//    int dim4th = 3;//This is the number of the components of the 3D Field (until now only 3)
//    vtkSmartPointer<vtkXMLImageDataReader> reader = vtkSmartPointer<vtkXMLImageDataReader>::New();
//    reader->SetFileName(filepath.c_str());
//    reader->Update(); 
//    vtkSmartPointer<vtkImageData> imageData = vtkSmartPointer<vtkImageData>::New();
//    imageData=reader->GetOutput();
//    std::cout<<imageData->GetNumberOfCells()<<std::endl;
//    std::cout<<imageData->GetNumberOfPoints()<<std::endl;
//    //std::cout<<reader->GetOutput();
//
//    int* dims = imageData->GetDimensions();
//    double* spacing = imageData->GetSpacing();
//    mesh=Mesh(dims[0],dims[1],dims[2],spacing[0],spacing[1],spacing[2]);
//    double* A_host = NULL;
//    A_host = new double[dim4th*imageData->GetNumberOfPoints()];
//    for (int z = 0; z < dims[2]; z++)
//      {
//      for (int y = 0; y < dims[1]; y++)
//        {
//        for (int x = 0; x < dims[0]; x++)
//          {
//          for (int im=0; im < dim4th; im++)
//            {
//            double* pixel = static_cast<double*>(imageData->GetScalarPointer(x,y,z));
//            A_host[x+dims[0]*(y+dims[1]*(z+ dims[2] * im))] = pixel[im] ;
//            }
//          }
//        }
//      }
//
//    af::array A(dim4th*imageData->GetNumberOfPoints(),1,1,1,A_host);
//    delete [] A_host;
//    A=af::moddims(A,af::dim4(dims[0],dims[1],dims[2],dim4th));
//    field=A;
//}
