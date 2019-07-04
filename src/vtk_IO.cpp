#include "vtk_IO.hpp"
//af_to_vti
#include <vtkSmartPointer.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkImageData.h>
//af_to_vtk
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkExtractRectilinearGrid.h>
#include <vtkPointData.h>
//vti_to_af
#include <vtkXMLImageDataReader.h>
#include <vtkDataSet.h>
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <vtkXMLRectilinearGridWriter.h>
#include <vtkXMLRectilinearGridReader.h>
#include <assert.h>
#include <string>

namespace magnumaf{



//3D vtkImageData vtkCellData writer
//Optimization should avoid generation of two vtkImageData objects
void implementation_vti_writer_micro(const af::array field, const Mesh& mesh, std::string outputname){

    try{
    double* host_a = field.host<double>();

    vtkSmartPointer<vtkImageData> imageDataPointCentered = vtkSmartPointer<vtkImageData>::New();
    imageDataPointCentered->SetDimensions(field.dims(0), field.dims(1), field.dims(2));
    imageDataPointCentered->SetSpacing(mesh.dx, mesh.dy, mesh.dz);
    imageDataPointCentered->SetOrigin(0, 0, 0);
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
            double* pixel = static_cast<double*>(imageDataPointCentered->GetScalarPointer(x, y, z));
            pixel[im] = host_a[x+dims[0]*(y+dims[1]*(z+ dims[2] * im))];
            }
          }
        }
      }

    vtkSmartPointer<vtkImageData> imageDataCellCentered = vtkSmartPointer<vtkImageData>::New();
    imageDataCellCentered->SetDimensions(field.dims(0)+1, field.dims(1)+1, field.dims(2)+1);
    imageDataCellCentered->SetOrigin(0, 0, 0);
    imageDataCellCentered->SetSpacing(mesh.dx, mesh.dy, mesh.dz);
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
    catch(const std::exception& e){
        printf("\33[1;31mWarning:\33[0m exception in vti_writer_micro:\n%s\n", e.what());
    }
}

void vti_writer_micro(const af::array field, const Mesh& mesh, std::string outputname){
    implementation_vti_writer_micro(field, mesh, outputname);
}

void pywrap_vti_writer_micro(const long int afarray_ptr, const double dx, const double dy, const double dz, const std::string outputname){
    af::array afarray = *(new af::array( *((void **) afarray_ptr)));
    Mesh mesh(afarray.dims(0), afarray.dims(1), afarray.dims(2), dx, dy, dz);
    implementation_vti_writer_micro(afarray, mesh, outputname);
}

//3D vtkImageData vtkPointData writer
void vti_writer_atom(const af::array field, const Mesh& mesh, std::string outputname){

    double* host_a = field.host<double>();

    vtkSmartPointer<vtkImageData> imageData = vtkSmartPointer<vtkImageData>::New();
    imageData->SetDimensions(field.dims(0), field.dims(1), field.dims(2));
    imageData->SetSpacing(mesh.dx, mesh.dy, mesh.dz);
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
            double* pixel = static_cast<double*>(imageData->GetScalarPointer(x, y, z));
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
    double* test_pixel = static_cast<double*>(imageData->GetScalarPointer(0, 0, 0));
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
        imageData->GetCellData()->GetScalars()->GetData(0, imageData->GetNumberOfCells()-1, 0, dim4th-1, temp);
        double* A_host = NULL;
        A_host = new double[dim4th*imageData->GetNumberOfCells()];
        //VERSION SORT WITH ARRAYFIRE
        for(int i=0; i < dim4th* imageData->GetNumberOfCells(); i++){
            A_host[i]=temp->GetValue(i);
        }
        af::array A(dim4th*imageData->GetNumberOfCells(), 1, 1, 1, A_host);
        delete [] A_host;
        A=af::moddims(A, af::dim4(dim4th, dims[0], dims[1], dims[2]));
        A=af::reorder(A, 1, 2, 3, 0);
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
        //af::array A(dim4th*imageData->GetNumberOfCells(), 1, 1, 1, A_host);
        //delete [] A_host;
        //A=af::moddims(A, af::dim4(dims_reduced[0], dims_reduced[1], dims_reduced[2], dim4th));
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
              double* pixel = static_cast<double*>(imageData->GetScalarPointer(x, y, z));
              for (int im=0; im < dim4th; im++)
                {
                A_host[x+dims[0]*(y+dims[1]*(z+ dims[2] * im))] = pixel[im] ;
                }
              }
            }
          }

        af::array A(dim4th*imageData->GetNumberOfPoints(), 1, 1, 1, A_host);
        delete [] A_host;
        A=af::moddims(A, af::dim4(dims[0], dims[1], dims[2], dim4th));
        field=A;
    }
    mesh=Mesh(dims[0], dims[1], dims[2], spacing[0], spacing[1], spacing[2]);
}


void vtr_writer(const af::array field, const Mesh& mesh, const std::vector<double> z_spacing, std::string outputname, const bool verbose){
    // writes af::array field as vtkRectilinearGrid to file outputname.vtr
    // Uses cell data, thus xyz dimensions are increased by 1:  field.dims(0)+1, field.dims(1)+1, field.dims(2)+1
    // supports arbitrary field.dims(3) (e.g. 3 for vector field, 1 for scalar field)

    // Creating vtk Object and setting dimensions obtained from af::array
    vtkSmartPointer<vtkRectilinearGrid> grid = vtkSmartPointer<vtkRectilinearGrid>::New();
    grid->SetDimensions(field.dims(0)+1, field.dims(1)+1, field.dims(2)+1);// Adding one node per dimension as we use cell data


    // declare xyz coordinate vectors
    vtkDataArray* coords[3];
    for(int i=0; i < 3; ++i){
        coords[i] = vtkDataArray::CreateDataArray(VTK_DOUBLE);
        coords[i]->SetNumberOfTuples(field.dims(i)+1);
    }

    // Calculate and populate coordinate vectors
    //x
    for(int j=0; j < field.dims(0)+1; ++j){
        double val = (double)j * mesh.dx;
        coords[0]->SetTuple(j, &val);
    }

    //y
    for(int j=0; j < field.dims(1)+1; ++j){
        double val = (double)j * mesh.dy;
        coords[1]->SetTuple(j, &val);
    }

    //z
    double add_val = 0;
    for(int j=0; j < field.dims(2)+1; ++j)
    {
        if ( j == 0){
            add_val = 0;
        }
        else{
          add_val += z_spacing.at(j-1);
        }
        coords[2]->SetTuple(j, &add_val);
    }

    grid->SetXCoordinates( coords[0] );
    grid->SetYCoordinates( coords[1] );
    grid->SetZCoordinates( coords[2] );

    coords[0]->Delete();
    coords[1]->Delete();
    coords[2]->Delete();

    // Write af::array data as vtk cell data
    const double* host_a = field.host<double>();// accesses af::array raw data
    const vtkIdType ncells = grid->GetNumberOfCells();

    vtkDoubleArray* data = vtkDoubleArray::New();
    data->SetName("m");
    data->SetNumberOfComponents(field.dims(3));
    data->SetNumberOfTuples( ncells );

    for(vtkIdType cellId=0; cellId < ncells; ++cellId){
        for(int i=0; i < field.dims(3); ++i){
            data->SetValue(field.dims(3) * cellId + i, host_a[cellId+i*ncells]);
        }
    }

    grid->GetCellData()->AddArray(data);

    // Writing output
    vtkXMLRectilinearGridWriter* writer = vtkXMLRectilinearGridWriter::New();
    if (std::string(outputname.end()-4, outputname.end()) != ".vtr") outputname.append(".vtr");// adding file extensions if needed

    if (verbose) std::cout<<"vtk_writer: writing array of dimension ["<< field.dims(0) << " " << field.dims(1) << \
        " " << field.dims(2) << " " << field.dims(3)<< "] to '" << outputname << "'" << std::endl;

    writer->SetFileName(outputname.c_str());
    writer->SetInputData( grid );
    writer->Write();


    writer->Delete();
    data->Delete();
    delete[] host_a;
}


void vtr_reader(af::array& field, Mesh& mesh, std::vector<double>& z_spacing, std::string filepath, const bool verbose){
    // Counterpart to vtr_writer()
    // Reads vktRectilinearGrid cell data from file
    // https://lorensen.github.io/VTKExamples/site/Cxx/IO/ReadRectilinearGrid/
    // https://vtk.org/gitweb?p=VTK.git;a=blob;f=Filters/Extraction/Testing/Cxx/TestExtractRectilinearGrid.cxx

    if (std::string(filepath.end()-4, filepath.end()) != ".vtr") filepath.append(".vtr");// adding file extension if needed

    // Obtain vktRectilinearGrid object from reader
    vtkSmartPointer<vtkXMLRectilinearGridReader> reader = vtkSmartPointer<vtkXMLRectilinearGridReader>::New();
    reader->SetFileName(filepath.c_str());
    reader->Update();
    vtkSmartPointer<vtkRectilinearGrid> output_data = vtkSmartPointer<vtkRectilinearGrid>::New();
    output_data=reader->GetOutput();

    // Converting coordinate vectors to spacing vectors
    // E.g. double[4] = {0, 1, 2, 3} -> double[3] = {1, 1, 1}
    double* xcoords = (double*) output_data->GetXCoordinates()->GetVoidPointer(0);
    double* ycoords = (double*) output_data->GetYCoordinates()->GetVoidPointer(0);
    double* zcoords = (double*) output_data->GetZCoordinates()->GetVoidPointer(0);

    // Calculating spacings from coordinate vectors
    std::vector<double> x_spacings;
    for (int i = 0; i < output_data->GetDimensions()[0]-1; i++){
        x_spacings.push_back(xcoords[i+1] - xcoords[i]);
    }

    std::vector<double> y_spacings;
    for (int i = 0; i < output_data->GetDimensions()[1]-1; i++){
        y_spacings.push_back(ycoords[i+1] - ycoords[i]);
    }

    std::vector<double> vec_z_spacing;
    for (int i = 0; i < output_data->GetDimensions()[2]-1; i++){
        vec_z_spacing.push_back(zcoords[i+1] - zcoords[i]);
    }

    // Copying vtkCellData to af::array
    int* grid_dims = output_data->GetDimensions();// equivalent to int[3] array. Note: this accesses the raw data
    vtkDoubleArray* xyz_data = vtkArrayDownCast<vtkDoubleArray>(output_data->GetCellData()->GetArray(0));///("xyz")
    const int data_dim = xyz_data->GetNumberOfComponents();
    double* xyz = static_cast<double*>( xyz_data->GetVoidPointer(0));
    double* A_host = NULL;
    A_host = new double[data_dim * output_data->GetNumberOfCells()];

    for(int i=0; i < data_dim * output_data->GetNumberOfCells(); i++){
        A_host[i] = xyz[i];
    }

    af::array A(data_dim * output_data->GetNumberOfCells(), 1, 1, 1, A_host);
    delete [] A_host;
    A=af::moddims(A, af::dim4(data_dim, grid_dims[0]-1, grid_dims[1]-1, grid_dims[2]-1));
    A=af::reorder(A, 1, 2, 3, 0);

    // Printing dimension info
    if(verbose) std::cout << "vtr_reader: read vtkCellData of dimension [" << grid_dims[0]-1 << ", " << grid_dims[1]-1 \
        << ", " << grid_dims[2]-1 << ", " << data_dim << "] from '" << filepath << "'" << std::endl;

    // Setting output variables
    field = A;
    z_spacing = vec_z_spacing;
    mesh=Mesh(grid_dims[0]-1, grid_dims[1]-1, grid_dims[2]-1, x_spacings[0], y_spacings[1], 0);//Note: dz is set to zero
    //TODO should be adapted with nonequi Mesh class
}
}// namespace magnumaf
