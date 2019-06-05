#ifndef VTK_WRITER_H
#define VTK_WRITER_H

#include <string>
#include "arrayfire.h"
#include "mesh.hpp"
#include <assert.h>

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

void vti_writer_micro(const af::array field, const Mesh& mesh, std::string outputname);//3D Image data
void pywrap_vti_writer_micro(const long int afarray_ptr, const double dx, const double dy, const double dz, const std::string outputname);
void vti_writer_atom(const af::array field, const Mesh& mesh, std::string outputname);//3D Image data

void vti_reader(af::array& field, Mesh& mesh, std::string filepath);
//void vti_reader_micro(af::array& field, Mesh& mesh, std::string filepath);
//void vti_reader_atom(af::array& field, Mesh& mesh, std::string filepath);

//void vtr_writer(const af::array field, const Mesh& mesh, std::string outputname);//Rectilinear grid writer
void vtr_writer(const af::array field, const Mesh& mesh, const std::vector<double> z_spacing, std::string outputname, const bool verbose = false);//Rectilinear grid writer
void vtr_reader(af::array& field, Mesh& mesh, std::vector<double>& z_spacing, std::string filepath, const bool verbose = true);
#endif
