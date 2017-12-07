#ifndef VTK_WRITER_H
#define VTK_WRITER_H

#include<string>
#include "arrayfire.h"
#include "mesh.hpp"
#include <assert.h>

//af_to_vti
#include <vtkSmartPointer.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkImageData.h>
//af_to_vtk
#include "vtkCellData.h"
#include "vtkDoubleArray.h"
#include "vtkExtractRectilinearGrid.h"
#include "vtkPointData.h"
#include "vtkRectilinearGridWriter.h"
//vti_to_af
#include <vtkXMLImageDataReader.h>
#include <vtkDataSet.h>
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>

void af_to_vtk(const af::array field, const Mesh& mesh, std::string outputname);//Rectilinear grid writer

void vti_writer(const af::array field, const Mesh& mesh, std::string outputname);//3D Image data
void atomistic_vti_writer(const af::array field, const Mesh& mesh, std::string outputname);//3D Image data

void vti_reader(af::array& field, Mesh& mesh, std::string filepath);
void atomistic_vti_reader(af::array& field, Mesh& mesh, std::string filepath);

#endif
