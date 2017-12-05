#ifndef VTK_WRITER_H
#define VTK_WRITER_H
#include<string>

#include "arrayfire.h"

#include "mesh.hpp"

//af_to_vti
#include <vtkSmartPointer.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkImageData.h>
//#include <vtkVersion.h>
//#include <vtkProperty.h>
//#include <vtkXMLImageDataReader.h>
//#include <vtkRenderWindowInteractor.h>
//#include <vtkRenderer.h>
//#include <vtkImageDataGeometryFilter.h>
//#include <vtkRenderWindow.h>
//#include <vtkActor.h>
//#include <vtkPolyDataMapper.h>

//af_to_vtk
#include "vtkCellData.h"
#include "vtkDoubleArray.h"
#include "vtkExtractRectilinearGrid.h"
#include "vtkPointData.h"
#include "vtkRectilinearGridWriter.h"
//#include "vtkDataArray.h"
//#include "vtkMathUtilities.h"
//#include "vtkNew.h"
//#include "vtkPointDataToCellData.h"
//#include "vtkRectilinearGrid.h"
//#include "vtkStructuredData.h"
//
//vti_to_af
#include <vtkXMLImageDataReader.h>
#include <vtkDataSet.h>

#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkProperty.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkXMLImageDataReader.h>
#include <vtkImageData.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkImageDataGeometryFilter.h>
#include <vtkRenderWindow.h>
#include <vtkActor.h>

#include <vtkDoubleArray.h>

void af_to_vtk(const af::array field, const Mesh& mesh, std::string outputname);//Rectilinear grid writer
void af_to_vti(const af::array field, const Mesh& mesh, std::string outputname);//3D Image data

void vti_to_af(const af::array& field, const Mesh& mesh, std::string filepath);
#endif
