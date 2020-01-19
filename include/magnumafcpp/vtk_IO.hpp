#pragma once
#include "mesh.hpp"
#include "nonequispaced_mesh.hpp"
#include "arrayfire.h"

namespace magnumafcpp
{

void vti_writer_micro(const af::array field, const Mesh &mesh, std::string outputname); //3D Image data
void pywrap_vti_writer_micro(const long int afarray_ptr, const double dx, const double dy, const double dz, const std::string outputname);
void vti_writer_atom(const af::array field, const Mesh &mesh, std::string outputname); //3D Image data

void vti_reader(af::array &field, Mesh &mesh, std::string filepath);
//void vti_reader_micro(af::array& field, Mesh& mesh, std::string filepath);
//void vti_reader_atom(af::array& field, Mesh& mesh, std::string filepath);

//void vtr_writer(const af::array field, const Mesh& mesh, std::string outputname);//Rectilinear grid writer
void vtr_writer(const af::array &field, const double dx, const double dy, const std::vector<double> z_spacing, std::string outputname, const bool verbose = false); //Rectilinear grid writer
void vtr_writer(const af::array &field, const NonequispacedMesh &nonequimesh, std::string outputname, const bool verbose = false);                                  //Rectilinear grid writer
void vtr_reader(af::array &field, NonequispacedMesh &mesh, std::string filepath, const bool verbose = false);
} // namespace magnumafcpp
