#include <boost/python.hpp>
#include "func.hpp"
#include "mesh.hpp"

using namespace boost::python;
using namespace magnumafcpp;

BOOST_PYTHON_MODULE(PYBINDTESTLIB)
{
    Py_Initialize();
    def("greet", greet);

    class_<Mesh>("Mesh", init<uint32_t, uint32_t, uint32_t, double, double, double>())
        .def(init<>())
        .def_readwrite("V", &Mesh::V)
        .def_readwrite("nx", &Mesh::n0)
        .def_readwrite("ny", &Mesh::n1)
        .def_readwrite("nz", &Mesh::n2)
        .def_readwrite("dx", &Mesh::dx)
        .def_readwrite("dy", &Mesh::dy)
        .def_readwrite("dz", &Mesh::dz)
    ;
}
