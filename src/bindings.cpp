// following: https://vsamy.github.io/en/blog/boost-python-cmake-build
#include "func.hpp"
//#include "magnum_af.hpp"
#include <boost/python.hpp>
// Include the headers of MyLib

BOOST_PYTHON_MODULE(pyMyLib)
{
    //Py_Initialize();

    boost::python::def("greet", greet);
    // Write the bindings here
}
