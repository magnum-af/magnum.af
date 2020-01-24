#include <boost/python.hpp>
#include "func.hpp"


BOOST_PYTHON_MODULE(PYBINDTESTLIB)
{
    Py_Initialize();
    using namespace boost::python;
    def("greet", greet);
}
