import arrayfire as af 
from libcpp.memory cimport shared_ptr
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp cimport bool

cdef extern from "../../src/llg_terms/LLGTerm.hpp":
  cdef cppclass LLGTerm

cdef extern from "<arrayfire.h>" namespace "af":
  cdef cppclass array:
    array()

cdef extern from "../../src/llg_terms/micro_exch.hpp":
  cdef cppclass ExchangeField:
    ExchangeField (Mesh mesh_in, Material param_in);
    double E(const State& state);
    double get_cpu_time();

cdef extern from "../../src/mesh.hpp":
  cdef cppclass Mesh:
    int n0,n1,n2;
    double dx,dy,dz;
    int n0_exp, n1_exp, n2_exp;
    Mesh (int, int, int, double, double, double)

cdef extern from "../../src/state.hpp":
  cdef cppclass State:
    State ()
    State (Mesh mesh_in, Material param_in, long int m_in);
    State (Mesh mesh_in, Material param_in, long int aptr, long int evaluate_mean_ptr);
    void Normalize();

    array m;
    void set_m(long int aptr);
    long int get_m_addr();
    array Ms;
    void set_micro_Ms_field(long int aptr);
    long int get_micro_Ms_field();
    array micro_A_field;
    void set_micro_A_field(long int aptr);
    long int get_micro_A_field();
    array micro_Ku1_field;
    void set_micro_Ku1_field(long int aptr);
    long int get_micro_Ku1_field();

    double t;
    Mesh mesh;
    Material material;

    void _vti_writer_micro(string outputname);
    void _vti_writer_micro_boolean(string outputname);
    void _vti_writer_atom (string outputname);
    void _vti_reader(string inputname);

    void _vtr_writer(string outputname);
    void _vtr_reader(string inputname);
    double meani(const int i);

cdef extern from "../../src/material.hpp":
  cdef cppclass Material:
    Material();
    Material(double alpha, double T, double ms, double A, double D, double Ku1, double D_axis_x, double D_axis_y, double D_axis_z, double Ku1_axis_x, double Ku1_axis_y, double Ku1_axis_z, double p, double J_atom, double D_atom, double K_atom, double D_atom_axis_x , double D_atom_axis_y, double D_atom_axis_z, double K_atom_axis_x, double K_atom_axis_y, double K_atom_axis_z, bool hexagonal_close_packed, int mode, bool afsync);
    double mu0;
    double gamma;
    double alpha;
    double T;

    double ms;
    double A;
    double D;
    double Ku1;
    double D_axis[3];
    double Ku1_axis[3];

    double p;
    double J_atom;
    double D_atom;
    double K_atom;
    double D_atom_axis[3];
    double K_atom_axis[3];
    bool  hexagonal_close_packed;

    int mode;
    bool afsync;

cdef extern from "../../src/llg_terms/atomistic_dmi.hpp":
  cdef cppclass AtomisticDmiField:
    AtomisticDmiField(Mesh, Material);
    double E(const State& state);
    double get_cpu_time();

cdef extern from "../../src/integrators/controller.hpp":
  cdef cppclass Controller:
    Controller(double hmin, double hmax, double atol, double rtol);

cdef extern from "../../src/integrators/new_llg.hpp":
  cdef cppclass LLGIntegrator:
    LLGIntegrator (vector[shared_ptr[LLGTerm]] vector_in, string mode, Controller);
    vector[shared_ptr[LLGTerm]] llgterms;
    void step(State& state);
    double E(const State& state);
    long int get_fheff_addr(const State& state);
cdef extern from "../../src/llg_terms/micro_demag.hpp":
  cdef cppclass DemagField:
    DemagField (Mesh mesh_in, Material param_in, bool verbose, bool caching);
    double E(const State& state);
    void print_Nfft();
    double get_cpu_time();

cdef extern from "../../src/llg_terms/micro_anisotropy.hpp":
  cdef cppclass UniaxialAnisotropyField:
    UniaxialAnisotropyField(Mesh, Material);
    double E(const State& state);
    double get_cpu_time();

cdef extern from "../../src/llg_terms/atomistic_demag.hpp":
  cdef cppclass AtomisticDipoleDipoleField:
    AtomisticDipoleDipoleField(Mesh);
    double E(const State& state);
    double get_cpu_time();

cdef extern from "../../src/llg_terms/atomistic_anisotropy.hpp":
  cdef cppclass AtomisticUniaxialAnisotropyField:
    AtomisticUniaxialAnisotropyField(Mesh, Material);
    double E(const State& state);
    double get_cpu_time();

cdef extern from "../../src/llg_terms/atomistic_exchange.hpp":
  cdef cppclass AtomisticExchangeField:
    AtomisticExchangeField(Mesh);
    double E(const State& state);
    double get_cpu_time();

cdef extern from "../../src/llg_terms/zee.hpp":
  cdef cppclass Zee:
    Zee (long int m_in);
    long int get_m_addr();
    double E(const State& state);
    double get_cpu_time();
    void set_xyz(const State&, const double x, const double y, const double z);

cdef extern from "../../src/solvers/lbfgs_minimizer.hpp":
  cdef cppclass LBFGS_Minimizer:
    LBFGS_Minimizer(double tolerance_ , size_t maxIter_ , int verbose_ );
    vector[shared_ptr[LLGTerm]] llgterms_;
    LBFGS_Minimizer(vector[shared_ptr[LLGTerm]] vector_in, double tolerance_, size_t maxIter_, int verbose_);
    double Minimize(State& state);
    double GetTimeCalcHeff();

cdef extern from "../../src/func.hpp":
  cdef cppclass WrappedArray:
    WrappedArray(array);
    WrappedArray(long int array_ptr);
    void set_array(long int array_ptr);
    long int get_array_addr();

cdef extern from "../../src/llg_terms/micro_spintransfertorque.hpp":
  cdef cppclass SpinTransferTorqueField:
    SpinTransferTorqueField (long int polarization_field_ptr, double nu_dampinglike, double nu_field, double j_e);
    double E(const State& state);
    double get_cpu_time();
    WrappedArray polarization_field;
