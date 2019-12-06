import arrayfire as af
from libcpp.memory cimport shared_ptr
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp cimport bool

cdef extern from "../../src/llg_terms/LLGTerm.hpp" namespace "magnumafcpp":
    cdef cppclass LLGTerm

cdef extern from "<arrayfire.h>" namespace "af":
    cdef cppclass array:
        array()

cdef extern from "../../src/llg_terms/micro_exch.hpp" namespace "magnumafcpp":
    cdef cppclass ExchangeField:
        ExchangeField (long int A_field_ptr);
        ExchangeField (double A);
        long int h_ptr(const State& state);
        double E(const State& state);
        double get_cpu_time();

cdef extern from "../../src/llg_terms/micro_exch_sparse.hpp" namespace "magnumafcpp":
    cdef cppclass SparseExchangeField:
        SparseExchangeField (long int A_exchange_field_ptr, Mesh mesh, bool verbose);
        SparseExchangeField (double A_exchange, Mesh mesh, bool verbose);
        long int h_ptr(const State& state);
        double E(const State& state);
        double get_cpu_time();

cdef extern from "../../src/llg_terms/micro_exch_nonequi.hpp" namespace "magnumafcpp":
    cdef cppclass NonequiExchangeField:
        NonequiExchangeField (long int A_exchange_field_ptr, NonequispacedMesh mesh, bool verbose);
        NonequiExchangeField (double A_exchange, NonequispacedMesh mesh, bool verbose);
        long int h_ptr(const State& state);
        double E(const State& state);
        double get_cpu_time();

cdef extern from "../../src/mesh.hpp" namespace "magnumafcpp":
    cdef cppclass Mesh:
        unsigned long n0,n1,n2;
        double dx,dy,dz;
        unsigned long n0_exp, n1_exp, n2_exp;
        Mesh (unsigned long, unsigned long, unsigned long, double, double, double)

cdef extern from "../../src/nonequispaced_mesh.hpp" namespace "magnumafcpp":
    cdef cppclass NonequispacedMesh:
        int nx, ny, nz;
        double dx, dy;
        vector[double] z_spacing;
        NonequispacedMesh (int, int, double, double, vector[double] z_spacing);

cdef extern from "../../src/state.hpp" namespace "magnumafcpp":
    cdef cppclass State:
        State ()
        State (Mesh mesh_in, double Ms, long int m_in, bool verbose, bool mute_warning);
        State (Mesh mesh_in, long int Ms_field_ptr, long int m_in, bool verbose, bool mute_warning);
        State (Mesh mesh_in, Material param_in, long int m_in);
        State (Mesh mesh_in, Material param_in, long int aptr, long int evaluate_mean_ptr);
        void Normalize();

        array m;
        void set_m(long int aptr);
        long int get_m_addr();

        void set_Ms_field(long int aptr);
        long int get_Ms_field();

        double t;
        unsigned long long steps;
        Mesh mesh;
        NonequispacedMesh nonequimesh;
        Material material;
        double Ms;

        void write_vti(string outputname);
        void _vti_writer_micro_boolean(string outputname);
        void _vti_writer_atom (string outputname);
        void _vti_reader(string inputname);

        void _vtr_writer(string outputname);
        void _vtr_reader(string inputname);
        double meani(const int i);

cdef extern from "../../src/material.hpp" namespace "magnumafcpp":
    cdef cppclass Material:
        Material();
        Material(double D, double D_axis_x, double D_axis_y, double D_axis_z, double p, double J_atom, double D_atom, double K_atom, double D_atom_axis_x , double D_atom_axis_y, double D_atom_axis_z, double K_atom_axis_x, double K_atom_axis_y, double K_atom_axis_z, bool hexagonal_close_packed);
        double D;
        double D_axis[3];

        double p;
        double J_atom;
        double D_atom;
        double K_atom;
        double D_atom_axis[3];
        double K_atom_axis[3];
        bool    hexagonal_close_packed;

cdef extern from "../../src/llg_terms/atomistic_dmi.hpp" namespace "magnumafcpp":
    cdef cppclass AtomisticDmiField:
        AtomisticDmiField(Mesh, Material);
        long int h_ptr(const State& state);
        double E(const State& state);
        double get_cpu_time();

cdef extern from "../../src/integrators/controller.hpp" namespace "magnumafcpp":
    cdef cppclass Controller:
        Controller(double hmin, double hmax, double atol, double rtol);

cdef extern from "../../src/integrators/new_llg.hpp" namespace "magnumafcpp":
    cdef cppclass LLGIntegrator:
        LLGIntegrator (double alpha, vector[shared_ptr[LLGTerm]] vector_in, string mode, Controller);
        vector[shared_ptr[LLGTerm]] llgterms;
        void step(State& state);
        double E(const State& state);
        void relax(State& state, double precision, const int iloop, const int iwritecout);
        long int h_addr(const State& state);
        double alpha;

cdef extern from "../../src/llg_terms/micro_demag.hpp" namespace "magnumafcpp":
    cdef cppclass DemagField:
        DemagField (Mesh mesh_in, bool verbose, bool caching, unsigned nthreads);
        long int h_ptr(const State& state);
        double E(const State& state);
        void print_Nfft();
        double get_cpu_time();

cdef extern from "../../src/llg_terms/micro_anisotropy.hpp" namespace "magnumafcpp":
    cdef cppclass UniaxialAnisotropyField:
        UniaxialAnisotropyField (long int Ku1_field, double Ku1_axis_0, double Ku1_axis_1, double Ku1_axis_2);
        UniaxialAnisotropyField (double Ku1, double Ku1_axis_0, double Ku1_axis_1, double Ku1_axis_2);
        UniaxialAnisotropyField (long int Ku1_field_ptr, long int Ku1_axis_field_ptr);
        UniaxialAnisotropyField (double Ku1, long int Ku1_axis_field_ptr);
        double E(const State& state);
        long int h_ptr(const State& state);
        double Ku1;
        double get_ku1_axis(int i);
        long int get_Ku1_field();
        double get_cpu_time();

cdef extern from "../../src/llg_terms/micro_anisotropy_nonequi.hpp" namespace "magnumafcpp":
    cdef cppclass NonequiUniaxialAnisotropyField:
        NonequiUniaxialAnisotropyField (long int Ku1_field, double Ku1_axis_0, double Ku1_axis_1, double Ku1_axis_2);
        NonequiUniaxialAnisotropyField (double Ku1, double Ku1_axis_0, double Ku1_axis_1, double Ku1_axis_2);
        NonequiUniaxialAnisotropyField (long int Ku1_field_ptr, long int Ku1_axis_field_ptr);
        NonequiUniaxialAnisotropyField (double Ku1, long int Ku1_axis_field_ptr);
        double E(const State& state);
        long int h_ptr(const State& state);
        double Ku1;
        double get_ku1_axis(int i);
        long int get_Ku1_field();
        double get_cpu_time();

cdef extern from "../../src/llg_terms/atomistic_demag.hpp" namespace "magnumafcpp":
    cdef cppclass AtomisticDipoleDipoleField:
        AtomisticDipoleDipoleField(Mesh);
        long int h_ptr(const State& state);
        double E(const State& state);
        double get_cpu_time();

cdef extern from "../../src/llg_terms/atomistic_anisotropy.hpp" namespace "magnumafcpp":
    cdef cppclass AtomisticUniaxialAnisotropyField:
        AtomisticUniaxialAnisotropyField(Mesh, Material);
        long int h_ptr(const State& state);
        double E(const State& state);
        double get_cpu_time();

cdef extern from "../../src/llg_terms/atomistic_exchange.hpp" namespace "magnumafcpp":
    cdef cppclass AtomisticExchangeField:
        AtomisticExchangeField(Mesh);
        long int h_ptr(const State& state);
        double E(const State& state);
        double get_cpu_time();

cdef extern from "../../src/llg_terms/zee.hpp" namespace "magnumafcpp":
    cdef cppclass ExternalField:
        ExternalField (long int m_in);
        long int get_m_addr();
        long int h_ptr(const State& state);
        double E(const State& state);
        double get_cpu_time();
        void set_homogeneous_field(const double x, const double y, const double z);

cdef extern from "../../src/solvers/lbfgs_minimizer.hpp" namespace "magnumafcpp":
    cdef cppclass LBFGS_Minimizer:
        LBFGS_Minimizer(double tolerance_ , size_t maxIter_ , int verbose );
        vector[shared_ptr[LLGTerm]] llgterms_;
        LBFGS_Minimizer(vector[shared_ptr[LLGTerm]] vector_in, double tolerance_, size_t maxIter_, int verbose);
        double Minimize(State& state);
        double GetTimeCalcHeff();

cdef extern from "../../src/func.hpp" namespace "magnumafcpp":
    cdef cppclass WrappedArray:
        WrappedArray(array);
        WrappedArray(long int array_ptr);
        void set_array(long int array_ptr);
        long int get_array_addr();

cdef extern from "../../src/llg_terms/micro_spintransfertorque.hpp" namespace "magnumafcpp":
    cdef cppclass SpinTransferTorqueField:
        SpinTransferTorqueField (long int polarization_field_ptr, double nu_dampinglike, double nu_field, double j_e);
        long int h_ptr(const State& state);
        double E(const State& state);
        double get_cpu_time();
        WrappedArray polarization_field;

cdef extern from "../../src/llg_terms/micro_exch_rkky.hpp" namespace "magnumafcpp":
    cdef cppclass RKKYExchangeField:
        RKKYExchangeField (long int rkky_values, long int exchange_values, Mesh mesh, long int rkky_indices, bool verbose);
        long int h_ptr(const State& state);
        double E(const State& state);

cdef extern from "../../src/vtk_IO.hpp" namespace "magnumafcpp":
    void pywrap_vti_writer_micro(const long int afarray_ptr, const double dx, const double dy, const double dz, string outputname);
