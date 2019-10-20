import arrayfire as af
from libcpp.memory cimport shared_ptr
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp cimport bool

cdef extern from "../../src/llg_terms/LLGTerm.hpp" namespace "magnumaf":
    cdef cppclass LLGTerm

cdef extern from "<arrayfire.h>" namespace "af":
    cdef cppclass array:
        array()

cdef extern from "../../src/llg_terms/micro_exch.hpp" namespace "magnumaf":
    cdef cppclass ExchangeField:
        ExchangeField (long int A_field_ptr);
        ExchangeField (float A);
        long int h_ptr(const State& state);
        float E(const State& state);
        float get_cpu_time();

cdef extern from "../../src/llg_terms/micro_exch_sparse.hpp" namespace "magnumaf":
    cdef cppclass SparseExchangeField:
        SparseExchangeField (long int A_exchange_field_ptr, Mesh mesh, bool verbose);
        SparseExchangeField (float A_exchange, Mesh mesh, bool verbose);
        long int h_ptr(const State& state);
        float E(const State& state);
        float get_cpu_time();

cdef extern from "../../src/llg_terms/micro_exch_nonequi.hpp" namespace "magnumaf":
    cdef cppclass NonequiExchangeField:
        NonequiExchangeField (long int A_exchange_field_ptr, NonequispacedMesh mesh, bool verbose);
        NonequiExchangeField (float A_exchange, NonequispacedMesh mesh, bool verbose);
        long int h_ptr(const State& state);
        float E(const State& state);
        float get_cpu_time();

cdef extern from "../../src/mesh.hpp" namespace "magnumaf":
    cdef cppclass Mesh:
        int n0,n1,n2;
        float dx,dy,dz;
        int n0_exp, n1_exp, n2_exp;
        Mesh (int, int, int, float, float, float)

cdef extern from "../../src/nonequispaced_mesh.hpp" namespace "magnumaf":
    cdef cppclass NonequispacedMesh:
        int nx, ny, nz;
        float dx, dy;
        vector[float] z_spacing;
        NonequispacedMesh (int, int, float, float, vector[float] z_spacing);

cdef extern from "../../src/state.hpp" namespace "magnumaf":
    cdef cppclass State:
        State ()
        State (Mesh mesh_in, float Ms, long int m_in, bool verbose, bool mute_warning);
        State (Mesh mesh_in, long int Ms_field_ptr, long int m_in, bool verbose, bool mute_warning);
        State (Mesh mesh_in, Material param_in, long int m_in);
        State (Mesh mesh_in, Material param_in, long int aptr, long int evaluate_mean_ptr);
        void Normalize();

        array m;
        void set_m(long int aptr);
        long int get_m_addr();

        void set_Ms_field(long int aptr);
        long int get_Ms_field();

        float t;
        unsigned long long steps;
        Mesh mesh;
        NonequispacedMesh nonequimesh;
        Material material;
        float Ms;

        void write_vti(string outputname);
        void _vti_writer_micro_boolean(string outputname);
        void _vti_writer_atom (string outputname);
        void _vti_reader(string inputname);

        void _vtr_writer(string outputname);
        void _vtr_reader(string inputname);
        float meani(const int i);

cdef extern from "../../src/material.hpp" namespace "magnumaf":
    cdef cppclass Material:
        Material();
        Material(float D, float D_axis_x, float D_axis_y, float D_axis_z, float p, float J_atom, float D_atom, float K_atom, float D_atom_axis_x , float D_atom_axis_y, float D_atom_axis_z, float K_atom_axis_x, float K_atom_axis_y, float K_atom_axis_z, bool hexagonal_close_packed);
        float D;
        float D_axis[3];

        float p;
        float J_atom;
        float D_atom;
        float K_atom;
        float D_atom_axis[3];
        float K_atom_axis[3];
        bool    hexagonal_close_packed;

cdef extern from "../../src/llg_terms/atomistic_dmi.hpp" namespace "magnumaf":
    cdef cppclass AtomisticDmiField:
        AtomisticDmiField(Mesh, Material);
        long int h_ptr(const State& state);
        float E(const State& state);
        float get_cpu_time();

cdef extern from "../../src/integrators/controller.hpp" namespace "magnumaf":
    cdef cppclass Controller:
        Controller(float hmin, float hmax, float atol, float rtol);

cdef extern from "../../src/integrators/new_llg.hpp" namespace "magnumaf":
    cdef cppclass LLGIntegrator:
        LLGIntegrator (float alpha, vector[shared_ptr[LLGTerm]] vector_in, string mode, Controller);
        vector[shared_ptr[LLGTerm]] llgterms;
        void step(State& state);
        float E(const State& state);
        void relax(State& state, float precision, const int iloop, const int iwritecout);
        long int h_addr(const State& state);
        float alpha;

cdef extern from "../../src/llg_terms/micro_demag.hpp" namespace "magnumaf":
    cdef cppclass DemagField:
        DemagField (Mesh mesh_in, bool verbose, bool caching, unsigned nthreads);
        long int h_ptr(const State& state);
        float E(const State& state);
        void print_Nfft();
        float get_cpu_time();

cdef extern from "../../src/llg_terms/micro_anisotropy.hpp" namespace "magnumaf":
    cdef cppclass UniaxialAnisotropyField:
        UniaxialAnisotropyField (long int Ku1_field, float Ku1_axis_0, float Ku1_axis_1, float Ku1_axis_2);
        UniaxialAnisotropyField (float Ku1, float Ku1_axis_0, float Ku1_axis_1, float Ku1_axis_2);
        UniaxialAnisotropyField (long int Ku1_field_ptr, long int Ku1_axis_field_ptr);
        UniaxialAnisotropyField (float Ku1, long int Ku1_axis_field_ptr);
        float E(const State& state);
        long int h_ptr(const State& state);
        float Ku1;
        float get_ku1_axis(int i);
        long int get_Ku1_field();
        float get_cpu_time();

cdef extern from "../../src/llg_terms/micro_anisotropy_nonequi.hpp" namespace "magnumaf":
    cdef cppclass NonequiUniaxialAnisotropyField:
        NonequiUniaxialAnisotropyField (long int Ku1_field, float Ku1_axis_0, float Ku1_axis_1, float Ku1_axis_2);
        NonequiUniaxialAnisotropyField (float Ku1, float Ku1_axis_0, float Ku1_axis_1, float Ku1_axis_2);
        NonequiUniaxialAnisotropyField (long int Ku1_field_ptr, long int Ku1_axis_field_ptr);
        NonequiUniaxialAnisotropyField (float Ku1, long int Ku1_axis_field_ptr);
        float E(const State& state);
        long int h_ptr(const State& state);
        float Ku1;
        float get_ku1_axis(int i);
        long int get_Ku1_field();
        float get_cpu_time();

cdef extern from "../../src/llg_terms/atomistic_demag.hpp" namespace "magnumaf":
    cdef cppclass AtomisticDipoleDipoleField:
        AtomisticDipoleDipoleField(Mesh);
        long int h_ptr(const State& state);
        float E(const State& state);
        float get_cpu_time();

cdef extern from "../../src/llg_terms/atomistic_anisotropy.hpp" namespace "magnumaf":
    cdef cppclass AtomisticUniaxialAnisotropyField:
        AtomisticUniaxialAnisotropyField(Mesh, Material);
        long int h_ptr(const State& state);
        float E(const State& state);
        float get_cpu_time();

cdef extern from "../../src/llg_terms/atomistic_exchange.hpp" namespace "magnumaf":
    cdef cppclass AtomisticExchangeField:
        AtomisticExchangeField(Mesh);
        long int h_ptr(const State& state);
        float E(const State& state);
        float get_cpu_time();

cdef extern from "../../src/llg_terms/zee.hpp" namespace "magnumaf":
    cdef cppclass ExternalField:
        ExternalField (long int m_in);
        long int get_m_addr();
        long int h_ptr(const State& state);
        float E(const State& state);
        float get_cpu_time();
        void set_homogeneous_field(const float x, const float y, const float z);

cdef extern from "../../src/solvers/lbfgs_minimizer.hpp" namespace "magnumaf":
    cdef cppclass LBFGS_Minimizer:
        LBFGS_Minimizer(float tolerance_ , size_t maxIter_ , int verbose );
        vector[shared_ptr[LLGTerm]] llgterms_;
        LBFGS_Minimizer(vector[shared_ptr[LLGTerm]] vector_in, float tolerance_, size_t maxIter_, int verbose);
        float Minimize(State& state);
        float GetTimeCalcHeff();

cdef extern from "../../src/func.hpp" namespace "magnumaf":
    cdef cppclass WrappedArray:
        WrappedArray(array);
        WrappedArray(long int array_ptr);
        void set_array(long int array_ptr);
        long int get_array_addr();

cdef extern from "../../src/llg_terms/micro_spintransfertorque.hpp" namespace "magnumaf":
    cdef cppclass SpinTransferTorqueField:
        SpinTransferTorqueField (long int polarization_field_ptr, float nu_dampinglike, float nu_field, float j_e);
        long int h_ptr(const State& state);
        float E(const State& state);
        float get_cpu_time();
        WrappedArray polarization_field;

cdef extern from "../../src/vtk_IO.hpp" namespace "magnumaf":
    void pywrap_vti_writer_micro(const long int afarray_ptr, const float dx, const float dy, const float dz, string outputname);
