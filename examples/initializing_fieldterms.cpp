#include "magnum_af.hpp"
#include <cassert>

// Demonstrating initialization of fieldterms
// Mainly with helper function templates fieldterm::to_vec and fieldterm::mv_to_vec

using namespace magnumafcpp;
int main() {
    // Parameter initialization
    const double x = 5.e-7, y = 1.25e-7, z = 3.e-9;
    const int nx = 100, ny = 25, nz = 1;
    const double A = 1.3e-11, Ms = 8e5, alpha = 0.2;

    // Generating Objects
    Mesh mesh(nx, ny, nz, x / nx, y / ny, z / nz);

    // Initial magnetic field
    af::array m = af::constant(0, nx, ny, nz, 3, f64);
    m(0, af::span, af::span, 1) = 1;
    m(af::seq(1, af::end - 1), af::span, af::span, 0) = 1;
    m(-1, af::span, af::span, 1) = 1;

    // State object
    State state(mesh, Ms, m);

    DemagField dmag(mesh, false, true, 0);
    ExchangeField exch(A);

    // copy fieldterms to LLGIntegrator, retaining copy
    LLGIntegrator llg1(alpha, fieldterm::to_vec(dmag, exch));
    // we can dmag,exch fieldterms afterwards:
    exch.h(state);
    dmag.h(state);

    // when we use std::move, fieldterm pointers are moved-from
    LLGIntegrator llg2(alpha, fieldterm::to_vec(std::move(dmag), std::move(exch)));
    // We can not use dmag, exch anylonger; would cause segfault:
    // exch.h(state); // segfaults
    // dmag.h(state); // segfaults

    // with fieldterm::mv_to_vec, elements are implicitly moved from, use with care
    DemagField dmag2(mesh, false, true, 0);
    ExchangeField exch2(A);
    LLGIntegrator llg3(alpha, fieldterm::mv_to_vec(dmag2, exch2));
    // dmag2.h(state); // segfaults
    // exch2.h(state); // segfaults

    // when fieldterms are const, nothing is moved
    {
        const DemagField dmag(mesh, false, true, 0);
        const ExchangeField exch(A);
        LLGIntegrator llg1(
            alpha, fieldterm::to_vec(std::move(dmag), std::move(exch))); // does not move as dmag, exch are const
        LLGIntegrator llg2(alpha, fieldterm::mv_to_vec(dmag, exch));     // does not move as dmag, exch are const
        // we can still use them
        exch.h(state);
        dmag.h(state);
    }

    // Various ways to construct LLGIntegrator containing a copy of fieldterms
    {
        DemagField dmag(mesh, false, true, 0);
        ExchangeField exch(A);

        // prefered way, using fieldterm::to_vec helper template function
        LLGIntegrator llg1(alpha, fieldterm::to_vec(dmag, exch));

        // copy the fieldterms into a vector-of-pointer using cp_to_uptr helper
        LLGIntegrator llg2(alpha, {magnumafcpp::fieldterm::cp_to_uptr(dmag), magnumafcpp::fieldterm::cp_to_uptr(exch)});

        // state explicitly what's going on
        LLGIntegrator llg3(alpha, {std::unique_ptr<FieldTerm>(std::make_unique<DemagField>(dmag)),
                                   std::unique_ptr<FieldTerm>(std::make_unique<ExchangeField>(exch))});

        // versions using new (should be avoided by using make_unique)
        LLGIntegrator llg4(alpha, {uptr_FieldTerm(new DemagField(dmag)), uptr_FieldTerm(new ExchangeField(exch))});
        LLGIntegrator llg5(alpha, {std::unique_ptr<FieldTerm>(new DemagField(dmag)),
                                   std::unique_ptr<FieldTerm>(new ExchangeField(exch))});

        // we can still use the fieldterms
        exch.h(state);
        dmag.h(state);
    }

    // Example creating uptr in advance and then moving from them
    {
        auto dmag = fieldterm::to_uptr<DemagField>(mesh, false, true, 0);
        auto exch = fieldterm::to_uptr<ExchangeField>(A);
        LLGIntegrator llg(alpha, {std::move(dmag), std::move(exch)});
        // after std::move, the unique_ptr is set to void and should not be used anymore
        // using it leads to a segfault
        assert(dmag == nullptr);
        assert(exch == nullptr);
    }

    // creating to_uptr in place:
    {
        LLGIntegrator llg(alpha);
        llg.llgterms.push_back(std::make_unique<DemagField>(mesh, false, true, 0));
        llg.llgterms.push_back(std::make_unique<ExchangeField>(A));
    }

    return 0;
}
