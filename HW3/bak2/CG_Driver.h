
//!#####################################################################
//! \file CG_Driver.h
//!#####################################################################
// Class CG_Driver
//######################################################################
#ifndef __CG_Driver__
#define __CG_Driver__

#include "CG_System.h"
#include <nova/Tools/Grids/Grid.h>
#include <nova/Tools/Krylov_Solvers/Conjugate_Gradient.h>

namespace Nova
{
template <class T, int d>
class CG_Driver
{
  private:
    using TV = Vector<T, d>;

    Grid<T, d> *grid;
    // CG_Storage<T,d> *storage;

    int size, m, n;

    double hx;

    // cg preference
    T cg_tolerance;
    int cg_min_iterations;
    int cg_max_iterations;

  public:
    CG_Driver(Grid<T, d> &grid)
        : grid(&grid)
    {
        size = grid.counts.Product();
        hx = grid.one_over_dX.Max();
    }
    ~CG_Driver()
    {
    }

    void SetCGPreference(const T tolerance, const int min_iterations, const int max_iterations)
    {
        cg_tolerance = tolerance;
        cg_min_iterations = min_iterations;
        cg_max_iterations = max_iterations;
    }

    void Execute()
    {

        Conjugate_Gradient<T> cg;

        cg.print_residuals = true;
        cg.print_diagnostics = true;
        cg.restart_iterations = 80;

        Array<TV> delta_X(size);
        Array<TV> rhs(size);

        Array<TV> temp_q(size), temp_s(size), temp_r(size), temp_k(size), temp_z(size);
        CG_Vector<T, d> cg_x(delta_X), cg_b(rhs), cg_q(temp_q), cg_s(temp_s), cg_r(temp_r), cg_k(temp_k), cg_z(temp_z);

        CG_System<T, d> cg_system();

        SetCGPreference(1e-5, 0, 500);

        // solve
        cg.Solve(cg_system, cg_x, cg_b, cg_q, cg_s, cg_r, cg_k, cg_z, cg_tolerance, cg_min_iterations, cg_max_iterations);
    }
};
} // namespace Nova

#endif
