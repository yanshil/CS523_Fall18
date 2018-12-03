//!#####################################################################
//! \file CG_Storage.h
//!#####################################################################
// Class CG_Storage
//######################################################################
#ifndef __CG_Storage__
#define __CG_Storage__
#include <nova/Tools/Grids/Grid.h>
#include <nova/Tools/Arrays/Array.h>

namespace Nova
{
template <class T, int d>
class CG_Storage
{
    using TV = Vector<T, d>;
    using T_INDEX = Vector<int, d>;

  public:
    T_INDEX counts;
    int size, m, n;
    int cg_iterations, cg_restart_iterations;
    T cg_tolerance;

    Array<TV> _Adiag;
    Array<TV> _Aplusi;
    Array<TV> _Aplusj;
    T one_over_dX_square;

    Array<TV> _rhs; // Right Hand Size
    int *_BF;       // Boundary Flag: -1 for Inertia and 0 for Dirichlet and 1 for Exteria

    // ==============================================================
    CG_Storage(Grid<T, d> &grid)
    {
        counts = grid.counts;
        size = counts.Product();
        // Specify convenience for 2D
        m = counts(0);
        n = counts(1);

        T hx = grid.one_over_dX.Max();
        one_over_dX_square = hx * hx;
    }

    ~CG_Storage()
    {
    }

    void setting(int cg_tolerance_input = 1e-5, int cg_iterations_input = 600, int cg_restart_iterations_input = 100)
    {
        cg_iterations = cg_iterations_input;
        cg_restart_iterations = cg_restart_iterations_input;
        cg_tolerance = cg_tolerance_input;
    }

    void initialize()
    {
        _Adiag.resize(size);
        _Aplusi.resize(size);
        _Aplusj.resize(size);
        _rhs.resize(size);
        
    }

    void calculateA()
    {
        for (int iy = 0; iy < n; iy++)
        {
            for (int ix = 0; ix < m; ix++)
            {
                int idx = iy * m + ix;
                if (ix < m - 1) // if u.next.valid() !exterior
                {
                    _Adiag(idx) += one_over_dX_square;     // _Adiag() + 1
                    _Adiag(idx + 1) += one_over_dX_square; // _AnextU.diag() + 1
                    _Aplusi(idx) = -one_over_dX_square;
                }
                else
                {
                    _Aplusi(idx) = 0;
                }

                if (iy < n - 1) //if v.next.valid()
                {
                    _Adiag(idx) += one_over_dX_square;
                    _Adiag(idx + m) += one_over_dX_square;
                    _Aplusj(idx) = -one_over_dX_square;
                }
                else
                {
                    _Aplusj(idx) = 0;
                }
            }
        }
    }

    
    void calculateRHS()
    {
        // memset(_rhs, 0, m * n * sizeof(double));
        for (int iy = 0; iy < n; iy++)
        {
            for (int ix = 0; ix < m; ix++)
            {
                int idx = iy * m + ix;
                _rhs(idx) = 0;
            }
        }

        // test
        _rhs(10) = 1;
    }
};
} // namespace Nova
#endif