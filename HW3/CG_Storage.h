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
    using T1 = Vector<T, 1>;
    using T_INDEX = Vector<int, d>;

  public:
    T_INDEX counts;
    int size, m, n;
    int cg_iterations, cg_restart_iterations;
    double dX, dY;
    T cg_tolerance;

    Array<T1> _Adiag;
    Array<T1> _Aplusi;
    Array<T1> _Aplusj;
    Array<T1> _trueValue;
    T one_over_dX_square;

    Array<T1> _rhs; // Right Hand Size
    int *_BF;       // Boundary Flag: -1 for Inertia and 0 for Dirichlet and 1 for Exteria

    // ==============================================================
    CG_Storage(Grid<T, d> &grid)
    {
        counts = grid.counts;
        size = counts.Product();
        // Specify convenience for 2D
        m = counts(0);
        n = counts(1);
        dX = grid.dX(0);
        dY = grid.dX(1);

        T hx = grid.one_over_dX.Max();
        one_over_dX_square = hx * hx;
    }

    ~CG_Storage()
    {
    }

    void setting(double cg_tolerance_input = 1e-5, int cg_iterations_input = 1000, int cg_restart_iterations_input = 1000)
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
        _trueValue.resize(size);
    }

    void printAdiag()
    {
        for (int i = 0; i < size; i++)
        {
            std::cout << _Adiag(i) << ", ";

            if ((i + 1) % this->m == 0)
            {
                std::cout << "\n";
            }
        }
        std::cout << std::endl;
    }

    void calculateA()
    {
        for (int iy = 0; iy < n; iy++)
        {
            for (int ix = 0; ix < m; ix++)
            {
                int idx = iy * m + ix;

                //======= That's for Exterior Problem =======================
                // if (ix < m - 1) // if u.next.valid() !exterior
                // {
                //     _Adiag(idx) += 1;     // _Adiag() + 1
                //     _Adiag(idx + 1) += 1; // _AnextU.diag() + 1
                //     _Aplusi(idx) = -1;
                // }
                // else
                // {
                //     _Aplusi(idx) = 0;
                // }

                // if (iy < n - 1) //if v.next.valid()
                // {
                //     _Adiag(idx) += 1;
                //     _Adiag(idx + m) += 1;
                //     _Aplusj(idx) = -1;
                // }
                // else
                // {
                //     _Aplusj(idx) = 0;
                // }
                //========================================

                _Adiag(idx) = 4;
                if (ix < m - 1) // if u.next.valid() !exterior
                {
                    _Aplusi(idx) = -1;
                }
                else
                {
                    _Aplusi(idx) = 0;
                }

                if (iy < n - 1) //if v.next.valid()
                {
                    _Aplusj(idx) = -1;
                }
                else
                {
                    _Aplusj(idx) = 0;
                }
                //=========================================

            }
        }
        _Adiag *= one_over_dX_square;
        _Aplusi *= one_over_dX_square;
        _Aplusj *= one_over_dX_square;
    }

    void calculateRHS()
    {
        // memset(_rhs, 0, m * n * sizeof(double));
        for (int idx = 0; idx < size; idx++)
        {
            _rhs(idx) = -4;
        }
    }

    void calculateTrueValue()
    {
        double radius_square = 0.25 * 0.25;
        for (int iy = 0; iy < n; iy++)
        {
            for (int ix = 0; ix < m; ix++)
            {
                int idx = iy * m + ix;
                // Get center of cell (ix, iy)
                double x = -0.5 + ix * dX + dX / 2.0, y = -0.5 + iy * dY + dY / 2.0;
                double phi_tmp = x * x + y * y - radius_square;

                // Inside Circle
                if (phi_tmp <= 0)
                {
                    _trueValue(idx) = phi_tmp;
                }
                else
                {
                    _trueValue(idx) = 0;
                }
            }
        }
    }

    //---------------------------------------------------
    void printAFromStorage()
    {
        for (int j = 0; j < size; j++)
        {
            for (int i = 0; i < size; i++)
            {
                int idx = j * size + i;
                double re;

                if (i == j)
                {
                    re = _Adiag(i).Max();
                }
                else if ((i + 1) == j)
                {
                    re = _Aplusi(i).Max();
                }
                else if ((i - 1) == j)
                {
                    re = _Aplusi(j).Max();
                }
                else if ((i + m) == j)
                {
                    re = _Aplusj(i).Max();
                }
                else if ((i - m) == j)
                {
                    re = _Aplusj(j).Max();
                }
                else
                {
                    re = 0;
                }

                printf("%3.0f", re);

                if ((idx + 1) % size == 0)
                {
                    std::cout << "\n";
                }
                else
                {
                    std::cout << ",";
                }
            }
        }

        std::cout << std::endl;
    }

    void testA()
    {
        _Adiag(10) = 1;
    }

    void testRHS()
    {
        for (int iy = 0; iy < n; iy++)
        {
            for (int ix = 0; ix < m; ix++)
            {
                int idx = iy * m + ix;
                _rhs(idx) = -1;
            }
        }
    }
};
} // namespace Nova
#endif