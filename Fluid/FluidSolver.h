#include "FluidQuantity.h"

namespace Nova
{
template <class T, int d>
class FluidSolver
{
    using TV = Vector<T, d>;
    using T_INDEX = Vector<int, d>;

    // FluidQuantity<T, d> *_d;
    // FluidQuantity<T, d> *_u;
    // FluidQuantity<T, d> *_v;

    // Density Field
    FluidQuantity<T, d> *_d;
    // Velocity Field
    FluidQuantity<T, d> *_v[d];

    FluidSimulator_Grid<T, d> *grid;

    int m;
    int n;

    T hx;
    T rho;

    T *_rhs;
    T *_p;
    T *_z;      /* Auxiliary vector */
    T *_s;      /* Search vector */
    T *_precon; /* Preconditioner */

    T *_Adiag;  /* Matrix diagonal */
    T *_Aplusi; /* Matrix off-diagonals */
    T *_Aplusj;

    void calculateRHS()
    {
        memset(_rhs, 0, m * n * sizeof(T));
        for (int iy = 0; iy < n; iy++)
        {
            for (int ix = 0; ix < m; ix++)
            {
                int idx = iy * m + ix;
                // T_INDEX index = offset2index(idx);

                // for(int axis = 0; axis < d; axis++)
                // {
                //     _rhs[idx] -= (_v[axis]->at(index) - _v[axis]->at(grid->Next_Cell(axis, index))) / hx;
                // }

                _rhs[idx] -= (_v[0]->at(ix + 1, iy) - _v[0]->at(ix, iy)) / hx;
                _rhs[idx] -= (_v[1]->at(ix, iy + 1) - _v[1]->at(ix, iy)) / hx;
            }
        }
    }

    /* Builds the pressure matrix. Since the matrix is very sparse and
     * symmetric, it allows for memory friendly storage.
     */
    void calculateA(T timestep)
    {
        T scale = timestep / rho / hx / hx;
        memset(_Adiag, 0, m * n * sizeof(T));

        for (int iy = 0; iy < n; iy++)
        {
            for (int ix = 0; ix < m; ix++)
            {
                int idx = iy * m + ix;
                if (ix < m - 1) // if u.next.valid()
                {
                    _Adiag[idx] += scale;     // A.diag() + 1
                    _Adiag[idx + 1] += scale; // A.nextU.diag() + 1
                    _Aplusi[idx] = -scale;
                }
                else
                {
                    _Aplusi[idx] = 0;
                }

                if (iy < n - 1) //if v.next.valid()
                {
                    _Adiag[idx] += scale;
                    _Adiag[idx + m] += scale;
                    _Aplusj[idx] = -scale;
                }
                else
                {
                    _Aplusj[idx] = 0;
                }
            }
        }
    }

    void applyPreconditioner(T *dst, T *src)
    {
        // Trying to solve M * z = r
        // Option 1: Set z = r ====> M = I (Do nothing)
        memcpy(dst, src, m * n * sizeof(T));

        // Option 2: M ~ A^T
        // MIC / IC / integrade
    }

    // x * y
    T InnerProduct(T *a, T *b)
    {
        T result = 0.0;
        for (int i = 0; i < m * n; i++)
            result += a[i] * b[i];
        return result;
    }

    // y <- A * x
    // A is a sparse matrix! stored in _Adiag, _Aplusi, _Aplusj
    // Also we only care about inertia area
    // A_ij,ij related to Adiag(ij), Aplusi(ij), Aplusj(ij), Aplusi(previousi), Aplusj(previousj)
    void Multiply(T *_dst, T *_src)
    {
        for (int iy = 0; iy < n; iy++)
        {
            for (int ix = 0; ix < m; ix++)
            {
                int idx = iy * m + ix;
                T sum = _Adiag[idx] * _src[idx];

                if (ix > 0)
                { // if u.previous valid
                    sum += _Aplusi[idx - 1] * _src[idx - 1];
                }

                if (iy > 0)
                {
                    sum += _Aplusj[idx - m] * _src[idx - m];
                }

                if (ix < m - 1)
                {
                    sum += _Aplusi[idx] * _src[idx + 1];
                }

                if (iy < n - 1)
                {
                    sum += _Aplusj[idx] * _src[idx + m];
                }

                _dst[idx] = sum;
            }
        }
    }

    // dst = x + y * k
    void saxpy(T *_dst, T *_x, T *_y, T k)
    {
        for (int i = 0; i < m * n; i++)
            _dst[i] = _x[i] + _y[i] * k;
    }

    // |x|_infinity
    // In fact should also cares about !inertior case
    T ConvergenceNorm(T *_x)
    {
        T maxNorm = __DBL_MIN__;
        for (int i = 0; i < m * n; i++)
        {
            maxNorm = std::max(maxNorm, fabs(_x[i]));
            // maxNorm = max(fabs(_x[i]), maxNorm);
        }
        return maxNorm;
    }

    void project_GS(int limit, T timestep, bool output = true)
    {
        T scale = timestep / rho / hx / hx;

        T maxDelta;

        for (int iteration = 0; iteration < limit; iteration++)
        {
            maxDelta = __DBL_MIN__;
            for (int iy = 0; iy < n; iy++)
            {
                for (int ix = 0; ix < m; ix++)
                {
                    int idx = iy * m + ix;
                    T Aii = 0, sum = 0;

                    if (ix > 0)
                    {
                        Aii += scale;
                        sum -= scale * _p[idx - 1]; // Previous u
                    }
                    if (iy > 0)
                    {
                        Aii += scale;
                        sum -= scale * _p[idx - m]; // Previous v
                    }
                    if (ix < m - 1)
                    {
                        Aii += scale;
                        sum -= scale * _p[idx + 1]; // Next u
                    }
                    if (iy < n - 1)
                    {
                        Aii += scale;
                        sum -= scale * _p[idx + m]; // Next v
                    }

                    T newP = (_rhs[idx] - sum) / Aii;

                    maxDelta = std::max(maxDelta, fabs(_p[idx] - newP));

                    _p[idx] = newP;
                }
            }
            if (maxDelta < 1e-5)
            {
                if (output)
                    printf("Converge with %d iteration, with Norm1 = %f\n", iteration, maxDelta);

                return;
            }
        }

        if (output)
            printf("Exceed Limit of %d, with Norm1 = %f\n", limit, maxDelta);
    }

    void project_CG(int limit, bool output = true)
    {
        // Set initial guess p = 0 and r = d(the RHS)
        memset(_p, 0, m * n * sizeof(T));
        // Set auxiliary vector z = applyPreconditioner(r)
        applyPreconditioner(_z, _rhs);

        // and search vector s = z;     dst, src, size
        memcpy(_s, _z, m * n * sizeof(T));

        T sigma = InnerProduct(_z, _rhs);
        // printf("sigma = %f\n", sigma);

        T norm1 = ConvergenceNorm(_rhs);
        if (norm1 < 1e-5)
            return;

        // Enter Iteration with limit
        for (int iteration = 0; iteration < limit; iteration++)
        {
            // Set Auxiliary vector z : Multiply(z, s)
            Multiply(_z, _s);

            T alpha = sigma / InnerProduct(_z, _s);
            // printf("alpha = %f\n", sigma);

            // Update p and r
            saxpy(_p, _p, _s, alpha);
            saxpy(_rhs, _rhs, _z, -alpha);

            // If with convergence-norm < tol, return
            norm1 = ConvergenceNorm(_rhs);
            if (norm1 < 1e-5)
            {
                if (output)
                    printf("Converge when iteration = %d, with 1Norm = %f\n", iteration, norm1);
                return;
            }

            // Set auxiliary vector z = applyPreconditioner(r)
            applyPreconditioner(_z, _rhs);
            // Sigma_new = dotproduct(z, r)
            T sigma_new = InnerProduct(_z, _rhs);
            // printf("sigma_new = %f\n", sigma_new);
            T beta = sigma_new / sigma;
            // Set search vector s = z + beta * s
            // saxpy(_s, beta, _s, _z);

            saxpy(_s, _z, _s, beta);
            sigma = sigma_new;
        }
        if (output)
            printf("The %d iteration Limit Exceeded with 1Norm = %f\n", limit, norm1);
    }

    void setBoundaryCondition()
    {
        for (int x = 0; x < m; x++)
        {
            _v[1]->at(x, 0) = 0.0;
            _v[1]->at(x, n) = 0.0;
        }

        for (int y = 0; y < n; y++)
        {
            _v[0]->at(m, y) = 0.0;
            _v[0]->at(0, y) = 0.0;
        }
    }

    void applyPressure(T timestep)
    {
        T scale = timestep / rho / hx;
        for (int iy = 0; iy < n; iy++)
        {
            for (int ix = 0; ix < m; ix++)
            {
                int idx = iy * m + ix;

                _v[0]->at(ix, iy) -= _p[idx] * scale;
                _v[0]->at(ix + 1, iy) += _p[idx] * scale;
                _v[1]->at(ix, iy) -= _p[idx] * scale;
                _v[1]->at(ix, iy + 1) += _p[idx] * scale;
            }
        }

        setBoundaryCondition();
    }

    int index2offset(const T_INDEX &index) const
    {
        return grid->index2offset(index, grid->counts);
    }

    T_INDEX offset2index(const int os) const
    {
        return grid->offset2index(os, grid->counts);
    }

  public:
    FluidSolver(FluidSimulator_Grid<T, d> &grid, T rho) : grid(&grid), rho(rho)
    {
        this->m = grid.counts[0];
        this->n = grid.counts[1];

        this->hx = grid.hx;

        _d = new FluidQuantity<T, d>(grid, -1);

        for (int axis = 0; axis < d; axis++)
            _v[axis] = new FluidQuantity<T, d>(grid, axis);

        _rhs = new T[m * n];
        _p = new T[m * n];
        _z = new T[m * n];
        _s = new T[m * n];
        _Adiag = new T[m * n];
        _Aplusi = new T[m * n];
        _Aplusj = new T[m * n];
        _precon = new T[m * n];
    }

    ~FluidSolver()
    {
        delete _d;
        for (int axis = 0; axis < d; axis++)
           delete  _v[axis];

        delete[] _rhs;
        delete[] _p;
        delete[] _z;
        delete[] _s;
        delete[] _Adiag;
        delete[] _Aplusi;
        delete[] _Aplusj;
        delete[] _precon;
    }

    void update(T timestep)
    {
        // Projection
        calculateRHS();
        //calculateA(timestep);
        //project_CG(1000);
        project_GS(1000, timestep);
        applyPressure(timestep);

        //Advection
        _d->advect(timestep, _v);
        _v[0]->advect(timestep, _v);
        _v[1]->advect(timestep, _v);

        // Flip
        _d->flip();
        _v[0]->flip();
        _v[1]->flip();
    }

    void addInflow(int ix0, int iy0, int ix1, int iy1, int axis, T value)
    {
        for (int y = std::max(iy0, 0); y < std::min(iy1, n); y++)
            for (int x = std::max(ix0, 0); x < std::min(ix1, m); x++)
            {
                int idx = y * m + x;
                //printf("idx = %d\n", idx);
                if (axis == -1)
                    _d->addInflow(x, y, value);
                if (axis == 0)
                    _v[0]->addInflow(x, y, value);
                if (axis == 1)
                    _v[1]->addInflow(x, y, value);
            }
    }

    /* Convert fluid density to RGB color scaled in 0 - 1 */
    T toRGB(int x, int y)
    {
        int idx = y * m + x;
        return std::max(std::min(1.0 - _d->Phi()[idx], 1.0), 0.0);
    }

    //===========================================================

    // Projection output as a tent-like plot
    void test_projection_tent(T timestep)
    {
        // Projection
        // calculateRHS();
        memset(_rhs, 0, m * n * sizeof(T));
        _rhs[7740] = 1;
        // calculateA(timestep);
        // calculatePreconditioner();
        project_GS(1000, timestep, false);
        // project_CG(1000, false);
        applyPressure(timestep);
    }
};
} // namespace Nova