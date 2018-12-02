
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <iostream>

using namespace std;

class FluidQuantity
{
    double *_Phi;
    double *_Phi_new;

    // m * n grid
    int m;
    int n;
    // Face indicator
    int axis;
    double ox, oy;
    // Grid Cell size
    double hx;

    // 1D Linear Interpolate between a and b in (0,1)
    double linp(double a, double b, double delta) const
    {
        return a * (1 - delta) + b * delta;
    }

  public:
    FluidQuantity(int w, int h, int axis, double hx)
        : m(w), n(h), axis(axis), hx(hx)
    {
        ox = (axis == 0) ? 0 : 0.5;
        oy = (axis == 1) ? 0 : 0.5;

        _Phi = new double[m * n];
        _Phi_new = new double[m * n];

        memset(_Phi, 0, m * n * sizeof(double));
    }

    ~FluidQuantity()
    {
        delete[] _Phi;
        delete[] _Phi_new;
    }

    void flip()
    {
        swap(_Phi, _Phi_new);
    }

    const double *Phi() const
    {
        return _Phi;
    }

    double at(int x, int y) const
    {
        return _Phi[x + y * m];
    }

    double &at(int x, int y)
    {
        return _Phi[x + y * m];
    }

    // Linear Interpolate on grid at index (x, y) (can be 0.5 if on face)
    double linp(double x, double y) const
    {
        // Clamp Coordinates
        // Extra 0.0001 for not Clamp to m+1 when x is exactly a int.
        x = max(x - ox, 0.0);
        x = min(x, m - 1.0001);
        y = max(y - oy, 0.0);
        y = min(y, n - 1.0001);

        // Get Index
        int ix = (int)x;
        int iy = (int)y;
        // Get the offset between (0,1)
        x -= ix;
        y -= iy;

        double x0 = linp(at(ix, iy), at(ix + 1, iy), x);
        double x1 = linp(at(ix, iy + 1), at(ix + 1, iy + 1), x);
        double re = linp(x0, x1, y);

        return re;
    }

    void advect(double timestep, const FluidQuantity &u, const FluidQuantity &v)
    {
        // Reduce Spatial locality

        for (int iy = 0; iy < n; iy++)
        {
            for (int ix = 0; ix < m; ix++)
            {
                int idx = iy * m + ix;

                double x = ix + ox;
                double y = iy + oy;

                // Compute velocity
                double velocity_u = u.linp(x, y) / hx;
                double velocity_v = v.linp(x, y) / hx;

                // Traceback
                x -= velocity_u * timestep;
                y -= velocity_v * timestep;

                _Phi_new[idx] = linp(x, y);
            }
        }
    }

    void addInflow(int ix, int iy, double value)
    {
        // if(ix >= 0 & ix < m & iy >=0 & iy < n)
        if (fabs(_Phi[iy * m + ix]) < fabs(value))
        {
            _Phi[iy * m + ix] = value;
            printf("Add Inflow at (%d, %d)\n", ix, iy);
        }
    }

    // /* Sets fluid quantity inside the given rect to value `v' */
    // void addInflow(double x0, double y0, double x1, double y1, double v)
    // {
    //     int ix0 = (int)(x0 / hx - ox);
    //     int iy0 = (int)(y0 / hx - oy);
    //     int ix1 = (int)(x1 / hx - ox);
    //     int iy1 = (int)(y1 / hx - oy);

    //     for (int y = max(iy0, 0); y < min(iy1, n); y++)
    //         for (int x = max(ix0, 0); x < min(ix1, n); x++)
    //             if (fabs(_Phi[x + y * m]) < fabs(v))
    //                 _Phi[x + y * m] = v;
    // }
};

class FluidSolver
{
    FluidQuantity *_d;
    FluidQuantity *_u;
    FluidQuantity *_v;

    int m;
    int n;

    double hx;
    double rho;

    double *_rhs;
    double *_p;
    double *_z;      /* Auxiliary vector */
    double *_s;      /* Search vector */
    double *_precon; /* Preconditioner */

    double *_Adiag;  /* Matrix diagonal */
    double *_Aplusi; /* Matrix off-diagonals */
    double *_Aplusj;

    void calculateRHS()
    {
        memset(_rhs, 0, m * n * sizeof(double));
        for (int iy = 0; iy < n; iy++)
        {
            for (int ix = 0; ix < m; ix++)
            {
                int idx = iy * m + ix;
                _rhs[idx] -= (_u->at(ix + 1, iy) - _u->at(ix, iy)) / hx;
                _rhs[idx] -= (_v->at(ix, iy + 1) - _v->at(ix, iy)) / hx;
            }
        }
    }

    /* Builds the pressure matrix. Since the matrix is very sparse and
     * symmetric, it allows for memory friendly storage.
     */
    void calculateA(double timestep)
    {
        double scale = timestep / rho / hx / hx;
        memset(_Adiag, 0, m * n * sizeof(double));

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

    // /* Builds the modified incomplete Cholesky preconditioner */
    // void calculatePreconditioner()
    // {
    //     const double tau = 0.97;
    //     const double sigma = 0.25;

    //     for (int y = 0, idx = 0; y < n; y++)
    //     {
    //         for (int x = 0; x < m; x++, idx++)
    //         {
    //             double e = _Adiag[idx];

    //             if (x > 0)
    //             {
    //                 double px = _Aplusi[idx - 1] * _precon[idx - 1];
    //                 double py = _Aplusj[idx - 1] * _precon[idx - 1];
    //                 e = e - (px * px + tau * px * py);
    //             }
    //             if (y > 0)
    //             {
    //                 double px = _Aplusi[idx - m] * _precon[idx - m];
    //                 double py = _Aplusj[idx - m] * _precon[idx - m];
    //                 e = e - (py * py + tau * px * py);
    //             }

    //             if (e < sigma * _Adiag[idx])
    //                 e = _Adiag[idx];

    //             _precon[idx] = 1.0 / sqrt(e);
    //         }
    //     }
    // }

    void applyPreconditioner(double *dst, double *src)
    {
        // Trying to solve M * z = r
        // Option 1: Set z = r ====> M = I (Do nothing)
        memcpy(dst, src, m * n * sizeof(double));

        // Option 2: M ~ A^T
        // MIC / IC / integrade
    }

    // x * y
    double InnerProduct(double *a, double *b)
    {
        double result = 0.0;
        for (int i = 0; i < m * n; i++)
            result += a[i] * b[i];
        return result;
    }

    // y <- A * x
    // A is a sparse matrix! stored in _Adiag, _Aplusi, _Aplusj
    // Also we only care about inertia area
    // A_ij,ij related to Adiag(ij), Aplusi(ij), Aplusj(ij), Aplusi(previousi), Aplusj(previousj)
    void Multiply(double *_dst, double *_src)
    {
        for (int iy = 0; iy < n; iy++)
        {
            for (int ix = 0; ix < m; ix++)
            {
                int idx = iy * m + ix;
                double sum = _Adiag[idx] * _src[idx];

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
    void saxpy(double *_dst, double *_x, double *_y, double k)
    {
        for (int i = 0; i < m * n; i++)
            _dst[i] = _x[i] + _y[i] * k;
    }

    // |x|_infinity
    // In fact should also cares about !inertior case
    double ConvergenceNorm(double *_x)
    {
        double maxNorm = __DBL_MIN__;
        for (int i = 0; i < m * n; i++)
        {
            maxNorm = max(maxNorm, fabs(_x[i]));
            // maxNorm = max(fabs(_x[i]), maxNorm);
        }
        return maxNorm;
    }

    void project_GS(int limit, double timestep, bool output = true)
    {
        double scale = timestep / rho / hx / hx;

        double maxDelta;

        for (int iteration = 0; iteration < limit; iteration++)
        {
            maxDelta = __DBL_MIN__;
            for (int iy = 0; iy < n; iy++)
            {
                for (int ix = 0; ix < m; ix++)
                {
                    int idx = iy * m + ix;
                    double Aii = 0, sum = 0;

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

                    double newP = (_rhs[idx] - sum) / Aii;

                    maxDelta = max(maxDelta, fabs(_p[idx] - newP));

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
        memset(_p, 0, m * n * sizeof(double));
        // Set auxiliary vector z = applyPreconditioner(r)
        applyPreconditioner(_z, _rhs);

        // printArray(_rhs);

        // printArray(_z);

        // and search vector s = z;     dst, src, size
        memcpy(_s, _z, m * n * sizeof(double));

        double sigma = InnerProduct(_z, _rhs);
        // printf("sigma = %f\n", sigma);

        double norm1 = ConvergenceNorm(_rhs);
        if (norm1 < 1e-5)
            return;

        // Enter Iteration with limit
        for (int iteration = 0; iteration < limit; iteration++)
        {
            // Set Auxiliary vector z : Multiply(z, s)
            Multiply(_z, _s);

            double alpha = sigma / InnerProduct(_z, _s);
            // printf("alpha = %f\n", sigma);

            // Update p and r
            // saxpy(_p, alpha, _s, _p);      // p += alpha * s
            // saxpy(_rhs, -alpha, _z, _rhs); // r -= alpha * z
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
            double sigma_new = InnerProduct(_z, _rhs);
            // printf("sigma_new = %f\n", sigma_new);
            double beta = sigma_new / sigma;
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
            _v->at(x, 0) = 0.0;
            _v->at(x, n) = 0.0;
        }

        for (int y = 0; y < n; y++)
        {
            _u->at(m, y) = 0.0;
            _u->at(0, y) = 0.0;
        }
    }

    void applyPressure(double timestep)
    {
        double scale = timestep / rho / hx;
        for (int iy = 0; iy < n; iy++)
        {
            for (int ix = 0; ix < m; ix++)
            {
                int idx = iy * m + ix;

                _u->at(ix, iy) -= _p[idx] * scale;
                _u->at(ix + 1, iy) += _p[idx] * scale;
                _v->at(ix, iy) -= _p[idx] * scale;
                _v->at(ix, iy + 1) += _p[idx] * scale;
            }
        }

        setBoundaryCondition();
    }

  public:
    FluidSolver(int w, int h, double density) : m(w), n(h), rho(density)
    {
        hx = 1.0 / min(w, h);

        // _d = new FluidQuantity(m, n, 0.5, 0.5, hx);
        // _u = new FluidQuantity(m + 1, n, 0.0, 0.5, hx);
        // _v = new FluidQuantity(m, n + 1, 0.5, 0.0, hx);

        _d = new FluidQuantity(m, n, -1, hx);
        _u = new FluidQuantity(m + 1, n, 0, hx);
        _v = new FluidQuantity(m, n + 1, 1, hx);

        _rhs = new double[m * n];
        _p = new double[m * n];
        _z = new double[m * n];
        _s = new double[m * n];
        _Adiag = new double[m * n];
        _Aplusi = new double[m * n];
        _Aplusj = new double[m * n];
        _precon = new double[m * n];
    }

    ~FluidSolver()
    {
        delete _d;
        delete _u;
        delete _v;

        delete[] _rhs;
        delete[] _p;
        delete[] _z;
        delete[] _s;
        delete[] _Adiag;
        delete[] _Aplusi;
        delete[] _Aplusj;
        delete[] _precon;
    }

    void update(double timestep)
    {
        // Projection
        calculateRHS();
        calculateA(timestep);
        project_CG(1000);
        // project_GS(1000, timestep);
        applyPressure(timestep);

        //Advection
        _d->advect(timestep, *_u, *_v);
        _u->advect(timestep, *_u, *_v);
        _v->advect(timestep, *_u, *_v);

        // Flip
        _d->flip();
        _u->flip();
        _v->flip();
    }

    void addInflow(int ix0, int iy0, int ix1, int iy1, int axis, double value)
    {
        for (int y = max(iy0, 0); y < min(iy1, n); y++)
            for (int x = max(ix0, 0); x < min(ix1, m); x++)
            {
                int idx = y * m + x;
                if (axis == -1)
                    _d->addInflow(x, y, value);
                if (axis == 0)
                    _u->addInflow(x, y, value);
                if (axis == 1)
                    _v->addInflow(x, y, value);
            }
    }

    // void addInflow(double x, double y, double w, double h, double d, double u, double v)
    // {
    //     _d->addInflow(x, y, x + w, y + h, d);
    //     _u->addInflow(x, y, x + w, y + h, u);
    //     _v->addInflow(x, y, x + w, y + h, v);
    // }

    /* Convert fluid density to RGB color scaled in 0 - 1 */
    double toRGB(int x, int y)
    {
        int idx = y * m + x;
        return max(min(1.0 - _d->Phi()[idx], 1.0), 0.0);
    }

    void printArray(double *src)
    {
        for (int i = 0; i < m * n; i++)
        {
            printf("%.3f", src[i]);

            if ((i + 1) % m == 0)
            {
                printf("\n");
            }
            else
            {
                printf(",");
            }
        }
    }

    //===========================================================

    // Projection output as a tent-like plot
    void test_projection_tent(double timestep)
    {
        // Projection
        calculateRHS();
        calculateA(timestep);
        // calculatePreconditioner();
        // project_GS(1000, timestep, false);
        project_CG(1000, false);
        applyPressure(timestep);

        printArray(_p);
    }
};
