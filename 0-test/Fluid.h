
#include <math.h>
#include <stdio.h>
#include <string.h>

#define TEST_CIRCLE false
#define TEST_CIRCLE_DENSITY false
#define TEST_SINGLEINFLOW true

using namespace std;

class FluidQuantity
{
  private:
    double *_Phi;
    double *_Phi_new;

    // m * n grid
    int m;
    int n;
    // Grid Cell size
    double dX;

    // Face indicator
    int axis;
    double ox, oy;

    // 1D Linear Interpolate between a and b in (0,1)
    double linp(double a, double b, double delta) const
    {
        return a * (1 - delta) + b * delta;
    }

  public:
    FluidQuantity(int m, int n, int axis)
        : m(n), n(n), axis(axis)
    {
        ox = (axis == 0) ? 0 : 0.5;
        oy = (axis == 1) ? 0 : 0.5;
        dX = 1.0 / min(m, n);

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
        std::swap(_Phi, _Phi_new);
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

    // Linear Interpolate on grid at (x, y)
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

    double advect(double timestep, const FluidQuantity &u, const FluidQuantity &v)
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
                // TODO / dX ????
                double velocity_u = u.linp(x, y) / dX;
                double velocity_v = v.linp(x, y) / dX;

                // Traceback
                x -= velocity_u * timestep;
                y -= velocity_v * timestep;

                _Phi_new[idx] = linp(x, y);
            }
        }
    }

    // TODO
    void addInflow(double x0, double y0, double x1, double y1, double v)
    {
        int ix0 = (int)(x0 / dX - ox);
        int ix1 = (int)(x1 / dX - ox);

        int iy0 = (int)(y0 / dX - oy);
        int iy1 = (int)(y1 / dX - oy);

        for (int y = max(iy0, 0); y < min(iy1, n); y++)
            for (int x = max(ix0, 0); x < min(ix1, n); x++)
                if (fabs(_Phi[x + y * m]) < fabs(v))
                    _Phi[x + y * m] = v;
    }
};

class FluidSolver
{
  private:
    FluidQuantity *_d;
    FluidQuantity *_u;
    FluidQuantity *_v;

    int m;
    int n;
    double dX;
    double rho; // rho

    double *_rhs;
    double *_p;

    void calculateRHS()
    {
        memset(_rhs, 0, m * n * sizeof(double));

        for (int iy = 0; iy < n; iy++)
        {
            for (int ix = 0; ix < m; ix++)
            {
                int idx = iy * m + ix;
                _rhs[idx] += (_u->at(ix + 1, iy) - _u->at(ix, iy)) / dX;
                _rhs[idx] += (_v->at(ix, iy + 1) - _v->at(ix, iy)) / dX;
            }
        }
    }

    // Circle Velocity Field?
    void test_initialize_CircleVelocityField()
    {
        for (int iy = 0; iy < n; iy++)
        {
            for (int ix = 0; ix < m; ix++)
            {
                double r_square = (ix - m / 2.0) * (ix - m / 2.0) + (iy - n / 2.0) * (iy - n / 2.0);

                if (r_square < 0.01)
                {
                    (_u->at(ix, iy) = 0);
                    (_v->at(ix, iy) = 0);
                }
                else
                {
                    double r = sqrt(r_square);
                    _u->at(ix, iy) = -(iy - n / 2.0) / r;
                    _v->at(ix, iy) = (ix - m / 2.0) / r;
                }
            }
        }
    }

    void test_initialize_circleDensityField()
    {
        for (int iy = 0; iy < n; iy++)
        {
            for (int ix = 0; ix < m; ix++)
            {
                double r_square = (ix - m / 2.0) * (ix - m / 2.0) + (iy - n / 2.0) * (iy - n / 2.0);

                if (r_square <= 900)
                {
                    _d->at(ix, iy) = 1;
                }

                if (r_square == 900)
                {
                    double r = sqrt(r_square);
                    _u->at(ix, iy) = -(iy - n / 2.0) / r;
                    _v->at(ix, iy) = (ix - m / 2.0) / r;
                }
            }
        }
    }

    void project(int limit, double timestep)
    {
        double scale = timestep / rho / dX / dX;

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
                printf("Converge with %d iteration, with Norm1 = %f\n", iteration, maxDelta);
                return;
            }
        }
        printf("Exceed Limit of %d, with Norm1 = %f\n", limit, maxDelta);
    }

    void setBoundaryCondition()
    {
        for (int x = 0; x < m; x++)
        {
            _v->at(x, 0) = 0.0;
            _v->at(x, n) = 0.0;
        }

        for (int y = 0; y < m; y++)
        {
            _u->at(m, y) = 0.0;
            _u->at(0, y) = 0.0;
        }
    }

    void applyPressure(double timestep)
    {
        double scale = timestep / rho / dX;
        for (int iy = 0; iy < n; iy++)
        {
            for (int ix = 0; ix < m; ix++)
            {
                int idx = iy * m + ix;

                _u->at(ix, iy) += _p[idx] * scale;
                _u->at(ix + 1, iy) -= _p[idx] * scale;
                _v->at(ix, iy) += _p[idx] * scale;
                _v->at(ix, iy + 1) -= _p[idx] * scale;
            }
        }

        setBoundaryCondition();
    }

  public:
    FluidSolver(int m, int n, double rho)
        : m(m), n(n), rho(rho)
    {
        _d = new FluidQuantity(m, n, -1);
        _u = new FluidQuantity(m + 1, n, 0);
        _v = new FluidQuantity(m, n + 1, 1);

        _rhs = new double[m * n];
        _p = new double[m * n];

        dX = 1.0 / min(m, n);

        memset(_p, 0, m * n * sizeof(double));
    }

    ~FluidSolver()
    {
        delete _d;
        delete _u;
        delete _v;

        delete[] _rhs;
        delete[] _p;
    }

    void test_projection_tent(double timestep)
    {
        // Projection
        calculateRHS();
        project(1000, timestep);
        applyPressure(timestep);
        
        
        for(int i = 0; i < m * n; i++)
        {
            printf("%.2f, ", _p[i]);
            
            if (i % m ==0) {
                printf("\n");
            }     
        }
    }

    void update(double timestep)
    {

        if (TEST_CIRCLE)
        {
            test_initialize_CircleVelocityField();
        }

        if (TEST_CIRCLE_DENSITY)
        {
            test_initialize_circleDensityField();
        }

        if (!TEST_CIRCLE)
        {
            // Projection
            calculateRHS();
            project(1000, timestep);
            applyPressure(timestep);
        }

        _d->advect(timestep, *_u, *_v);
        _u->advect(timestep, *_u, *_v);
        _v->advect(timestep, *_u, *_v);

        _d->flip();
        _u->flip();
        _v->flip();
    }

    void addInflow(double x, double y, double x1, double y1, double d, double u, double v)
    {
        _d->addInflow(x, y, x1, y1, d);
        _u->addInflow(x, y, x1, y1, u);
        _v->addInflow(x, y, x1, y1, v);
    }

    /* Convert fluid density to RGB color scaled in 0 - 1 */
    double toRGB(int x, int y)
    {
        int idx = y * m + x;
        return max(min(1.0 - _d->Phi()[idx], 1.0), 0.0);
    }
};
