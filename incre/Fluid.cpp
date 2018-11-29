
#include <math.h>
#include <stdio.h>
#include "../lodepng/lodepng.h"

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
    FluidQuantity(int m, int n, int axis, double dX)
        : m(n), n(n), axis(axis), dX(dX)
    {
        ox = (axis == 0) ? 0 : 0.5;
        oy = (axis == 1) ? 0 : 0.5;
        dX = 1.0 / min(m, n);

        _Phi = new double[m * n];
        _Phi_new = new double[m * n];
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
        printf("at: x = %d, y = %d \n", x, y);
        return _Phi[x + y * m];
    }

    double &at(int x, int y)
    {
        printf("&at: x = %d, y = %d \n", x, y);
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
                printf("Advect: x = %f, y = %f, ox = %f \n", x, y, ox);
                double velocity_u = u.linp(x, y) / dX;
                double velocity_v = v.linp(x, y) / dX;

                // Traceback
                x -= velocity_u * timestep;
                y -= velocity_v * timestep;
                printf("TraceBack: x = %f, y = %f, timestep = %f \n", x, y, timestep);

                _Phi_new[idx] = linp(x, y);
            }
        }
    }

    // TODO
    void addInflow(double x0, double y0, double x1, double y1, double v)
    {
        int ix0 = (int)(x0 / dX - ox);
        int iy0 = (int)(y0 / dX - oy);
        int ix1 = (int)(x1 / dX - ox);
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
    double density;

    double *_rhs;
    double *_p;

    void calculateRHS()
    {
        for (int iy = 0; iy < n; iy++)
        {
            for (int ix = 0; ix < m; ix++)
            {
                int idx = iy * m + ix;
                _rhs[idx] = 0;
                _rhs[idx] -= (_u->at(ix + 1, iy) - _u->at(ix, iy)) / dX;
                _rhs[idx] -= (_v->at(ix, iy + 1) - _v->at(ix, iy)) / dX;
            }
        }
    }

    void project(int limit, double timestep)
    {
        double scale = 1.0 / dX / dX * timestep;

        double maxDelta;

        for (int iteration = 0; iteration < limit; iteration++)
        {
            maxDelta = __DBL_MIN__;
            for (int iy = 0; iy < n; iy++)
            {
                for (int ix = 0; ix < m; ix++)
                {
                    // printf("x = %d, y = %d \n", ix, iy);

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
                        sum -= scale * _p[idx + 1]; // Previous u
                    }
                    if (iy < n - 1)
                    {
                        Aii += scale;
                        sum -= scale * _p[idx + m]; // Previous u
                    }

                    double newP = (_rhs[idx] - sum) / Aii;

                    maxDelta = max(maxDelta, fabs(_p[idx] - newP));

                    _p[idx] = newP;
                }
            }
            if (maxDelta < 1e-5)
            {
                printf("Exiting solver after %d iterations, maximum change is %f\n", iteration, maxDelta);
                return;
            }
        }
        printf("Exceeded budget of %d iterations, maximum change was %f\n", limit, maxDelta);
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
        double scale = 1.0 / dX * timestep;
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
    FluidSolver(int m, int n, double density)
        : m(m), n(n), density(density)
    {
        _d = new FluidQuantity(m, n, -1, dX);
        _u = new FluidQuantity(m + 1, n, 0, dX);
        _v = new FluidQuantity(m, n + 1, 1, dX);

        _rhs = new double[m * n];
        _p = new double[m * n];
    }
    ~FluidSolver()
    {
        delete _d;
        delete _u;
        delete _v;

        delete[] _rhs;
        delete[] _p;
    }

    void update(double timestep)
    {
        calculateRHS();
        printf("Start Projection");
        project(600, timestep);
        printf("Start Apply Pressure");
        applyPressure(timestep);

        printf("Start Advection");
        _d->advect(timestep, *_u, *_v);
        _u->advect(timestep, *_u, *_v);
        _v->advect(timestep, *_u, *_v);

        _d->flip();
        _u->flip();
        _v->flip();
    }

    // TODO
    void addInflow(double x, double y, double w, double h, double d, double u, double v)
    {
        _d->addInflow(x, y, x + w, y + h, d);
        _u->addInflow(x, y, x + w, y + h, u);
        _v->addInflow(x, y, x + w, y + h, v);
    }

    // TODO
    /* Convert fluid density to RGBA image */
    void toImage(unsigned char *rgba)
    {
        for (int i = 0; i < m * n; i++)
        {
            int shade = (int)((1.0 - _d->Phi()[i]) * 255.0);
            shade = max(min(shade, 255), 0);

            rgba[i * 4 + 0] = shade;
            rgba[i * 4 + 1] = shade;
            rgba[i * 4 + 2] = shade;
            rgba[i * 4 + 3] = 0xFF;
        }
    }
};

main(int argc, char const *argv[])
{
    /* Play with these constants, if you want */
    const int sizeX = 128;
    const int sizeY = 128;

    const double density = 0.1;
    const double timestep = 0.005;

    unsigned char *image = new unsigned char[sizeX * sizeY * 4];

    FluidSolver *solver = new FluidSolver(sizeX, sizeY, density);

    double time = 0.0;
    int iterations = 0;

    while (time < 8.0)
    {
        /* Use four substeps per iteration */
        for (int i = 0; i < 4; i++)
        {
            solver->addInflow(0.45, 0.2, 0.1, 0.01, 1.0, 0.0, 3.0);
            solver->update(timestep);
            time += timestep;
            fflush(stdout);
        }

        solver->toImage(image);

        char path[256];
        sprintf(path, "Frame%05d.png", iterations++);
        lodepng_encode32_file(path, image, sizeX, sizeY);
    }
    return 0;
}
