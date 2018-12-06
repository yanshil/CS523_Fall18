#include <string.h>
#include <math.h>
#include <iostream>
#include <stdio.h>
#include "FluidSimulator_Grid.h"

namespace Nova
{
template <class T, int d>
class FluidQuantity
{
    using TV = Vector<T, d>;
    using T_INDEX = Vector<int, d>;

    T *_Phi;
    T *_Phi_new;

    FluidSimulator_Grid<T, d> *grid;
    T_INDEX simulation_counts;
    TV faceOffset;
    int size;
    // Face indicator
    int axis;
    // Grid Cell size
    T hx;

    // size grid
    int m;
    int n;

    T ox, oy;

    // 1D Linear Interpolate between a and b in (0,1)
    T linp(T a, T b, T delta) const
    {
        return a * (1 - delta) + b * delta;
    }

  public:
    FluidQuantity(FluidSimulator_Grid<T, d> &grid, int axis)
        : grid(&grid), axis(axis)
    {
        m = grid.counts[0];
        n = grid.counts[1];

        if (axis == 0)
            m++;
        if (axis == 1)
            n++;

        this->simulation_counts = T_INDEX(grid.counts);
        this->hx = grid.hx;
        this->faceOffset = TV(0.5);

        if (axis != -1)
        {
            this->faceOffset(axis) = 0;
            this->simulation_counts(axis) += 1;
        }

        this->size = simulation_counts.Product();

        ox = (axis == 0) ? 0 : 0.5;
        oy = (axis == 1) ? 0 : 0.5;

        _Phi = new T[size];
        _Phi_new = new T[size];

        memset(_Phi, 0, size * sizeof(T));
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

    const T *Phi() const
    {
        return _Phi;
    }

    T at(int x, int y) const
    {
        return _Phi[x + y * m];
    }

    T at(T_INDEX &index) const
    {
        return _Phi[index2offset(index)];
    }

    T &at(int x, int y)
    {
        return _Phi[x + y * m];
    }

    T &at(T_INDEX &index)
    {
        return _Phi[index2offset(index)];
    }

    // Linear Interpolate on grid at index (x, y) (can be 0.5 if on face)
    T linp(T x, T y) const
    {
        // Clamp Coordinates
        // Extra 0.0001 for not Clamp to m+1 when x is exactly a int.
        x = std::max(x - ox, 0.0);
        x = std::min(x, m - 1.0001);
        y = std::max(y - oy, 0.0);
        y = std::min(y, n - 1.0001);

        // Get Index
        int ix = (int)x;
        int iy = (int)y;
        // Get the offset between (0,1)
        x -= ix;
        y -= iy;

        T x0 = linp(at(ix, iy), at(ix + 1, iy), x);
        T x1 = linp(at(ix, iy + 1), at(ix + 1, iy + 1), x);
        T re = linp(x0, x1, y);

        // if (oy < 0.0001) {
        //     printf("x0 = %0.3f, x1 = %0.3f\n", x0, x1);
        // printf("at(%d, %d) = %f, Next_0: at(%d, %d) = %f\n", ix, iy, at(ix, iy),ix+1, iy, at(ix+1, iy));
        // }

        return re;
    }

    void advect(T timestep, FluidQuantity *_v[d])
    {
        // Reduce Spatial locality

        // for (int idx = 0; idx < size; idx++)
        // {
        //     T_INDEX index = offset2index(idx);
        //     // TV location = (axis == -1) ? grid->Center(index) : grid->Face(axis, index);
        //     // TV location = (TV)index + faceOffset - TV(1);
        //     TV location = TV{index(0)-1, index(1)-1} + faceOffset;

        //     for (int axis = 0; axis < d; axis++)
        //     {
        //         T vtmp = _v[axis]->linp(location(0), location(1)) / hx;
        //         // Traceback
        //         location[axis] -= vtmp * timestep;
        //     }

        //     _Phi_new[idx] = linp(location(0), location(1));
        // }

        for (int iy = 0; iy < n; iy++)
        {
            for (int ix = 0; ix < m; ix++)
            {
                int idx = iy * m + ix;
                T_INDEX index = offset2index(idx);
                // printf("(%d, %d) = idx = %d, 2index = (%d, %d)\n", ix, iy, idx, index(0), index(1));

                // TV location = TV{ix, iy} + faceOffset;
                TV location = (TV)index - TV(1) + faceOffset;

                // Compute velocity

                for (int axis = 0; axis < d; axis++)
                {
                    T vtmp = _v[axis]->linp(location(0), location(1)) / hx;
                    // Traceback
                    location[axis] -= vtmp * timestep;
                }

                _Phi_new[idx] = linp(location(0), location(1));
            }
        }
    }

    void addInflow(int ix, int iy, T value)
    {
        // if(ix >= 0 & ix < m & iy >=0 & iy < n)
        if (fabs(_Phi[iy * m + ix]) < fabs(value))
        {
            _Phi[iy * m + ix] = value;
            //printf("Add Inflow at (%d, %d)\n", ix, iy);
        }
    }

    int index2offset(const T_INDEX &index) const
    {
        return grid->index2offset(index, simulation_counts);
    }

    T_INDEX offset2index(const int os) const
    {
        return grid->offset2index(os, simulation_counts);
    }
};
} // namespace Nova