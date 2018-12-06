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

    T at(const T_INDEX &index) const
    {
        return _Phi[index2offset(index)];
    }

    T &at(int x, int y)
    {
        return _Phi[x + y * m];
    }

    T &at(const T_INDEX &index)
    {
        return _Phi[index2offset(index)];
    }

    // Linear Interpolate on grid at index (x, y) (can be 0.5 if on face)
    T linp(TV location)
    {
        // Clamp Coordinates
        // Extra 0.0001 for not Clamp to m+1 when x is exactly a int.
        location -= faceOffset;

        // TODO
        for (int axis = 0; axis < d; axis++)
        {
            location(axis) = std::max(location(axis), 0.0);
            location(axis) = std::min(location(axis), simulation_counts(axis) - 1.0001);
        }

        T_INDEX ixyz(0);

        for (int axis = 0; axis < d; axis++)
        {
            ixyz(axis) = (int)location(axis);
        }

        TV offset = location - (TV)ixyz;
        T_INDEX c000 = ixyz + T_INDEX(1);

        T x0 = linp(at(c000), at(grid->Next_Cell(0, c000)), offset[0]);
        T x1 = linp(at(grid->Next_Cell(1, c000)), at(grid->Next_Cell(1, grid->Next_Cell(0, c000))), offset[0]);
        T re = linp(x0, x1, offset[1]);

        return re;

        // TODO: 3D
    }

    void advect(T timestep, FluidQuantity *_v[d])
    {
        // Reduce Spatial locality

        for (int idx = 0; idx < size; idx++)
        {
            T_INDEX index = offset2index(idx);
            TV location = (TV)index - TV(1) + faceOffset;

            for (int axis = 0; axis < d; axis++)
            {
                T vtmp = _v[axis]->linp(location) / hx;
                // Traceback
                location[axis] -= vtmp * timestep;
            }

            _Phi_new[idx] = linp(location);
        }
    }

    void addInflow(const T_INDEX &index, T value)
    {
        // if(ix >= 0 & ix < m & iy >=0 & iy < n)
        if (fabs(at(index)) < fabs(value))
            at(index) = value;

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