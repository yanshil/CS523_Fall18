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
    Range<T, d> simulation_range;

    TV faceOffset;
    int size;
    // Face indicator
    int axis;
    // Grid Cell size
    T hx;

    // BC Flag
    int *_BCFlag;

    // 1D Linear Interpolate between a and b in (0,1)
    T linp(T a, T b, T delta) const
    {
        return a * (1 - delta) + b * delta;
    }

  public:
    FluidQuantity(FluidSimulator_Grid<T, d> &grid, int axis)
        : grid(&grid), axis(axis)
    {
        this->simulation_counts = T_INDEX(grid.counts);
        this->hx = grid.hx;
        this->faceOffset = TV(0.5);
        this->simulation_range = grid.domain;

        if (axis != -1)
        {
            this->faceOffset(axis) = 0;
            this->simulation_counts(axis) += 1;
            this->simulation_range.max_corner(axis) += grid.dX(axis);
        }

        this->size = simulation_counts.Product();

        _Phi = new T[size];
        _Phi_new = new T[size];
        _BCFlag = new int[size];

        memset(_Phi, 0, size * sizeof(T));
    }

    ~FluidQuantity()
    {
        delete[] _Phi;
        delete[] _Phi_new;
        delete[] _BCFlag;
    }

    void flip()
    {
        std::swap(_Phi, _Phi_new);
    }

    const T *Phi() const
    {
        return _Phi;
    }

    T at(const T_INDEX &index) const
    {
        return _Phi[index2offset(index)];
    }

    T &at(const T_INDEX &index)
    {
        return _Phi[index2offset(index)];
    }

    void calculateBCFlag()
    {
        // 0 indicate Inertia, 1 is exteria
        memset(_BCFlag, 0, size * sizeof(int));

        if (axis != -1)
        {
            for (int idx = 0; idx < size; idx++)
            {
                T_INDEX index = offset2index(idx);

                if ((index(axis) == 1) | (index(axis) == simulation_counts(axis)))
                {
                    _BCFlag[idx] = 1;
                }
            }
        }
    }

    void setBoundaryValue()
    {

        for (int idx = 0; idx < size; idx++)
        {
            T_INDEX index = offset2index(idx);

            if (_BCFlag[idx] == 1)
            {
                at(index) = 0.0;
            }
        }
    }
    //-------------------Interpolate by Index ----------------------------
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

        T_INDEX c000, c100, c010, c110;
        c000 = ixyz + T_INDEX(1);
        c100 = grid->Next_Cell(0, c000);
        c010 = grid->Next_Cell(1, c000);
        c110 = grid->Next_Cell(0, c010);

        T x00 = linp(at(c000), at(c100), offset[0]);
        T x10 = linp(at(c010), at(c110), offset[0]);
        T y0 = linp(x00, x10, offset[1]);

        if (d == 2)
        {
            return y0;
        }

        // ----------------- If d = 3----------------------

        T_INDEX c001, c101, c011, c111;

        c001 = grid->Next_Cell(2, c000);
        c101 = grid->Next_Cell(2, c100);
        c011 = grid->Next_Cell(2, c010);
        c111 = grid->Next_Cell(2, c110);

        T x01 = linp(at(c001), at(c101), offset[0]);
        T x11 = linp(at(c011), at(c111), offset[0]);
        T y1 = linp(x01, x11, offset[1]);

        T z = linp(y0, y1, offset[2]);

        return z;
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

    //------------------------Interpolate by Location ------------------------
    // Linear Interpolate on grid at index (x, y) (can be 0.5 if on face)
    // T linp(TV location)
    // {
    //     // Clamp Coordinates
    //     // Extra 0.0001 for not Clamp to m+1 when x is exactly a int.
    //     location -= faceOffset * grid->dX;

    //     location = grid->domain.Clamp(location);

    //     T_INDEX c000, c100, c010, c110;

    //     c000 = grid->Clamp_To_Cell(location);

    //     // TV location2 = (axis == -1)? grid->Center(c000):grid->Face(axis, c000);
    //     TV location2 = grid->Node(c000);

    //     TV offset = (location - location2) * grid->one_over_dX;

    //     c100 = grid->Next_Cell(0, c000);
    //     c010 = grid->Next_Cell(1, c000);
    //     c110 = grid->Next_Cell(0, c010);

    //     T x00 = linp(at(c000), at(c100), offset[0]);
    //     T x10 = linp(at(c010), at(c110), offset[0]);
    //     T y0 = linp(x00, x10, offset[1]);

    //     if (d == 2)
    //     {
    //         return y0;
    //     }

    //     // ----------------- If d = 3----------------------

    //     T_INDEX c001, c101, c011, c111;

    //     c001 = grid->Next_Cell(2, c000);
    //     c101 = grid->Next_Cell(2, c100);
    //     c011 = grid->Next_Cell(2, c010);
    //     c111 = grid->Next_Cell(2, c110);

    //     T x01 = linp(at(c001), at(c101), offset[0]);
    //     T x11 = linp(at(c011), at(c111), offset[0]);
    //     T y1 = linp(x01, x11, offset[1]);

    //     T z = linp(y0, y1, offset[2]);

    //     return z;
    // }

    // void advect(T timestep, FluidQuantity *_v[d])
    // {
    //     // Reduce Spatial locality

    //     for (int idx = 0; idx < size; idx++)
    //     {
    //         T_INDEX index = offset2index(idx);
    //         TV location = (axis == -1) ? grid->Center(index) : grid->Face(axis, index);

    //         for (int axis = 0; axis < d; axis++)
    //         {
    //             T vtmp = _v[axis]->linp(location);
    //             // Traceback
    //             location[axis] -= vtmp * timestep;
    //         }

    //         _Phi_new[idx] = linp(location);
    //     }
    // }
    //----------------------------------------------------------
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