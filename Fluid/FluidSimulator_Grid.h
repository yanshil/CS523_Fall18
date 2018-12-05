
//!#####################################################################
//! \file FluidSimulator_Grid.h
//!#####################################################################
// Class FluidSimulator_Grid
//######################################################################
#ifndef __FluidSimulator_Grid__
#define __FluidSimulator_Grid__

#include <nova/Tools/Grids/Grid.h>

namespace Nova
{
template <class T, int d>
class FluidSimulator_Grid : public Grid<T, d>
{
  private:
    using Base = Grid<T, d>;
    using TV = Vector<T, d>;
    using T_INDEX = Vector<int, d>;

  public:
    // Max Grid Cell size
    T hx;

    FluidSimulator_Grid()
        : Base()
    {
    }

    FluidSimulator_Grid(const T_INDEX &counts_input, const Range<T, d> &domain_input)
        : Base(counts_input, domain_input)
    {
        this->hx = this->dX.Max();
    }
    ~FluidSimulator_Grid()
    {
    }

    T_INDEX Next_Cell(const int axis, const T_INDEX &index) const
    {
        T_INDEX n_index = T_INDEX(index);
        n_index(axis) += 1;
        return n_index;
    }

    T_INDEX Previous_Cell(const int axis, const T_INDEX &index) const
    {
        T_INDEX p_index = T_INDEX(index);
        p_index(axis) -= 1;
        return p_index;
    }

    int index2offset(const T_INDEX &index, T_INDEX dim) const
    {
        // Becuase index in the grid start from (1,1)...
        T_INDEX tmp_index = index;
        tmp_index += -T_INDEX(1);

        int os = tmp_index[1] * dim[0] + tmp_index[0];
        if (d == 3)
            os += tmp_index[2] * dim[0] * dim[1];
        return os;
    }

    T_INDEX offset2index(const int os, T_INDEX dim) const
    {
        // 3D: os = z * m * n + y * m + x
        // 2D: os = y * m + x
        T_INDEX tmp_index = T_INDEX();

        // x <- os mod m
        tmp_index[0] = os % dim[0];

        if (d == 2)
            // y <- (os - x) / m
            tmp_index[1] = (os - tmp_index[0]) / dim[1];
        else
        {
            // y <- (os - x) mod n
            tmp_index[1] = (os - tmp_index[0]) % dim[1];

            // z <- (os - x - y * m) / (m*n)
            tmp_index[2] = (os - tmp_index[0] - tmp_index[1] * dim[0]) / dim[0] / dim[1];
        }

        // Becuase index in the grid start from (1,1)...
        tmp_index += T_INDEX(1);
        return tmp_index;
    }
};

} // namespace Nova
#endif