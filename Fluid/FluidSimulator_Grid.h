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
    using TV = Vector<T, d>;
    using Base = Grid<T, d>;
    using T_INDEX = Vector<int, d>;
    T_INDEX whole_domain;
    int number_of_ghost_cells;

  public:
    FluidSimulator_Grid()
        : Base()
    {
        this->whole_domain = this->counts + T_INDEX(number_of_ghost_cells * 2);
    }
    FluidSimulator_Grid(const T_INDEX &counts_input, const Range<T, d> &domain_input, const int number_of_ghost_cells = 0)
        : Base(counts_input, domain_input), number_of_ghost_cells(number_of_ghost_cells)
    {
        this->whole_domain = this->counts + T_INDEX(number_of_ghost_cells * 2);
        std::cout << "whole_domain = " << whole_domain << std::endl;
    }
    ~FluidSimulator_Grid()
    {
    }

    //######################################################################
    /////////////////////////////////////////////////
    /// Get Next / Previous Cell for given axis
    /////////////////////////////////////////////////
    T_INDEX Next_Cell(const int axis, const T_INDEX &index);
    T_INDEX Previous_Cell(const int axis, const T_INDEX &index);

    T_INDEX offset2index(const int os);
    int index2offset(const T_INDEX &index);
};

template <typename T, int d>
Vector<int, d> FluidSimulator_Grid<T, d>::Next_Cell(const int axis, const T_INDEX &index)
{
    T_INDEX shifted_index(index);
    shifted_index(axis) += 1;

    return shifted_index;
}
// TODO
template <typename T, int d>
Vector<int, d> FluidSimulator_Grid<T, d>::Previous_Cell(const int axis, const T_INDEX &index)
{
    T_INDEX shifted_index(index);
    shifted_index(axis) -= 1;

    return shifted_index;
}

/**
 * index2offset
 * Get exact offset for Column-based 1D array
 * Index start from (1,1) to _counts_
 * return (z * xSize * ySize) + (y * xSize) + x;
 */
template <typename T, int d>
int FluidSimulator_Grid<T, d>::index2offset(const T_INDEX &index)
{
    // Becuase index in the grid start from (1,1)...
    T_INDEX tmp_index = index;
    tmp_index += T_INDEX(number_of_ghost_cells) - T_INDEX(1);

    int os = tmp_index[1] * whole_domain[0] + tmp_index[0];
    if (d == 3)
        os += tmp_index[2] * whole_domain[0] * whole_domain[1];
    return os;
}
/**
 * offset2index
 * return cell index for given offset
 * Index start from (1,1) to _counts_
 */
template <typename T, int d>
Vector<int, d> FluidSimulator_Grid<T, d>::offset2index(const int os)
{
    // 3D: os = z * m * n + y * m + x
    // 2D: os = y * m + x
    T_INDEX tmp_index = T_INDEX();

    // x <- os mod m
    tmp_index[0] = os % whole_domain[0];

    if (d == 2)
        // y <- (os - x) / m
        tmp_index[1] = (os - tmp_index[0]) / whole_domain[1];
    else
    {
        // y <- (os - x) mod n
        tmp_index[1] = (os - tmp_index[0]) % whole_domain[1];

        // z <- (os - x - y * m) / (m*n)
        tmp_index[2] = (os - tmp_index[0] - tmp_index[1] * whole_domain[0]) / whole_domain[0] / whole_domain[1];
    }

    // Becuase index in the grid start from (1,1)...
    tmp_index += T_INDEX(1) - T_INDEX(number_of_ghost_cells);
    return tmp_index;
}

} // namespace Nova
#endif