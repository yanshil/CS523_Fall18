//!#####################################################################
//! \file FluidQuantity.cpp
//!#####################################################################
#include "FluidQuantity.h"
#include <nova/Tools/Utilities/Range_Iterator.h>

using namespace Nova;

/**
 * Linear INterPolate for current fluid quantity in a specific location
 */
template <typename T, int d>
T FluidQuantity<T, d>::linp(TV location) const
{
    // Fix with face offset
    TV faceOffset = TV(0.5);
    if (axis != -1)
        faceOffset(axis) = 0;
    faceOffset *= (grid->dX);

    // Fix location with Face/Scalar Indicator. True faceOffset = 0.5 * dX for each adjustification.
    location -= faceOffset;

    // Clamp location -> Inside Simulation Domain
    // TODO: What if clamped to a location that is the edge of simulation_domain and clamp_to_cell is out of bound?
    TV clamped_location = simulation_domain.Clamp(location);

    T_INDEX index = grid->Clamp_To_Cell(clamped_location);

    // TV c000 = (axis == -1) ? grid->Center(index) : grid->Face(axis, index);
    TV c000 = grid->Node(index);

    TV offsets = (clamped_location - c000) * grid->one_over_dX;

    T x0 = linp(at(index), at(grid->Next_Cell(0, index)), offsets[0]);
    T x1 = linp(at(grid->Next_Cell(1, index)), at(grid->Next_Cell(1, grid->Next_Cell(0, index))), offsets[0]);
    T re = linp(x0, x1, offsets[1]);

    return re;
}
/**
 * Advection 
 */
template <typename T, int d>
void FluidQuantity<T, d>::advect(T timestep, FluidQuantity *_v[d])
{
    for (int idx = 0; idx < size; idx++)
    {
        T_INDEX index = offset2index(idx);
        TV velocity;
        TV location = (axis == -1) ? grid->Center(index) : grid->Face(axis, index);

        for (int axis = 0; axis < d; axis++)
        {
            velocity(axis) = _v[axis]->linp(location);
        }

        // Traceback
        location -= velocity * timestep;

        _Phi_new[idx] = linp(location);
    }
}
/**
 * Add inflow 
 * 
 */
template <typename T, int d>
void FluidQuantity<T, d>::addInflow(const T_INDEX &min_corner, const T_INDEX &max_corner, const T input)
{
    Range<int, d> range(min_corner, max_corner);

    T_INDEX currIndex;
    for (Range_Iterator<d> iterator(range); iterator.Valid(); iterator.Next())
    {
        currIndex = T_INDEX() + iterator.Index();

        if (std::fabs(at(currIndex)) < std::fabs(input))
        {
            // std::cout << "addinflow(" << currIndex << ") = " << input << std::endl;
            at(currIndex) = input;
        }
    }
}
//######################################################################
template class Nova::FluidQuantity<float, 2>;
template class Nova::FluidQuantity<float, 3>;
#ifdef COMPILE_WITH_DOUBLE_SUPPORT
template class Nova::FluidQuantity<double, 2>;
template class Nova::FluidQuantity<double, 3>;
#endif