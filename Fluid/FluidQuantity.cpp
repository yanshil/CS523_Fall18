//!#####################################################################
//! \file FluidQuantity.cpp
//!#####################################################################
#include "FluidQuantity.h"

using namespace Nova;

/**
 * Linear INterPolate for current fluid quantity in a specific location
 */
template <typename T, int d>
T FluidQuantity<T, d>::linp(TV location) const
{
    // // Fix with face offset
    // TV faceOffset = TV(0.5);
    // if (axis != -1)
    //     faceOffset(axis) = 0;
    // faceOffset *= (grid->dX);

    // // Fix location with Face/Scalar Indicator. True faceOffset = 0.5 * dX for each adjustification.
    // location -= faceOffset;

    // // Clamp location -> Inside Simulation Domain
    // // TODO: What if clamped to a location that is the edge of simulation_domain and clamp_to_cell is out of bound?
    // // TODO: Fake Shrink for the domain in .h
    // TV clamped_location = simulation_domain.Clamp(location);

    // T_INDEX index = grid->Clamp_To_Cell(clamped_location);

    // // TV c000 = (axis == -1) ? grid->Center(index) : grid->Face(axis, index);
    // TV c000 = grid->Node(index);

    // TV offsets = (clamped_location - c000) * grid->one_over_dX;

    // T x0 = linp(at(index), at(grid->Next_Cell(0, index)), offsets[0]);
    // T x1 = linp(at(grid->Next_Cell(1, index)), at(grid->Next_Cell(1, grid->Next_Cell(0, index))), offsets[0]);
    // T re = linp(x0, x1, offsets[1]);

    // return re;

    // TODO: 3D

    double x = location(0);
    double y = location(1);

    TV faceOffset = TV(0.5);
    if (axis != -1)
        faceOffset(axis) = 0;

    T ox = faceOffset(0);
    T oy = faceOffset(1);

    int m = simulation_counts(0);
    int n = simulation_counts(1);

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

    T_INDEX index = T_INDEX{ix, iy} + T_INDEX(1);

    T x0 = linp(at(index), at(grid->Next_Cell(0, index)), x);
    T x1 = linp(at(grid->Next_Cell(1, index)), at(grid->Next_Cell(1, grid->Next_Cell(0, index))), x);
    T re = linp(x0, x1, y);

    // double x0 = linp(at(ix, iy), at(ix + 1, iy), x);
    // double x1 = linp(at(ix, iy + 1), at(ix + 1, iy + 1), x);
    // double re = linp(x0, x1, y);
    
    // if (axis == 1) {
    //     printf("x0 = %0.3f, x1 = %0.3f\n", x0, x1);

    //     std::cout << "at("<<index <<") = " << at(index) << ", ";
    //     std::cout << "Next_0: at("<< grid->Next_Cell(0, index) <<") = " << at(grid->Next_Cell(0, index)) <<std::endl;
    // }
    

    return re;
}
/**
 * Advection 
 */
template <typename T, int d>
void FluidQuantity<T, d>::advect(T timestep, FluidQuantity *_v[d])
{
    // for (int idx = 0; idx < size; idx++)
    // {
    //     T_INDEX index = offset2index(idx);
    //     TV velocity;
    //     // TV location = (axis == -1) ? grid->Center(index) : grid->Face(axis, index);

    //     TV faceOffset = TV(0.5);
    //     if (axis != -1)
    //         faceOffset(axis) = 0;

    //     TV location = (TV)index + faceOffset - TV(1);

    //     for (int axis = 0; axis < d; axis++)
    //     {
    //         velocity(axis) = _v[axis]->linp(location) / hx;
    //     }

    //     // Traceback
    //     location -= velocity * timestep;

    //     _Phi_new[idx] = linp(location);
    // }

    // Reduce Spatial locality

    TV faceOffset = TV(0.5);
    if (axis != -1)
        faceOffset(axis) = 0;

    T ox = faceOffset(0);
    T oy = faceOffset(1);

    int m = simulation_counts(0);
    int n = simulation_counts(1);

    for (int iy = 0; iy < n; iy++)
    {
        for (int ix = 0; ix < m; ix++)
        {
            int idx = iy * m + ix;

            T x = ix + ox;
            T y = iy + oy;

            // Compute velocity
            T velocity_u = _v[0]->linp(TV{x, y}) / hx;
            T velocity_v = _v[1]->linp(TV{x, y}) / hx;

            // Traceback
            x -= velocity_u * timestep;
            y -= velocity_v * timestep;

            _Phi_new[idx] = linp(TV{x, y});
            // printf("ox=%0.2f, oy = %0.2f, v = (%0.3f, %0.3f), linp(%f, %f) = %f\n", ox, oy,velocity_u, velocity_v, x, y, _Phi_new[idx]);
        }
    }
}
/**
 * Add inflow 
 * 
 */
template <typename T, int d>
void FluidQuantity<T, d>::addInflow(const T_INDEX &index, const T input)
{
    if (std::fabs(at(index)) < std::fabs(input))
    {
        at(index) = input;
        // printf("Add Inflow at (%d, %d)\n", index(0) - 1, index(1) - 1);
    }
}
//######################################################################
template class Nova::FluidQuantity<float, 2>;
template class Nova::FluidQuantity<float, 3>;
#ifdef COMPILE_WITH_DOUBLE_SUPPORT
template class Nova::FluidQuantity<double, 2>;
template class Nova::FluidQuantity<double, 3>;
#endif