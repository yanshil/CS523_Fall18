//!#####################################################################
//! \file FluidQuantity.cpp
//!#####################################################################
#include "FluidQuantity.h"
using namespace Nova;
//######################################################################
// Constructor
//######################################################################
template <typename T, int d>
FluidQuantity<T, d>::FluidQuantity(Grid<T, d> &grid, int axis, int number_of_ghost_cells)
{
    std::cout << "Correct Constructor!" << std::endl;

    this->grid = &grid;
    this->number_of_ghost_cells = number_of_ghost_cells;
    this->axis = axis;

    storing_counts = T_INDEX(grid.counts);
    if (axis != -1)
        storing_counts(axis) += 1;

    Phi = new T[storing_counts.Product()];
    Phi_new = new T[storing_counts.Product()];
}
//######################################################################
// Destructor
//######################################################################
template <typename T, int d>
FluidQuantity<T, d>::~FluidQuantity()
{
    delete[] Phi;
    delete[] Phi_new;
}

//######################################################################
// fill
//######################################################################
template <typename T, int d>
void FluidQuantity<T, d>::fill(T content)
{
    for (int i = 0; i < storing_counts.Product(); i++)
    {
        Phi[i] = content;
    }
}

//######################################################################
// computeVelocity
//######################################################################
template <typename T, int d>
Vector<T, d> FluidQuantity<T, d>::computeVelocity(const T_INDEX &index, FluidQuantity *velocityField[d])
{
    TV velocity = TV(0.0);

    // If Scalar (for Density)
    if (axis == -1)
    {
        for (int i = 0; i < d; i++)
        {
            velocity[i] += 0.5 * velocityField[i]->at(index);
            velocity[i] += 0.5 * velocityField[i]->at(Next_Cell(i, index));
        }
    }
    else
    {
        for (int i = 0; i < d; i++)
        {
            if (i == axis)
                velocity[i] = velocityField[i]->at(index);
            else
            {
                velocity[i] += 0.25 * velocityField[i]->at(index);
                velocity[i] += 0.25 * velocityField[i]->at(Next_Cell(i, index));
                velocity[i] += 0.25 * velocityField[i]->at(Previous_Cell(axis, index));
                velocity[i] += 0.25 * velocityField[i]->at(Next_Cell(i, Previous_Cell(axis, index)));
            }
        }
    }

    return velocity;
}
//######################################################################
// TODO
//######################################################################
template <typename T, int d> 
Vector<T, d> FluidQuantity<T, d>::Clamp_To_Domain(const TV &location)
{
    TV tmp_location = TV(location);

    // Clamp to domain if not in domain
    T_INDEX tmp_index = grid->Cell(tmp_location, number_of_ghost_cells);

    if (!(*grid).Inside_Domain(tmp_index))
    {
        tmp_location = (TV)(*grid).domain.Clamp(tmp_location);
    }

    return tmp_location;
}

/*!
 * Linear Interpolator for TV{i, j, k} on grid
 * Coordinates will be clamped to lie in simulation domain
*/
//######################################################################
// TODO
//######################################################################
template <typename T, int d> 
T FluidQuantity<T, d>::linter(TV &location)
{
    /*
         * (0.5, 0.5, 0.5) for center     axis = -1
         * (0, 0.5, 0.5) for face on X    axis = 0
         * (0.5, 0, 0.5) for face on Y    axis = 1
         * (0.5, 0.5, 0) for face on Z    axis = 2
         */

    T_INDEX locationIndex = (*grid).Clamp_To_Cell(location, number_of_ghost_cells);
    TV clocation = (*grid).Node(locationIndex);

    TV faceOffset = TV(0.5);

    if (axis != -1)
        faceOffset(axis) = 0;

    faceOffset *= ((*grid).dX);

    // Fix location with Face/Scalar Indicator. True faceOffset = 0.5 * dX for each adjustification.
    TV fixed_location = location - faceOffset;

    // Find the true linear interpolate cell.
    T_INDEX c000, c100, c010, c110;
    c000 = (*grid).Clamp_To_Cell(fixed_location, number_of_ghost_cells);

    TV c000_fixlocation = (axis == -1) ? (*grid).Center(c000) : (*grid).Face(axis, c000);

    /* Project offset to (0, 1) and apply Linear Interpolate */
    TV offset = location - c000_fixlocation;
    offset *= (*grid).one_over_dX;

    c100 = Next_Cell(0, c000);
    c010 = Next_Cell(1, c000);
    c110 = Next_Cell(0, c010);
    T px00 = linter(at(c000), at(c100), offset[0]);
    T px10 = linter(at(c010), at(c110), offset[0]);
    T py0 = linter(px00, px10, offset[1]);

    if (d == 2)
    {
        return py0;
    }

    // ----------------- If d = 3----------------------

    T_INDEX c001, c101, c011, c111;

    c001 = Next_Cell(2, c000);
    c101 = Next_Cell(2, c100);
    c011 = Next_Cell(2, c010);
    c111 = Next_Cell(2, c110);

    T px01 = linter(at(c001), at(c101), offset[0]);
    T px11 = linter(at(c011), at(c111), offset[0]);
    T py1 = linter(px01, px11, offset[1]);

    T pz = linter(py0, py1, offset[2]);

    return pz;
}
//######################################################################
// TODO
//######################################################################
template <typename T, int d> 
void FluidQuantity<T, d>::advect(const T_INDEX &index, T timestep, FluidQuantity *velocityField[d])
{
    TV velocity = computeVelocity(index, velocityField);

    TV location = (axis == -1) ? (*grid).Center(index) : (*grid).Face(axis, index);

    TV location_traceback = Clamp_To_Domain(location - timestep * velocity);

    new_at(index) = linter(location_traceback);
}
//######################################################################
// TODO
//######################################################################
template <typename T, int d> 
void FluidQuantity<T, d>::advect(T timestep, FluidQuantity *velocityField[d])
{
    T_INDEX currIndex;
    for (Range_Iterator<d> iterator(Range<int, d>(T_INDEX(1), storing_counts)); iterator.Valid(); iterator.Next())
    {
        currIndex = T_INDEX() + iterator.Index();
        advect(currIndex, timestep, velocityField);
    }
}
//######################################################################
template class Nova::FluidQuantity<float,2>;
template class Nova::FluidQuantity<float,3>;
#ifdef COMPILE_WITH_DOUBLE_SUPPORT
template class Nova::FluidQuantity<double,2>;
template class Nova::FluidQuantity<double,3>;
#endif