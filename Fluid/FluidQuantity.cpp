//!#####################################################################
//! \file FluidQuantity.cpp
//!#####################################################################
#include "FluidQuantity.h"
using namespace Nova;
////////////////////////////////////////////////////////////////////////
/// Fundamental Read & Write
////////////////////////////////////////////////////////////////////////
/// Constructor
/// (The one should not be called because grid is not initialized)
////////////////////////////////////////////////////////////////////////
template <typename T, int d>
FluidQuantity<T, d>::FluidQuantity()
{
    grid = nullptr;
    std::cout << "Wrong Constructor!" << std::endl;
}
////////////////////////////////////////////////////////////////////////
/// Constructor
////////////////////////////////////////////////////////////////////////
template <typename T, int d>
FluidQuantity<T, d>::FluidQuantity(Grid<T, d> &grid, int axis, int number_of_ghost_cells)
    : grid(&grid), axis(axis), number_of_ghost_cells(number_of_ghost_cells)
{
    std::cout << "Correct Constructor!" << std::endl;

    storing_counts = T_INDEX(grid.counts);
    if (axis != -1)
        storing_counts(axis) += 1;

    Phi = new T[storing_counts.Product()];
    Phi_new = new T[storing_counts.Product()];
}
////////////////////////////////////////////////////////////////////////
/// Destructor
////////////////////////////////////////////////////////////////////////
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
        Phi[i] = content;
}
//######################################################################
// at
//######################################################################
template <typename T, int d>
T FluidQuantity<T, d>::at(const T_INDEX &index)
{
    // e.g. For n = 4
    // Density: 4 * 4 * 4;      velocityU: 5 * 4 * 4;   velocity V: 4 * 5 * 4;  velocity W: 4 * 4 * 5
    // this->Inside_Domain are checking if is in (1,1,1) to (5, 4, 4) [For U]
    // grid->Inside_Domain(index) is checking if is in (1,1,1) to (4,4,4)
    // grid->Inside_Domain(index, num_of_ghost_cell = 1) is checking if is in (0,0,0) to (5,5,5)

    //----------------- Simulation Field Domain (d/u/v domain) ------------------
    if (!Inside_Domain(index))
    {
        TV location = grid->Center(index);
        T_INDEX clamped_index = grid->Clamp_To_Cell(location, number_of_ghost_cells);
        return Phi[index2offset(clamped_index)];
    }
    return Phi[index2offset(index)];
}
/**
 * modify_at
 * Read and Write access in the quantity field for Phi
 */
template <typename T, int d>
T &FluidQuantity<T, d>::modify_at(const T_INDEX &index)
{
    if (!Inside_Domain(index))
    {
        std::cout << "index = " << index << std::endl;
        // Raise exception
        throw std::runtime_error("Try to write at an Out_of_domain area");
    }
    // if (!grid->Inside_Domain(index))
    // {
    //     std::cout << "index = " << index << std::endl;
    //     // Raise exception
    //     throw std::runtime_error("Try to write at an Out_of_domain area");
    // }
    return Phi[index2offset(index)];
}
//######################################################################
// new_at: Read & Write access in the quantity field for Phi_new
//######################################################################
template <typename T, int d>
T &FluidQuantity<T, d>::new_at(const T_INDEX &index)
{
    if (!Inside_Domain(index))
    {
        // Raise exception
        throw std::runtime_error("Try to write at an Out_of_domain area");
    }
    // if (!grid->Inside_Domain(index))
    // {
    //     // Raise exception
    //     throw std::runtime_error("Try to write at an Out_of_domain area");
    // }
    return Phi_new[index2offset(index)];
}
//######################################################################
// rgb_at: Get RGB color scaled to (0,1)
//######################################################################
template <typename T, int d>
T FluidQuantity<T, d>::rgb_at(const T_INDEX &index)
{
    return std::max(std::min(1.0 - at(index), 1.0), 0.0);
}
/**
 * linter
 * 1D Linear Interpolate etween a and b for x in (0, 1)
 */
template <typename T, int d>
T FluidQuantity<T, d>::linter(T a, T b, T x)
{
    return (1.0 - x) * a + x * b;
}
/**
 * Inside_Simulation_Domain
 * Check if quantity inside its storing (simulation) domain
 */
template <typename T, int d>
bool FluidQuantity<T, d>::Inside_Domain(const T_INDEX &index)
{
    T_INDEX diff = storing_counts - index;

    for (int i = 0; i < d; i++)
    {
        if (diff(i) < 0)
        {
            return false;
        }
    }
    return true;
}
/**
 * index2offset
 * Get exact offset for Column-based 1D array
 * Index start from (1,1) to _counts_
 * return (z * xSize * ySize) + (y * xSize) + x;
 */
template <typename T, int d>
int FluidQuantity<T, d>::index2offset(const T_INDEX &index)
{
    // Becuase index in the grid start from (1,1)...
    T_INDEX tmp_index = index - T_INDEX(1);

    int os = tmp_index[1] * storing_counts[0] + tmp_index[0];
    if (d == 3)
        os += tmp_index[2] * storing_counts[0] * storing_counts[1];
    return os;
}
/**
 * offset2index
 * return cell index for given offset
 * Index start from (1,1) to _counts_
 */
template <typename T, int d>
Vector<int, d> FluidQuantity<T, d>::offset2index(const int os)
{
    // 3D: os = z * m * n + y * m + x
    // 2D: os = y * m + x
    T_INDEX tmp_index = T_INDEX();

    // x <- os mod m
    tmp_index[0] = os % storing_counts[0];

    if (d == 2)
        // y <- (os - x) / m
        tmp_index[1] = (os - tmp_index[0]) / storing_counts[1];
    else
    {
        // y <- (os - x) mod n
        tmp_index[1] = (os - tmp_index[0]) % storing_counts[1];

        // z <- (os - x - y * m) / (m*n)
        tmp_index[2] = (os - tmp_index[0] - tmp_index[1] * storing_counts[0]) / storing_counts[0] / storing_counts[1];
    }

    // Becuase index in the grid start from (1,1)...
    tmp_index += T_INDEX(1);
    return tmp_index;
}
/////////////////////////////////////////////////
/// Get Next / Previous Cell for given axis
/////////////////////////////////////////////////
template <typename T, int d>
Vector<int, d> FluidQuantity<T, d>::Next_Cell(const int axis, const T_INDEX &index)
{
    T_INDEX shifted_index(index);
    shifted_index(axis) += 1;

    return shifted_index;
}
template <typename T, int d>
Vector<int, d> FluidQuantity<T, d>::Previous_Cell(const int axis, const T_INDEX &index)
{
    T_INDEX shifted_index(index);
    shifted_index(axis) -= 1;

    return shifted_index;
}
//######################################################################
// computeVelocity:
//######################################################################
template <typename T, int d>
Vector<T, d> FluidQuantity<T, d>::computeVelocity(const TV &location, FluidQuantity *velocityField[d])
{
    TV velocity = TV(0.0);
    // TV location = (axis == -1) ? grid->Center(index) : grid->Face(axis, index);

    for (int i = 0; i < d; i++)
    {
        velocity[i] = velocityField[i]->linter(location);
    }

    // If Scalar (for Density)
    // if (axis == -1)
    // {
    //     for (int i = 0; i < d; i++)
    //     {
    //         velocity[i] += 0.5 * velocityField[i]->at(index);
    //         if(!Inside_Domain(Next_Cell(i, index)))
    //             velocity[i] += 0.5 * velocityField[i]->at(index);
    //         else
    //             velocity[i] += 0.5 * velocityField[i]->at(Next_Cell(i, index));
    //     }
    // }
    // else
    // {
    //     for (int i = 0; i < d; i++)
    //     {
    //         if (i == axis)
    //             velocity[i] = velocityField[i]->at(index);
    //         else
    //         {
    //             velocity[i] += 0.25 * velocityField[i]->at(index);
    //             if(!Inside_Domain(Next_Cell(i, index)))
    //                 velocity[i] += 0.25 * velocityField[i]->at(index);
    //             else
    //                 velocity[i] += 0.25 * velocityField[i]->at(Next_Cell(i, index));
    //             if(!Inside_Domain(Previous_Cell(axis, index)))
    //                 velocity[i] += 0.25 * velocityField[i]->at(index);
    //             else
    //                 velocity[i] += 0.25 * velocityField[i]->at(Previous_Cell(axis, index));
    //             if(!Inside_Domain(Next_Cell(i, Previous_Cell(axis, index))))
    //                 velocity[i] += 0.25 * velocityField[i]->at(index);
    //             else
    //                 velocity[i] += 0.25 * velocityField[i]->at(Next_Cell(i, Previous_Cell(axis, index)));
    //         }
    //     }
    // }

    // // If Scalar (for Density)
    // if (axis == -1)
    // {
    //     for (int i = 0; i < d; i++)
    //     {
    //         velocity[i] += 0.5 * velocityField[i]->at(index);
    //         velocity[i] += 0.5 * velocityField[i]->at(Next_Cell(i, index));
    //     }
    // }
    // else
    // {
    //     for (int i = 0; i < d; i++)
    //     {
    //         if (i == axis)
    //             velocity[i] = velocityField[i]->at(index);
    //         else
    //         {
    //             velocity[i] += 0.25 * velocityField[i]->at(index);
    //             velocity[i] += 0.25 * velocityField[i]->at(Next_Cell(i, index));
    //             velocity[i] += 0.25 * velocityField[i]->at(Previous_Cell(axis, index));
    //             velocity[i] += 0.25 * velocityField[i]->at(Next_Cell(i, Previous_Cell(axis, index)));
    //         }
    //     }
    // }

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

    if (!grid->Inside_Domain(tmp_index))
    {
        tmp_location = (TV)grid->domain.Clamp(tmp_location);
    }

    return tmp_location;
}

/*!
 * Linear Interpolator for TV{i, j, k} on grid
 * Coordinates will be clamped to lie in simulation domain
 */
template <typename T, int d>
T FluidQuantity<T, d>::linter(const TV &location)
{
    /*
         * (0.5, 0.5, 0.5) for center     axis = -1
         * (0, 0.5, 0.5) for face on X    axis = 0
         * (0.5, 0, 0.5) for face on Y    axis = 1
         * (0.5, 0.5, 0) for face on Z    axis = 2
         */

    TV faceOffset = TV(0.5);
    if (axis != -1)
        faceOffset(axis) = 0;
    faceOffset *= (grid->dX);

    // Fix location with Face/Scalar Indicator. True faceOffset = 0.5 * dX for each adjustification.
    TV fixed_location = location - faceOffset;

    // Find the true linear interpolate cell.
    T_INDEX c000, c100, c010, c110;
    c000 = grid->Clamp_To_Cell(fixed_location, number_of_ghost_cells);

    TV c000_fixlocation = (axis == -1) ? grid->Center(c000) : grid->Face(axis, c000);

    /* Project offset to (0, 1) and apply Linear Interpolate */
    TV offset = location - c000_fixlocation;
    offset *= grid->one_over_dX;

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
    TV location = (axis == -1) ? grid->Center(index) : grid->Face(axis, index);
    TV velocity = computeVelocity(location, velocityField);

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
template class Nova::FluidQuantity<float, 2>;
template class Nova::FluidQuantity<float, 3>;
#ifdef COMPILE_WITH_DOUBLE_SUPPORT
template class Nova::FluidQuantity<double, 2>;
template class Nova::FluidQuantity<double, 3>;
#endif