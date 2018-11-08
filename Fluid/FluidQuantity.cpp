#include "FluidQuantity.h"

using namespace Nova;

T_INDEX Next_Cell(const int axis, const T_INDEX &index);
T_INDEX Previous_Cell(const int axis, const T_INDEX &index);

FluidQuantity::FluidQuantity(Grid<T, d> &grid, int axis)
{

    Phi = new float[grid.counts.Product()];
    Phi_new = new float[grid.counts.Product()];
    this->grid = &grid;
    this->number_of_ghost_cells = 1;
    this->axis = axis;
}

FluidQuantity::~FluidQuantity()
{
    delete[] Phi;
    delete[] Phi_new;
}

TV FluidQuantity::computeVelocity(T_INDEX &index, FluidQuantity *FluidVelocity[d])
{
    TV velocity = TV(0.0);

    // If Scalar
    if (axis == -1)
    {
        for (int i = 0; i < d; i++)
        {
            // T_INDEX next_cell_index(index);
            // T_INDEX next_cell_index2(Next_Cell(i, next_cell_index));

            velocity[i] += 0.5 * FluidVelocity[i]->at(index);
            // velocity[i] += 0.5 * FluidVelocity[i]->at(next_cell_index2);
            velocity[i] += 0.5 * FluidVelocity[i]->at(Next_Cell(i, index));
        }
    }
    else
    {
        for (int i = 0; i < d; i++)
        {
            if (i == axis)
                velocity[i] = FluidVelocity[i]->at(index);
            else
            {
                T_INDEX next_cell_index_onI(index);
                T_INDEX next_cell_index_onI2(Next_Cell(i, next_cell_index_onI));

                T_INDEX previous_cell_index_onAxis(index);
                T_INDEX previous_cell_index_onAxis2(Previous_Cell(axis, previous_cell_index_onAxis));

                T_INDEX next_pre_index(previous_cell_index_onAxis);
                T_INDEX next_pre_index2(Next_Cell(i, next_pre_index));

                velocity[i] += 0.25 * FluidVelocity[i]->at(index);
                velocity[i] += 0.25 * FluidVelocity[i]->at(next_cell_index_onI2);
                velocity[i] += 0.25 * FluidVelocity[i]->at(previous_cell_index_onAxis2);
                velocity[i] += 0.25 * FluidVelocity[i]->at(next_pre_index2);
            }
        }
    }
}

TV FluidQuantity::traceBack(TV &location, float timestep, TV &velocity)
{
    TV tmp_location = TV(location);
    tmp_location -= timestep * velocity;

    // Clamp to domain if not in domain
    T_INDEX tmp_index = (*grid).Cell(tmp_location, number_of_ghost_cells);

    if (!(*grid).Inside_Domain(tmp_index, number_of_ghost_cells))
    {
        tmp_location = (TV)(*grid).domain.Clamp(tmp_location);
        // tmp_location = (TV)(*grid).Clamp_To_Cell(location, number_of_ghost_cells);
    }

    return tmp_location;
}

/*!
 * Linear Interpolator for TV{i, j, k} on grid
 * Coordinates will be clamped to lie in simulation domain
*/
float FluidQuantity::linter(TV &location)
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

    TV fixed_location = location - faceOffset;

    T_INDEX c000, c100, c010, c110;
    c000 = (*grid).Clamp_To_Cell(fixed_location, number_of_ghost_cells);

    /* Project offset to (0, 1) and apply Linear Interpolate */
    TV offset = fixed_location - (TV)c000;

    for (size_t i = 0; i < d; i++)
    {
        offset[i] = offset[i] * (*grid).one_over_dX[i];
    }

    c100 = Next_Cell(0, c000);
    c010 = Next_Cell(1, c000);
    c110 = Next_Cell(0, c010);

    T_INDEX c001, c101, c011, c111;
    if(d == 3)
    {
        c001 = Next_Cell(2, c000);
        c101 = Next_Cell(2, c100);
        c011 = Next_Cell(2, c010);
        c111 = Next_Cell(2, c110);
    }

    float px00 = linter(at(c000), at(c100), offset[0]);
    float px10 = linter(at(c010), at(c110), offset[0]);
    float py0 = linter(px00, px10, offset[1]);

    if(d==2)
        return py0;
    else
    {
        float px01 = linter(at(c001), at(c101), offset[0]);
        float px11 = linter(at(c011), at(c111), offset[0]);
        float py1 = linter(px00, px10, offset[1]);

        float pz = linter(py0, py1, offset[2]);

        return pz;
    }

}


void FluidQuantity::advection(float timestep, FluidQuantity * velocityField[d])
{
    // for all cell in domain with &index
    //
    T_INDEX index = T_INDEX(1);

    TV velocity = computeVelocity(index, velocityField);
    TV location;
    TV location_traceback = traceBack(location, timestep, velocity);
    Phi_new[index2offset(index)] = linter(location_traceback);
}

/*!
 * Return Index of Next Cell for the given cell index on axis = {0, 1, (,2)}
 */
T_INDEX Next_Cell(const int axis, const T_INDEX &index)
{
    T_INDEX shifted_index(index);
    shifted_index(axis) += 1;
    return shifted_index;
}

/*!
 * Return Index of Previous Cell for the given cell index on axis = {0, 1, (,2)}
 */
T_INDEX Previous_Cell(const int axis, const T_INDEX &index)
{

    T_INDEX shifted_index(index);
    shifted_index(axis) -= 1;
    return shifted_index;
}

// /*!
//  * Auxiliary Function: Get exact offset for Column-based 1D array
//  * return (z * xSize * ySize) + (y * xSize) + x;
//  */
// int index2offset(const T_INDEX &index, const T_INDEX &counts)
// {
//     int os = index[1] * counts[0] + index[0];
//     if (d == 3)
//         os += index[2] * counts[0] * counts[1];
//     return os;
// }