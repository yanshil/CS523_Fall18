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
            T_INDEX next_cell_index(index);
            T_INDEX next_cell_index2(Next_Cell(i, next_cell_index));

            velocity[i] += 0.5 * FluidVelocity[i]->at(index);
            velocity[i] += 0.5 * FluidVelocity[i]->at(next_cell_index2);
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

TV traceBack(TV &location, float timestep, TV &velocity)
{
    TV tmp_location = TV(location);

    return tmp_location -= timestep * velocity;
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

/*!
 * Auxiliary Function: Get exact offset for Column-based 1D array 
 * return (z * xSize * ySize) + (y * xSize) + x;
 */
int index2offset(const T_INDEX &index, const T_INDEX &counts)
{
    int os = index[1] * counts[0] + index[0];
    if (d == 3)
        os += index[2] * counts[0] * counts[1];
    return os;
}