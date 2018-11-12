#include "Fluid.h"

using namespace Nova;

T_INDEX Next_Cell(const int axis, const T_INDEX &index);
T_INDEX Previous_Cell(const int axis, const T_INDEX &index);

// ------------------------------------------------------

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

/* Calculate (b-a) / delta */
double slope(double b, double a, double delta)
{
    return (b - a) / delta;
}

// ------------------------------------------------------

// ------------------ Fluid Quantity-------------------

FluidQuantity::FluidQuantity(Grid<T, d> &grid, int axis)
{
    std::cout<<"Correct Constructor!"<<std::endl;

    Phi = new double[grid.counts.Product()];
    Phi_new = new double[grid.counts.Product()];
    this->grid = &grid;
    this->number_of_ghost_cells = 1;
    this->axis = axis;
}

FluidQuantity::~FluidQuantity()
{
    delete[] Phi;
    delete[] Phi_new;
}

TV FluidQuantity::computeVelocity(const T_INDEX &index, FluidQuantity *velocityField[d])
{
    TV velocity = TV(0.0);

    // If Scalar
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
}

TV FluidQuantity::Clamp_To_Domain(const TV &location)
{
    TV tmp_location = TV(location);

    // Clamp to domain if not in domain
    T_INDEX tmp_index = (*grid).Cell(tmp_location, number_of_ghost_cells);

    if (!(*grid).Inside_Domain(tmp_index, number_of_ghost_cells))
    {
        tmp_location = (TV)(*grid).domain.Clamp(tmp_location);
    }

    return tmp_location;
}

TV FluidQuantity::traceBack(const TV &location, double timestep, TV &velocity)
{
    return Clamp_To_Domain(location - timestep * velocity);
}

/*!
 * Linear Interpolator for TV{i, j, k} on grid
 * Coordinates will be clamped to lie in simulation domain
*/
double FluidQuantity::linter(TV &location)
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
    if (d == 3)
    {
        c001 = Next_Cell(2, c000);
        c101 = Next_Cell(2, c100);
        c011 = Next_Cell(2, c010);
        c111 = Next_Cell(2, c110);
    }

    double px00 = linter(at(c000), at(c100), offset[0]);
    double px10 = linter(at(c010), at(c110), offset[0]);
    double py0 = linter(px00, px10, offset[1]);

    if (d == 2)
        return py0;
    else
    {
        double px01 = linter(at(c001), at(c101), offset[0]);
        double px11 = linter(at(c011), at(c111), offset[0]);
        double py1 = linter(px00, px10, offset[1]);

        double pz = linter(py0, py1, offset[2]);

        return pz;
    }
}

void FluidQuantity::advect(const T_INDEX &index, double timestep, FluidQuantity *velocityField[d])
{
    TV location;

    if (axis != -1)
        location = (*grid).Center(index);
    else
        location = (*grid).Face(axis, index);

    TV velocity = computeVelocity(index, velocityField);

    TV location_traceback = traceBack(location, timestep, velocity);
    new_at(index) = linter(location_traceback);
}

// /** Do Projection for given cell
//  *
//  */
// void FluidQuantity::projection(double timestep, T_INDEX &index, FluidQuantity *velocityField[d])
// {
//     //[p(i+1, j) - p(i,j)] / dX
//     double np0, np1, np = 0;

//     for (int i = 0; i < d; i++)
//     {
//         np0 = slope(new_at(Next_Cell(i, index)), new_at(index), (*grid).dX(i));
//         np1 = slope(new_at(index), new_at(Previous_Cell(i, index)), (*grid).dX(i));
//         np += slope(np0, np1, (*grid).dX(i));
//     }

//     double tildep = timestep * np;

//     double rhs = 0;

//     for (int i = 0; i < d; i++)
//     {
//         // [u_{i+1/2, j} - u_{i-1/2, j}] / Delta x
//         rhs += slope((*velocityField[i]).new_at(Next_Cell(i, index)),
//                      (*velocityField[i]).new_at(index), (*grid).dX(i));
//     }
// }

// -------- Fluid Solver Implementation -----------------------

void FluidSolver::advection(double timestep)
{
    T_INDEX currIndex;

    // For Each Cell in Grid
    for (Range_Iterator<d> iterator(Range<int, d>(T_INDEX(1), (*grid).Number_Of_Cells())); iterator.Valid(); iterator.Next())
    {
        // For each Fluid Velocity
        for (int i = 0; i < d; i++)
        {
            // const T_INDEX &index, double timestep, FluidQuantity *velocityField[d]
            (*velocityField[i]).advect(currIndex, timestep, velocityField);
        }
    }
}

void FluidSolver::calculateRHS()
{
    T_INDEX currIndex;

    // Min Corner:(1,1,1) To Max Corner (T_INDEX counts)
    for (Range_Iterator<d> iterator(Range<int, d>(T_INDEX(1), (*grid).Number_Of_Cells())); iterator.Valid(); iterator.Next())
    {
        currIndex = T_INDEX() + iterator.Index();
        rhs[index2offset(currIndex)] = 0;

        for (int i = 0; i < d; i++)
        { // [u_{i+1/2, j} - u_{i-1/2, j}] / Delta x
            rhs[index2offset(currIndex)] += slope((*velocityField[i]).new_at(Next_Cell(i, currIndex)),
                                                  (*velocityField[i]).new_at(currIndex), (*grid).dX(i));
        }
    }
}

void FluidSolver::projection()
{
    
}

// TODO: Timestep get involve??
void FluidSolver::updateVelocity(T_INDEX &index, double timestep)
{
    for (int i = 0; i < d; i++)
    {
        double np;
        // p.at(index) - p.at(PreviousCell(i, index)) / DeltaX(i)
        np = slope(pressure_solution[index2offset(index)],
                   pressure_solution[index2offset(Previous_Cell(i, index))], (*grid).dX(i));

        (*velocityField[i]).new_at(index) -= np * timestep;
    }
}

void FluidSolver::updateVelocity(double timestep)
{
    T_INDEX currIndex;

    // Min Corner:(1,1,1) To Max Corner (T_INDEX counts)
    for (Range_Iterator<d> iterator(Range<int, d>(T_INDEX(1), (*grid).Number_Of_Cells())); iterator.Valid(); iterator.Next())
    {
        updateVelocity(currIndex, timestep);
    }
}

void FluidSolver::flip()
{
    (*density_field).flip();

    for (int i = 0; i < d; i++)
    {
        (*velocityField)[i].flip();
    }
}


void FluidSolver::update(double timestep)
{
    // Set rhs
    calculateRHS();

    advection(timestep);
    projection();
    updateVelocity(timestep);
    flip();

}
