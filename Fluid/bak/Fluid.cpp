#include "Fluid.h"

using namespace Nova;

// ------------------ Fluid Quantity-------------------

FluidQuantity::FluidQuantity(Grid<T, d> &grid, int axis, int number_of_ghost_cells)
{
    std::cout << "Correct Constructor!" << std::endl;

    this->grid = &grid;
    this->number_of_ghost_cells = number_of_ghost_cells;
    this->axis = axis;

    storing_counts = T_INDEX(grid.counts);
    if (axis != -1)
        storing_counts(axis) += 1;

    Phi = new double[storing_counts.Product()];
    Phi_new = new double[storing_counts.Product()];
}

FluidQuantity::~FluidQuantity()
{
    delete[] Phi;
    delete[] Phi_new;
}

void FluidQuantity::fill(double content)
{
    for (int i = 0; i < storing_counts.Product(); i++)
    {
        Phi[i] = content;
    }
}

TV FluidQuantity::computeVelocity(const T_INDEX &index, FluidQuantity *velocityField[d])
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

TV FluidQuantity::Clamp_To_Domain(const TV &location)
{
    TV tmp_location = TV(location);

    // Clamp to domain if not in domain
    T_INDEX tmp_index = (*grid).Cell(tmp_location, number_of_ghost_cells);

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
double FluidQuantity::linter(TV &location)
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
    double px00 = linter(at(c000), at(c100), offset[0]);
    double px10 = linter(at(c010), at(c110), offset[0]);
    double py0 = linter(px00, px10, offset[1]);

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

    double px01 = linter(at(c001), at(c101), offset[0]);
    double px11 = linter(at(c011), at(c111), offset[0]);
    double py1 = linter(px01, px11, offset[1]);

    double pz = linter(py0, py1, offset[2]);

    return pz;
}

void FluidQuantity::advect(const T_INDEX &index, double timestep, FluidQuantity *velocityField[d])
{
    TV velocity = computeVelocity(index, velocityField);

    TV location = (axis == -1) ? (*grid).Center(index) : (*grid).Face(axis, index);

    TV location_traceback = Clamp_To_Domain(location - timestep * velocity);

    new_at(index) = linter(location_traceback);
}

void FluidQuantity::advect(double timestep, FluidQuantity *velocityField[d])
{
    T_INDEX currIndex;
    for (Range_Iterator<d> iterator(Range<int, d>(T_INDEX(1), storing_counts)); iterator.Valid(); iterator.Next())
    {
        currIndex = T_INDEX() + iterator.Index();
        advect(currIndex, timestep, velocityField);
    }
}

// ---------------- Fluid Solver----------------------
//
//
//
// ---------------------------------------------------

double FluidSolver::getRGBcolorDensity(T_INDEX &index)
{
    return (*density_field).rgb_at(index);
}

// TODO
void FluidSolver::initialize()
{

    // for (int i = 0; i < d; i++)
    // {
    //     (*velocityField[i]).fill(0.5);
    // }
    (*velocityField[0]).fill(0);
    (*velocityField[1]).fill(0);
    // (*velocityField[2]).fill(0.5);

    (*density_field).fill(0);

}

void FluidSolver::advection(double timestep)
{
    (*density_field).advect(timestep, velocityField);

    for (int i = 0; i < d; i++)
    {
        (*velocityField[i]).advect(timestep, velocityField);
    }
}

// TODO
void copy_(double src[], double dst[], int size)
{
    for (int i = 0; i < size; i++)
    {
        dst[i] = src[i];
    }
}

// TODO
void printArray(double arr[], int size)
{
    std::cout << "----\n";
    for (int i = 0; i < size; i++)
    {

        if (i % 16 == 0)
        {
            std::cout << "\n";
        }

        std::cout << arr[i] << ", ";
    }
    std::cout << std::endl;
}

void FluidSolver::calculateDivergence()
{

    T_INDEX currIndex;
    for (Range_Iterator<d> iterator(Range<int, d>(T_INDEX(1), storing_counts)); iterator.Valid(); iterator.Next())
    {
        currIndex = T_INDEX() + iterator.Index();
        int idx = index2offset(currIndex);
        divG[idx] = 0;

        for (int axis = 0; axis < d; axis++)
        {
            T_INDEX n_index = Next_Cell(axis, currIndex);
            if ((*grid).Inside_Domain(n_index))
            {
                divG[idx] += ((*velocityField[axis]).at(n_index) - (*velocityField[axis]).at(currIndex)) * (*grid).one_over_dX[axis];
                // divG[idx] += ((*velocityField[axis]).new_at(n_index) - (*velocityField[axis]).new_at(currIndex)) * (*grid).one_over_dX[axis];
            }
        }
    }

}

double calculate1Norm(double a[], double b[], int size)
{
    double max = __DBL_MIN__;
    for (int i = 0; i < size; i++)
        max = std::max(std::abs(a[i] - b[i]), max);

    return max;
}

void FluidSolver::Project(int limit)
{
    double scale = (*grid).one_over_dX[0] * (*grid).one_over_dX[0];
    double maxDelta;

    for (int iteration = 0; iteration < limit; iteration++)
    {
        maxDelta = __DBL_MIN__;
        T_INDEX index;
        for (Range_Iterator<d> iterator(Range<int, d>(T_INDEX(1), (*grid).Number_Of_Cells())); iterator.Valid(); iterator.Next())
        {
            index = T_INDEX() + iterator.Index();

            double Aii = 0;
            double AijSum = 0;
            int os = index2offset(index);

            for (int axis = 0; axis < d; axis++)
            {
                T_INDEX p_index = Previous_Cell(axis, index);
                T_INDEX n_index = Next_Cell(axis, index);
                if ((*grid).Inside_Domain(p_index))
                {
                    Aii -= scale;
                    AijSum += scale * pressure_solution[index2offset(p_index)];
                }
                if ((*grid).Inside_Domain(n_index))
                {
                    Aii -= scale;
                    AijSum += scale * pressure_solution[index2offset(n_index)];
                }
            }

            double newP = (divG[os] - AijSum) / Aii;

            maxDelta = std::max(maxDelta, fabs(pressure_solution[os] - newP));

            pressure_solution[os] = newP;
        }

        if (maxDelta < 1e-5)
        {
            printf("Exiting solver after %d iterations, maximum change is %f\n", iteration, maxDelta);
            return;
        }
    }
    printf("Exceeded budget of %d iterations, maximum change was %f\n", limit, maxDelta);
}

void FluidSolver::updateVelocity(double timestep)
{
    // Update Velocity
    T_INDEX index;
    for (Range_Iterator<d> iterator(Range<int, d>(T_INDEX(1), (*grid).Number_Of_Cells())); iterator.Valid(); iterator.Next())
    {
        index = T_INDEX() + iterator.Index();
        // pressure_solution[index2offset(index)] = size * divG[index2offset(index)];

        for (int axis = 0; axis < d; axis++)
        {
            if ((*grid).Inside_Domain(Previous_Cell(axis, index)))
            {
                (*velocityField[axis]).new_at(index) -= timestep * (pressure_solution[index2offset(index)] - pressure_solution[index2offset(Previous_Cell(axis, index))]) *
                                                        (*grid).one_over_dX(axis);
            }
        }
    }
}

void FluidSolver::flip()
{
    (*density_field).flip();

    for (int i = 0; i < d; i++)
    {
        (*velocityField[i]).flip();
    }
}

void FluidSolver::addInflow(const T_INDEX &index, const double density, const TV &velocity)
{
    (*density_field).modify_at(index) = density;

    for (int i = 0; i < d; i++)
    {
        (*velocityField[i]).modify_at(index) = velocity[i];
    }
}

void FluidSolver::SetPressureBoundary()
{
    //1.m -> m,n
    Range<int, d> vTop(T_INDEX{1, (*grid).counts[0]}, (*grid).counts);
    // 1,1 -> 1, m
    Range<int, d> uLeft(T_INDEX(1), T_INDEX{1, (*grid).counts[0]});
    // n,1 -> m,n
    Range<int, d> uRight(T_INDEX{(*grid).counts[1], 1}, (*grid).counts);

    T_INDEX currIndex;
    for (Range_Iterator<d> iterator(vTop); iterator.Valid(); iterator.Next())
    {
        currIndex = T_INDEX() + iterator.Index();
        pressure_solution[index2offset(currIndex)] = pressure_solution[index2offset(Previous_Cell(1, currIndex))];
    }

    for (Range_Iterator<d> iterator(uLeft); iterator.Valid(); iterator.Next())
    {
        currIndex = T_INDEX() + iterator.Index();
        pressure_solution[index2offset(currIndex)] = pressure_solution[index2offset(Next_Cell(0, currIndex))];
    }

    for (Range_Iterator<d> iterator(uRight); iterator.Valid(); iterator.Next())
    {
        currIndex = T_INDEX() + iterator.Index();
        pressure_solution[index2offset(currIndex)] = pressure_solution[index2offset(Previous_Cell(0, currIndex))];
    }
}

void FluidSolver::SetVelocityBoundary()
{
    //1.m -> m,n
    Range<int, d> vTop(T_INDEX{1, (*grid).counts[0]}, (*grid).counts);
    // 1,1 -> 1, m
    Range<int, d> uLeft(T_INDEX(1), T_INDEX{1, (*grid).counts[0]});
    // n,1 -> m,n
    Range<int, d> uRight(T_INDEX{(*grid).counts[1], 1}, (*grid).counts);

    T_INDEX currIndex;
    for (Range_Iterator<d> iterator(vTop); iterator.Valid(); iterator.Next())
    {
        currIndex = T_INDEX() + iterator.Index();
        (*velocityField[0]).new_at(currIndex) = (*velocityField[0]).new_at(Previous_Cell(1, currIndex));
    }

    for (Range_Iterator<d> iterator(uLeft); iterator.Valid(); iterator.Next())
    {
        currIndex = T_INDEX() + iterator.Index();
        (*velocityField[1]).new_at(currIndex) = (*velocityField[1]).new_at(Next_Cell(0, currIndex));
    }

    for (Range_Iterator<d> iterator(uRight); iterator.Valid(); iterator.Next())
    {
        currIndex = T_INDEX() + iterator.Index();
        (*velocityField[1]).new_at(currIndex) = (*velocityField[1]).new_at(Previous_Cell(0, currIndex));
    }
}

void FluidSolver::SetVelocityBoundary_as0()
{
    //1.m -> m,n
    Range<int, d> vTop(T_INDEX{1, (*grid).counts[0]}, (*grid).counts);
    // 1,1 -> 1, m
    Range<int, d> uLeft(T_INDEX(1), T_INDEX{1, (*grid).counts[0]});
    // n,1 -> m,n
    Range<int, d> uRight(T_INDEX{(*grid).counts[1], 1}, (*grid).counts);

    T_INDEX currIndex;
    for (Range_Iterator<d> iterator(vTop); iterator.Valid(); iterator.Next())
    {
        currIndex = T_INDEX() + iterator.Index();
        (*velocityField[1]).new_at(currIndex) = 0;
    }

    for (Range_Iterator<d> iterator(uLeft); iterator.Valid(); iterator.Next())
    {
        currIndex = T_INDEX() + iterator.Index();
        (*velocityField[0]).new_at(currIndex) = 0;
    }

    for (Range_Iterator<d> iterator(uRight); iterator.Valid(); iterator.Next())
    {
        currIndex = T_INDEX() + iterator.Index();
        (*velocityField[0]).new_at(currIndex) = 0;
    }
}

void FluidSolver::SetDivBoundary()
{
    //1.m -> m,n
    Range<int, d> vTop(T_INDEX{1, (*grid).counts[0]}, (*grid).counts);
    // 1,1 -> 1, m
    Range<int, d> uLeft(T_INDEX(1), T_INDEX{1, (*grid).counts[0]});
    // n,1 -> m,n
    Range<int, d> uRight(T_INDEX{(*grid).counts[1], 1}, (*grid).counts);

    T_INDEX currIndex;
    for (Range_Iterator<d> iterator(vTop); iterator.Valid(); iterator.Next())
    {
        currIndex = T_INDEX() + iterator.Index();
        divG[index2offset(currIndex)] = divG[index2offset(Previous_Cell(1, currIndex))];
    }

    for (Range_Iterator<d> iterator(uLeft); iterator.Valid(); iterator.Next())
    {
        currIndex = T_INDEX() + iterator.Index();
        divG[index2offset(currIndex)] = divG[index2offset(Next_Cell(0, currIndex))];
    }

    for (Range_Iterator<d> iterator(uRight); iterator.Valid(); iterator.Next())
    {
        currIndex = T_INDEX() + iterator.Index();
        divG[index2offset(currIndex)] = divG[index2offset(Previous_Cell(0, currIndex))];
    }
}

void FluidSolver::update(double timestep)
{
    calculateDivergence();
    advection(timestep);

    Project(1000);
    updateVelocity(timestep);
    SetVelocityBoundary_as0();

    flip();

    for (int i = 0; i < size; i++)
    {
        pressure_solution[i] = 0;
    }

}
