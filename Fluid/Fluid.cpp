#include "Fluid.h"

using namespace Nova;

/* Calculate (b-a) / delta */
double slope(double b, double a, double delta)
{
    return (b - a) / delta;
}

// ------------------ Fluid Quantity-------------------

FluidQuantity::FluidQuantity(Grid<T, d> &grid, int axis, int number_of_ghost_cells)
{
    std::cout << "Correct Constructor!" << std::endl;

    int size = grid.counts.Product();

    Phi = new double[size];
    Phi_new = new double[size];
    this->grid = &grid;
    this->number_of_ghost_cells = number_of_ghost_cells;
    this->axis = axis;
}

FluidQuantity::~FluidQuantity()
{
    delete[] Phi;
    delete[] Phi_new;
}

void FluidQuantity::fill(double content)
{
    for (int i = 0; i < (*grid).counts.Product(); i++)
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

    if (!(*grid).Inside_Domain(tmp_index, number_of_ghost_cells))
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

    // T_INDEX locationTB = (*grid).Clamp_To_Cell(location_traceback, number_of_ghost_cells);

    // int tmp2 = (index - locationTB).Sum();
    // double tmp = linter(location_traceback);

    // if (axis == -1)
    // {
    //     std::cout << "index = " << index << std::endl;
    //     std::cout << "axis = " << axis << std::endl;
    //     std::cout << "location = " << location << std::endl;
    //     std::cout << "velocity = " << velocity << std::endl;
    //     std::cout << "timestep = " << timestep << std::endl;
    //     std::cout << "location_traceback = " << location_traceback << std::endl;
    //     std::cout << "phi_new at = " << tmp << std::endl;
    //     std::cout << "---------------- " << std::endl;
    // }
}

// -----------------------

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
    (*velocityField[1]).fill(0.1);
    // (*velocityField[2]).fill(0.5);

    (*density_field).fill(0);
}

void FluidSolver::advection(double timestep)
{
    T_INDEX currIndex;
    // Advect Density
    for (Range_Iterator<d> iterator(Range<int, d>(T_INDEX(1), (*grid).Number_Of_Cells())); iterator.Valid(); iterator.Next())
    {
        currIndex = T_INDEX() + iterator.Index();
        (*density_field).advect(currIndex, timestep, velocityField);
    }

    // (*density_field).printPhi_new();
    // std::cout<<"advect density Finished..." << std::endl;

    // For each Fluid Velocity
    for (int i = 0; i < d; i++)
    {
        // For Each Cell in Grid
        for (Range_Iterator<d> iterator(Range<int, d>(T_INDEX(1), (*grid).Number_Of_Cells())); iterator.Valid(); iterator.Next())
        {
            currIndex = T_INDEX() + iterator.Index();
            (*velocityField[i]).advect(currIndex, timestep, velocityField);
        }
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
        std::cout << arr[i] << ", ";
    }
    std::cout << std::endl;
}

void FluidSolver::calculateRHS()
{

    // // Initialize
    // for (int i = 0; i < size; i++)
    //     rhs[i] = 0;

    for (int i = 0; i < size; i++)
    {
        rhs[i] = 0;
        T_INDEX index = offset2index(i);

        // double one_over_dXSquare = (*grid).one_over_dX[0] * (*grid).one_over_dX[0];
        // for(int axis = 0; axis < d; axis++)
        // {
        //     // Usually should be ONE if NEXT CELL exist
        //     if ((Next_Cell(axis, index) - index).Sum() != 0) {
        //         rhs[i] += one_over_dXSquare * pressure_solution[index2offset(Next_Cell(axis, index))];
        //         rhs[i] -= one_over_dXSquare * pressure_solution[i];
        //     }
        //     if ((Previous_Cell(axis, index) - index).Sum() != 0) {
        //         rhs[i] += one_over_dXSquare * pressure_solution[index2offset(Previous_Cell(axis, index))];
        //         rhs[i] -= one_over_dXSquare * pressure_solution[i];
        //     }
        // }

        // [u_{i+1/2, j, k} - u_{i-1/2, j, k}] / Delta x
        for (int axis = 0; axis < d; axis++)
        {
            rhs[i] += slope((*velocityField[axis]).at(Next_Cell(axis, index)),
                            (*velocityField[axis]).at(index), (*grid).dX[axis]);
        }
    }
}

void FluidSolver::pressure_solution_Jacobi()
{
}

double calculate1Norm(double a[], double b[], int size)
{
    double max = __DBL_MIN__;
    for (int i = 0; i < size; i++)
    {
        max = std::max(std::abs(a[i] - b[i]), max);
    }

    return max;
}

void FluidSolver::calculateNeighborI(int os, Vector<int, (2 * d)> &neighborI)
{
    // int os_Aijne[2 * d];
    T_INDEX index = offset2index(os);

    neighborI[0] = index2offset(Previous_Cell(0, index));
    neighborI[1] = index2offset(Next_Cell(0, index));
    neighborI[2] = index2offset(Previous_Cell(1, index));
    neighborI[3] = index2offset(Next_Cell(1, index));

    if (d == 3)
    {
        neighborI[4] = index2offset(Previous_Cell(2, index));
        neighborI[5] = index2offset(Next_Cell(2, index));
    }
}

double *FluidSolver::getSparseA(int offset)
{
    double spareA_i[(2 * d + 1)];

    T_INDEX index = offset2index(offset);
    double one_over_dXSquare = (*grid).one_over_dX[0] * (*grid).one_over_dX[0];

    for (int i = 0; i < (2 * d + 1); i++)
    {
        for (int axis = 0; axis < d; axis++)
        {
            // Lower
            if ((Previous_Cell(axis, index) - index).Sum() != 0)
            {
                int tmp = axis;
                spareA_i[tmp] += one_over_dXSquare * pressure_solution[index2offset(Previous_Cell(axis, index))];
                spareA_i[(2 * d)] -= one_over_dXSquare * pressure_solution[offset];
            }
            // Upper
            // Usually should be ONE if NEXT CELL exist
            if ((Next_Cell(axis, index) - index).Sum() != 0)
            {
                int tmp = axis + d;
                spareA_i[tmp] += one_over_dXSquare * pressure_solution[index2offset(Next_Cell(axis, index))];
                spareA_i[(2 * d)] -= one_over_dXSquare * pressure_solution[offset];
            }
        }
    }

    // for(int axis = 0; axis < d; axis++)
    // {
    //     // Usually should be ONE if NEXT CELL exist
    //     if ((Next_Cell(axis, index) - index).Sum() != 0) {
    //         rhs[i] += one_over_dXSquare * pressure_solution[index2offset(Next_Cell(axis, index))];
    //         rhs[i] -= one_over_dXSquare * pressure_solution[i];
    //     }
    //     if ((Previous_Cell(axis, index) - index).Sum() != 0) {
    //         rhs[i] += one_over_dXSquare * pressure_solution[index2offset(Previous_Cell(axis, index))];
    //         rhs[i] -= one_over_dXSquare * pressure_solution[i];
    //     }
    // }

    return spareA_i;
}

// TODO
void FluidSolver::projection()
{
    double Aij = (*grid).one_over_dX[0] * (*grid).one_over_dX[0];
    double Aii = (-4) * Aij;

    double x[size];
    double x_new[size];
    Vector<int, (2 * d)> neighborI;

    int iteration = 0;
    for (int i = 0; i < size; i++)
    {
        x[i] = 0;
        x_new[i] = 0;
    }

    while (true)
    {
        // // Jacobi
        // for (int i = 0; i < size; i++)
        // {
        //     x_new[i] = rhs[i] / Aii;

        //     int sum = 0;
        //     calculateNeighborI(i, neighborI);

        //     for (int j = 0; j < (2 * d); j++)
        //     {
        //         sum += x[neighborI[j]] * Aij;
        //     }

        //     x_new[i] -= sum / Aii;
        // }

        // Gauss-Seidel
        for (int i = 0; i < size; i++)
        {
            x_new[i] = rhs[i] / Aii;
            int sum = 0;
            calculateNeighborI(i, neighborI);

            for (int j = 0; j < (2 * d); j++)
            {
                if (neighborI[j] < i)
                    sum += Aij * x_new[neighborI[j]];
                if (neighborI[j] > i)
                    sum += Aij * x[neighborI[j]];
            }
            x_new[i] -= sum / Aii;
        }

        double norm1 = calculate1Norm(x, x_new, size);
        if (norm1 < 0.000001)
        {
            std::cout << "Natual Break with 1Norm = " << norm1 << std::endl;
            std::cout << "Iteration: " << iteration << std::endl;
            printArray(x_new, size);
            copy_(x_new, pressure_solution, size);
            break;
        }

        if (iteration >= 1000)
        {
            std::cout << "break because of large iteration." << std::endl;
            copy_(x_new, pressure_solution, size);
            printArray(x_new, size);
            break;
        }

        copy_(x_new, x, size);

        iteration++;
    }
}

// TODO: Timestep get involve??
void FluidSolver::updateVelocity(T_INDEX &index, double timestep)
{
    for (int i = 0; i < d; i++)
    {
        double np = nablapOnI(index, i);
        (*velocityField[i]).new_at(index) -= np * timestep;
    }
}

void FluidSolver::updateVelocity(double timestep)
{

    for (int i = 0; i < size; i++)
    {
        T_INDEX index = offset2index(i);
        updateVelocity(index, timestep);
    }

    // T_INDEX currIndex;

    // // Min Corner:(1,1,1) To Max Corner (T_INDEX counts)
    // for (Range_Iterator<d> iterator(Range<int, d>(T_INDEX(1), (*grid).Number_Of_Cells())); iterator.Valid(); iterator.Next())
    // {
    //     currIndex = T_INDEX() + iterator.Index();
    //     updateVelocity(currIndex, timestep);
    // }
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
    (*density_field).at(index) = density;

    for (int i = 0; i < d; i++)
    {
        (*velocityField[i]).at(index) = velocity[i];
    }
}

void FluidSolver::update(double timestep)
{
    // Set rhs
    calculateRHS();
    printArray(rhs, size);

    advection(timestep);

    projection();
    updateVelocity(timestep);

    flip();

    // std::cout << "----- END UPDATE ---" << std::endl;
}
