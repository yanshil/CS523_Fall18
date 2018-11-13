#include "Fluid.h"

using namespace Nova;

/*-------- Global typedef & varaible------------*/

using T_INDEX = Vector<int, d>;

//----------------------------------------------------

main(int argc, char *argv[])
{
    double timestep = 0.12;
    double density = 1;

    int iterations = 0;

    Grid<T, d> grid(T_INDEX(16), Range<T, d>::Unit_Box());
    // FluidSolver *solver = new FluidSolver(grid, 0, 1);

    int iteration = 1;

    // solver->initialize();

    FluidQuantity *velocityField[d];

    for (int i = 0; i < d; i++)
    {
        velocityField[i] = new FluidQuantity(grid, i, 1);
        (*velocityField[i]).fill(0.5);
    }

    // Scalar
    FluidQuantity *density_field = new FluidQuantity(grid, -1, 1);
    (*density_field).fill(0);

    (*velocityField[0]).printPhi();

    T_INDEX index = T_INDEX(1);

    //Test Compute Velocity
    TV velocity = (*velocityField[0]).computeVelocity(index, velocityField);

    std::cout << "velocity[0] at (" << index << ") = " << (*velocityField[0]).at(index) << std::endl;
    std::cout << "velocity[1] at (" << index << ") = " << (*velocityField[1]).at(index) << std::endl;
    std::cout << "velocity at (" << index << ") = " << velocity << std::endl;

    return 0;
}
