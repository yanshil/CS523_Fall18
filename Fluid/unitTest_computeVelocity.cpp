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

    Grid<T, d> grid(T_INDEX(8), Range<T, d>::Unit_Box());
    FluidSolver *solver = new FluidSolver(grid, 0, 1);

    int iteration = 1;

    solver->initialize();

    solver->addInflow(T_INDEX{3, 1}, density, TV{0, 1});
    solver->update(0.3);
    std::cout << "--------------------------" << std::endl;
    solver->update(0.3);
    

    return 0;
}
