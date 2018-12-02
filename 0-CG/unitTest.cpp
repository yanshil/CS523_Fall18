#include <iostream>
#include "Fluid.h"

FluidSolver *solver = new FluidSolver(128, 128, 1);

void unitTest_Projection_Tent(double timestep)
{
    solver->addInflow(0.3, 0.3, 0.01, 0.01, 0, 3, 0);
    solver->test_projection_tent(timestep);
}

int main(int argc, char const *argv[])
{
    // Output as a csv file and plot
    unitTest_Projection_Tent(0.01);
    return 0;
}
