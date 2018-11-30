#include <iostream>
#include "Fluid.h"

int main(int argc, char const *argv[])
{
    FluidSolver *solver = new FluidSolver(128, 128, 1);
    solver->addInflow(0.5, 0.5, 0.51, 0.51, 0, 10, 0);
    solver->test_projection_tent(0.01);
    return 0;
}
