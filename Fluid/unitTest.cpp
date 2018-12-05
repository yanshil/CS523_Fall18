#include <GL/glut.h>
#include <iostream>
#include "FluidSolver.h"
#define CELLCOUNTS 16

using namespace Nova;

enum
{
    d = 2
};
using T = double;
using TV = Vector<T, d>;
using T_INDEX = Vector<int, d>;

int main(int argc, char const *argv[])
{
    FluidSimulator_Grid<T, d> grid(T_INDEX(CELLCOUNTS), Range<T, d>::Unit_Box());
    FluidSolver<T, d> *solver = new FluidSolver<T, d>(grid, 1);

    solver->unitTest_tentProjection(0.01);
    return 0;
}
