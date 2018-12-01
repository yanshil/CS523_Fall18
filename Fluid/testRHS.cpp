#include <iostream>
#include "FluidSolver.h"

using namespace Nova;
int main(int argc, char const *argv[])
{
    enum
    {
        d = 2
    };
    using T = double;
    using T_INDEX = Vector<int, d>;
    using TV = Vector<T, d>;

    FluidSimulator_Grid<T, d> grid(T_INDEX(32), Range<T, d>::Unit_Box());
    FluidSolver<T, d> *solver = new FluidSolver<T, d>(grid, 1);

    double timestep = 0.01;
    
    solver->addInflow(T_INDEX{15, 15}, T_INDEX{17, 17}, 0.0, TV{5, 0});
    solver->unitTest_tentProjection(0.01);
    return 0;
}
