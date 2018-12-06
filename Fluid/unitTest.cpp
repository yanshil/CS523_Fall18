#include <GL/glut.h>
#include <iostream>
#include "FluidSolver.h"
#define CELLCOUNTS 128

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


    // Projection output as a tent-like plot
    T timestep = 0.01;
    int size = grid.counts.Product();

    T* ptr_rhs = solver->rhs();
    memset(ptr_rhs, 0, size * sizeof(T));
    solver->rhs()[7740] = 1;
    solver->project_GS(1000, timestep, false);
    solver->applyPressure(timestep);
    
    // Print Pressure Solution
    for(int i = 0; i < size; i++)
    {
        printf("%f", solver->p()[i]);

        int os = (i + 1) % grid.counts[0];

        if (os == 0) {
            printf("\n");
        }
        else
        {
            printf(",");
        }
    }
    return 0;
}



