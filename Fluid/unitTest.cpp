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

void test_tent();
void test_grid();

int main(int argc, char const *argv[])
{
    // test_tent();
    test_grid();
    return 0;
}
//----------------------------------

void test_grid()
{
    Range<T, d> range(TV(), TV(128));
    FluidSimulator_Grid<T, d> grid(T_INDEX(CELLCOUNTS), range);
    FluidSolver<T, d> *solver = new FluidSolver<T, d>(grid, 1);

    std::cout << "grid.dx= " << grid.dX << std::endl;

    // Projection output as a tent-like plot
    T timestep = 0.01;
    int size = grid.counts.Product();
    // T_INDEX index = T_INDEX{60, 60, 60};
    // z * m * n + y * m + x

    T *ptr_rhs = solver->rhs();
    memset(ptr_rhs, 0, size * sizeof(T));
    // solver->rhs()[59 * 128 * 128 + 59 * 128 + 59] = 100;
    solver->rhs()[7740] = 1;
    solver->output = false;
    solver->project_GS(1000, timestep);
    solver->applyPressure(timestep);

    // // Print Pressure Solution
    // for (int i = 0; i < size; i++)
    // {
    //     printf("%f", solver->p()[i]);

    //     int os = (i + 1) % grid.counts[0];

    //     if (os == 0)
    //     {
    //         printf("\n");
    //     }
    //     else
    //     {
    //         printf(",");
    //     }
    // }
}

void test_tent()
{
    FluidSimulator_Grid<T, d> grid(T_INDEX(CELLCOUNTS), Range<T, d>::Unit_Box());
    FluidSolver<T, d> *solver = new FluidSolver<T, d>(grid, 1);

    // Projection output as a tent-like plot
    T timestep = 0.01;
    int size = grid.counts.Product();
    T_INDEX index = T_INDEX{60, 60};
    // z * m * n + y * m + x
    T *ptr_rhs = solver->rhs();
    memset(ptr_rhs, 0, size * sizeof(T));
    // solver->rhs()[59 * 128 * 128 + 59 * 128 + 59] = 100;
    solver->rhs()[7740] = 1;
    solver->output = false;
    solver->project_GS(1000, timestep);
    solver->applyPressure(timestep);

    // Print Pressure Solution
    for (int i = 0; i < size; i++)
    {
        printf("%f", solver->p()[i]);

        int os = (i + 1) % grid.counts[0];

        if (os == 0)
        {
            printf("\n");
        }
        else
        {
            printf(",");
        }
    }
}
