#include <nova/Tools/Parsing/Parse_Args.h>
#include <nova/Tools/Utilities/File_Utilities.h>

#include <iostream>

#include "Fluid.h"
#include "lodepng.h"

using namespace Nova;

/*-------- Global typedef & varaible------------*/

using T_INDEX = Vector<int, d>;

int main(int argc, char *argv[])
{

    double timestep = 0.01;
    double density = 1;

    double time = 0.0;
    int iterations = 0;

    unsigned char *image = new unsigned char[32 * 32 * 4];

    Grid<T, d> grid(T_INDEX(16), Range<T, d>::Unit_Box());
    FluidSolver *solver = new FluidSolver(grid, 0, 1);

    solver->initialize();

    while (time < 8.0)
    {
        /* Use four substeps per iteration */
        for (int i = 0; i < 4; i++)
        {
            solver->addInflow(T_INDEX{8, 8}, density, TV{0, 0.2});
            solver->update(timestep);
            time += timestep;
            fflush(stdout);
        }

        solver->toImage(image);

        char path[256];
        sprintf(path, "Frame%03d.png", iterations++);
        lodepng_encode32_file(path, image, 32, 32);
    }

    return 0;
}
