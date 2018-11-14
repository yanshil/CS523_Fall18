#include <nova/Tools/Grids/Grid.h>
#include "Fluid.h"

using namespace Nova;

typedef double T;
typedef Vector<T, d> TV;
typedef Vector<int, d> T_INDEX;

int num_of_ghosts = 1;

Grid<T, d> grid(T_INDEX(8), Range<T, d>::Unit_Box());

void linterTest(TV &location, FluidQuantity &fq)
{
    std::cout << "axis = " << fq.axis << "\t";
    double re = fq.linter(location);
    // std::cout << "output = " << re << std::endl;
    std::cout << "\n-----------" << std::endl;
}

TV printlocationInfo(const TV &location)
{
    std::cout << "-------------Start------------ " << std::endl;
    std::cout << "location = " << location << "\t";
    T_INDEX cindex = grid.Clamp_To_Cell(location, num_of_ghosts);
    std::cout << "cindex = " << cindex << std::endl;
    std::cout << "---------------- " << std::endl;

    return location;
}

main(int argc, char *argv[])
{

    FluidQuantity *u = new FluidQuantity(grid, 0, num_of_ghosts);
    FluidQuantity *v = new FluidQuantity(grid, 1, num_of_ghosts);
    // FluidQuantity *w = new FluidQuantity(grid, 2, num_of_ghosts);
    FluidQuantity *d = new FluidQuantity(grid, -1, num_of_ghosts);

    std::cout << grid.dX << std::endl;

    TV location = printlocationInfo(TV(0.21));

    linterTest(location, *d);
    linterTest(location, *u);
    linterTest(location, *v);
    // linterTest(location, *w);

    location = printlocationInfo(TV(0.77));

    linterTest(location, *d);
    linterTest(location, *u);
    linterTest(location, *v);
    // linterTest(location, *w);

    location = printlocationInfo(TV(0.4));

    linterTest(location, *d);
    linterTest(location, *u);
    linterTest(location, *v);
    // linterTest(location, *w);

    location = printlocationInfo(TV(0.0));

    linterTest(location, *d);
    linterTest(location, *u);
    linterTest(location, *v);
    // linterTest(location, *w);

    return 0;
}
