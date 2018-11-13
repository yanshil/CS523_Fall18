#include <nova/Tools/Parsing/Parse_Args.h>
#include <nova/Tools/Utilities/File_Utilities.h>
#include <iostream>

#include "Fluid.h"

using namespace Nova;

/*-------- Global typedef & varaible------------*/

using T_INDEX = Vector<int, d>;

// -----------Parse Argument -------------------


T_INDEX parseArgument(int argc, char *argv[])
{
    /* Customize Grid Size with Parsing Arguments */
    Parse_Args parse_args;
    if (d == 2)
        parse_args.Add_Vector_2D_Argument("-size", Vector<double, 2>(16), "", "Grid resolution");
    else
        parse_args.Add_Vector_3D_Argument("-size", Vector<double, 3>(16), "", "Grid resolution");
    parse_args.Add_String_Argument("-o", ".", "", "Output directory");
    parse_args.Parse(argc, argv);

    T_INDEX counts;
    if (d == 2)
    {
        auto counts_2d = parse_args.Get_Vector_2D_Value("-size");
        for (int v = 0; v < d; ++v)
            counts(v) = counts_2d(v);
    }
    else
    {
        auto counts_3d = parse_args.Get_Vector_3D_Value("-size");
        for (int v = 0; v < d; ++v)
            counts(v) = counts_3d(v);
    }
    std::string output_directory = parse_args.Get_String_Value("-o");
    
    // Log::cout << "Counts: " << counts << std::endl;

    return counts;
}

// -------------------------------------------------

int main(int argc, char **argv)
{
    T_INDEX counts = parseArgument(argc, argv);

    Grid<T, d> grid(counts, Range<T, d>::Unit_Box());

    double timestep = 0.12;
    double density = 0.5;


    FluidSolver *solver = new FluidSolver(grid, density);

    double time = 0.0;
    int iterations = 0;

    solver->initialize();

    while(time < 8)
    {
        // addInflow(T_INDEX &index, double density, TV &velocity);
        solver->addInflow(T_INDEX(1), density, TV{0,0.2});
        solver->update(timestep);
        time += timestep;
    }
    
    return 0;
}