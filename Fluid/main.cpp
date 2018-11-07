#include <nova/Tools/Parsing/Parse_Args.h>
#include <nova/Tools/Utilities/File_Utilities.h>
#include <iostream>

#include "FluidQuantity.h"

using namespace Nova;

/*-------- Global typedef & varaible------------*/

using T_INDEX = Vector<int, d>;

int main(int argc, char **argv)
{
    enum
    {
        d = 2
    };

    /* Customize Grid Size with Parsing Arguments */
    Parse_Args parse_args;
    if (d == 2)
        parse_args.Add_Vector_2D_Argument("-size", Vector<double, 2>(500), "", "Grid resolution");
    else
        parse_args.Add_Vector_3D_Argument("-size", Vector<double, 3>(500), "", "Grid resolution");
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

    Grid<T, d> grid(counts, Range<T, d>::Unit_Box());
    // File_Utilities::Write_To_File(output_directory+"/grid.grid",grid);

    return 0;
}