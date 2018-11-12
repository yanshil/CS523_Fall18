#include <nova/Tools/Grids/Grid.h>
#include <nova/Tools/Parsing/Parse_Args.h>
#include <nova/Tools/Utilities/File_Utilities.h>

using namespace Nova;

enum
{
    d = 2
};

typedef float T;
typedef Vector<T, d> TV;
typedef Vector<int, d> T_INDEX;

Grid<T, d> *grid;
int number_of_ghost_cells = 0;

main(int argc, char *argv[])
{
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
    

    return 0;
}
