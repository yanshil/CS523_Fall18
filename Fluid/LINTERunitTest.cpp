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

/*!
     * 1D Linear Interpolate etween a and b for x in (0, 1)
     */
float linter(float a, float b, float x)
{
    return (1.0 - x) * a + x * b;
}

/*!
     * Linear Interpolator for TV{i, j, k} on grid
     * Coordinates will be clamped to lie in simulation domain
     */
float linter(TV &location)
{
    int axis = -1;
    T_INDEX clamp_index;
    clamp_index = (*grid).Clamp_To_Cell(location, number_of_ghost_cells);

    TV offset = location - (TV)clamp_index;

    float _o[d];
    
    for(size_t i = 0; i < d; i++)
    {

        if (axis == i) {
            _o[i] = 0.5;
        }
        else
        {
            _o[i] = 0;
        }

        
        if (axis == -1) {
            _o[i] = 0.5;
        }
         
    }
    

    float c00, c10, c01, c11;

    // this->at(clamp_index);
}

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
    // File_Utilities::Write_To_File(output_directory+"/grid.grid",grid);



    return 0;
}
