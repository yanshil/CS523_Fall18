#include "CG_Driver.h"
#include <nova/Tools/Parsing/Parse_Args.h>

// 16×16 ,32×32, 64×64, 128×128 and 256×256
// Error:
// Iterations:  25, 51, 102, 207, 420

using namespace Nova;

int main(int argc, char **argv)
{
    enum
    {
        d = 2
    };
    using T = double;
    using TV = Vector<T, d>;
    using T_INDEX = Vector<int, d>;

    Parse_Args parse_args;
    parse_args.Add_Vector_2D_Argument("-size", Vector<double, 2>(16), "", "Grid resolution");
    parse_args.Parse(argc, argv);

    T_INDEX counts;
    auto counts_2d = parse_args.Get_Vector_2D_Value("-size");
    for (int v = 0; v < d; ++v)
        counts(v) = counts_2d(v);

    Range<T, d> range(TV(-0.5), TV(0.5));    
    Grid<T, d> *grid = new Grid<T,d>(counts, range);
    CG_Storage<T, d> storage(*grid);

    // void setting(tolerance, iterations, restart_iteration)
    storage.setting(1e-05, 1000, 1000);
    storage.initialize();
    storage.calculateA();
    storage.calculateRHS();

    CG_Driver<T, d> driver(storage);

    driver.Execute();

    //Print Info
    Log::cout<<"Gird Counts: "<<counts<<std::endl;

    return 0;
}
