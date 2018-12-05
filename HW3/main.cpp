#include "CG_Driver.h"
#define CELLSIZE 16

// 16×16 ,32×32, 64×64, 128×128 and 256×256
// 25, 51, 102, 207, 420

using namespace Nova;

int main(int argc, char const *argv[])
{

    enum
    {
        d = 2
    };
    using T = double;
    using TV = Vector<T, d>;
    using T_INDEX = Vector<int, d>;

    Range<T,d> range(TV(-0.5), TV(0.5));
    Grid<T,d> * grid = new Grid<T,d>(T_INDEX(CELLSIZE), range);
    CG_Storage<T,d> storage(*grid);

    // void setting(tolerance, iterations, restart_iteration)
    storage.setting(1e-05, 1000, 1000);
    storage.initialize();
    storage.calculateA();
    storage.calculateRHS();
    
    CG_Driver<T,d> driver(storage);

    driver.Execute();

    return 0;
}
