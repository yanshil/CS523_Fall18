
#include "CG_Driver.h"

using namespace Nova;

int main(int argc, char const *argv[])
{
    enum{d=2};
    using T = double;
    using TV = Vector<T,d>;
    using T_INDEX = Vector<int,d>;

    Range<T,d> range(TV(-0.5), TV(0.5));

    // (const T_INDEX &counts_input, const Range<T, d> &domain_input)
    // TODO: Parse Arg

    int counts = 16;
    Grid<T,d> * grid = new Grid<T,d>(T_INDEX(counts), range);
    
    CG_Driver<T,d> driver(*grid);

    driver.Execute();
    
    return 0;
}
