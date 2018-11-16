#include <nova/Tools/Grids/Grid.h>
#include <nova/Tools/Utilities/Range_Iterator.h>

using namespace Nova;

enum
{
    d = 2
};

typedef double T;
typedef Vector<T, d> TV;
typedef Vector<int, d> T_INDEX;


int main(int argc, char const *argv[])
{
    Grid<T, d> grid(T_INDEX(3), Range<T, d>::Unit_Box());

    
    return 0;
}
