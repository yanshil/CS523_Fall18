#include "FluidQuantity.h"

using namespace Nova;

enum
{
    d = 2
};
using T = double;
using T_INDEX = Vector<int, d>;
using TV = Vector<T, d>;

int number_of_ghost = 1;
FluidSimulator_Grid<T, d> grid(T_INDEX(32), Range<T, d>::Unit_Box(), number_of_ghost);

/*-------- Global typedef & varaible------------*/

using T_INDEX = Vector<int, d>;

void testIndex2Offset(const T_INDEX &index, FluidQuantity<T,d>& u)
{
    std::cout<<"Original index = "<< index << std::endl;
    int os = grid.index2offset(index);
    T_INDEX index_new = grid.offset2index(os);
    std::cout<<"Middle offset = "<< os << std::endl;
    std::cout<<"New index = " << index_new << std::endl;
    std::cout<<"------"<<std::endl;

}

void testOffset2Index(const int os, FluidQuantity<T,d>& u)
{
    std::cout<<"Original offset = "<< os << std::endl;
    T_INDEX index = grid.offset2index(os);
    int os_new = grid.index2offset(index);
    
    std::cout<<"Middle Index = "<< index << std::endl;
    std::cout<<"New offset = " << os_new << std::endl;
    std::cout<<"------"<<std::endl;

}

main(int argc, char const *argv[])
{
    

    FluidQuantity<T,d> *u = new FluidQuantity<T,d>(grid, 1, number_of_ghost);

    std::cout<<"================="<<std::endl;
    testIndex2Offset(T_INDEX{1,1}, *u);
    testIndex2Offset(T_INDEX{4,5}, *u);
    testIndex2Offset(T_INDEX{2,2}, *u);
    testIndex2Offset(T_INDEX{5,4}, *u);
    testIndex2Offset(T_INDEX{16,16}, *u);

    std::cout<<"================="<<std::endl;
    
    testOffset2Index(100, *u);
    testOffset2Index(1, *u);
    testOffset2Index(35, *u);
    testOffset2Index(256, *u);
    testOffset2Index(255, *u);

    // testIndex2Offset(T_INDEX{2,2}, *u);
    // testIndex2Offset(T_INDEX{2,1}, *u);
    // testIndex2Offset(T_INDEX{1,2}, *u);
    // testIndex2Offset(T_INDEX{2,3}, *u);
    // testIndex2Offset(T_INDEX{3,2}, *u);
    testIndex2Offset(T_INDEX{32,32}, *u);
    testIndex2Offset(T_INDEX{1,32}, *u);
    testIndex2Offset(T_INDEX{0,32}, *u);
    testIndex2Offset(T_INDEX{0,33}, *u);
    testIndex2Offset(T_INDEX{33,33}, *u);

    delete u;
    
    return 0;
}
