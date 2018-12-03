#include "CG_Driver.h"

using namespace Nova;

enum
{
    d = 2
};
using T = double;
using TV = Vector<T, d>;
using T_INDEX = Vector<int, d>;

void printArray(Array<TV> array, int size)
{
    for (int i = 0; i < size; i++)
    {
        TV tmp = array(i);
        std::cout << tmp(0) << ", ";

        if ((i + 1) % 16 == 0)
        {
            std::cout << "\n";
        }
    }
    std::cout << std::endl;
}

void test1(CG_System<T, d> &cg_system, Array<TV> &result, int size);

int main(int argc, char const *argv[])
{

    int counts = 16;
    Range<T, d> range(TV(-0.5), TV(0.5));
    Grid<T, d> *grid = new Grid<T, d>(T_INDEX(counts), range);
    CG_Storage<T, d> storage(*grid);

    storage.setting();
    storage.initialize();
    storage.calculateA();
    storage.testRHS();

    // CG_Driver<T,d> driver(storage);
    // driver.Execute();
    CG_System<T, d> cg_system(storage);

    Array<TV> result(storage.size);
    // Will Print A
    test1(cg_system, result, storage.size);
    // storage.printAdiag();

    std::cout << "one_over_dX_square = " << storage.one_over_dX_square << std::endl;
    std::cout << "grid->one_over_dX = " << grid->one_over_dX << std::endl;

    return 0;
}

void unitVector(Array<TV> &e, int j, int size)
{
    e.resize(size);
    e.Fill(TV());
    e(j) = 1;
}

void test1(CG_System<T, d> &cg_system, Array<TV> &result, int size)
{
    CG_Vector<T, d> cg_result(result);
    for (int i = 0; i < size; i++)
    {

        for (int j = 0; j < size; j++)
        {
            int idx = j * 16 + i;
            Array<TV> e1, e2;
            unitVector(e1, i, size);
            unitVector(e2, j, size);

            CG_Vector<T, d> cg_e1(e1), cg_e2(e2);

            cg_system.Multiply(cg_e1, cg_result);
            double re = cg_system.Inner_Product(cg_e2, cg_result);
            // std::cout <<"re(" << i << ")= " << re;
            printf("%3.1f", re);
            if ((idx + 1) % 256 == 0)
            {
                std::cout << "\n";
            }
            else{
                std::cout <<",";
            }
        }
    }
    std::cout << std::endl;
}
