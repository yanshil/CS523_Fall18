#include "CG_Driver.h"
#define CELLSIZE 4

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

        if ((i + 1) % CELLSIZE == 0)
        {
            std::cout << "\n";
        }
    }
    std::cout << std::endl;
}

void test1(CG_System<T, d> &cg_system, Array<TV> &result, int size);

int main(int argc, char const *argv[])
{
    Range<T, d> range(TV(-0.5), TV(0.5));
    Grid<T, d> *grid = new Grid<T, d>(T_INDEX(CELLSIZE), range);
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

    for (int j = 0; j < size; j++)
    {
        for (int i = 0; i < size; i++)
        {
            int idx = j * size + i;
            
            Array<TV> e1, e2;
            unitVector(e1, i, size);
            unitVector(e2, j, size);

            CG_Vector<T, d> cg_e1(e1), cg_e2(e2);

            cg_system.Multiply(cg_e1, cg_result);
            double re = cg_system.Inner_Product(cg_e2, cg_result);


            printf("%3.0f", re);

            if ((idx + 1) % (CELLSIZE * CELLSIZE) == 0)
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

void test2()
{
    // Random x and y

    // x^T A y = y^T A x
    // if A is symmetric, compute x^T A y an y^T A x and check if equal
}
